import csv
import importlib
import importlib.util
import io
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path


class TestKojakImport(unittest.TestCase):
    fixture = (
        Path(__file__).resolve().parents[1]
        / "fixtures/kojak_pxd015044"
    )

    def setUp(self):
        self.expected = json.loads(
            (self.fixture / "expected.json").read_text(encoding="utf-8")
        )
        self.results = self.fixture / self.expected["result_filename"]
        self.fasta = self.fixture / "search.fasta"

    def api(self):
        name = "claudio.dataset_catalog.kojak"
        self.assertIsNotNone(
            importlib.util.find_spec(name),
            "Missing Kojak adapter",
        )
        module = importlib.import_module(name)
        self.assertTrue(callable(module.read_kojak_csms))
        return module

    def read_bundle(self, results=None, fasta=None):
        return self.api().read_kojak_csms(
            results or self.results,
            fasta or self.fasta,
            dataset_accession="PXD015044",
            decoy_prefix="decoy_",
        )

    def paired(self, bundle):
        return [
            record for record in bundle["records"]
            if record["record_kind"] == "paired_crosslink_candidate"
        ]

    def test_import_preserves_every_source_row(self):
        bundle = self.read_bundle()
        self.assertEqual(bundle["format"], "kojak_csm")
        self.assertEqual(len(bundle["source_columns"]), 16)
        self.assertEqual(len(bundle["records"]), 344)
        self.assertEqual(
            Counter(record["record_kind"] for record in bundle["records"]),
            {
                "paired_crosslink_candidate": 149,
                "single_peptide_record": 30,
                "no_peptide_reported": 157,
                "loop_link_candidate": 8,
            },
        )
        self.assertEqual(bundle["parser"]["backend"], "native")
        self.assertIsNone(bundle["parser"].get("library"))
        self.assertEqual(
            len({record["source_record_id"] for record in bundle["records"]}),
            344,
        )
        for actual, expected in zip(
            bundle["records"], self.expected["records"]
        ):
            for field in (
                "source_record_id", "source_record_number",
                "source_fields", "source_values", "record_kind",
            ):
                self.assertEqual(actual[field], expected[field])
        self.assertTrue(all(
            record["record_status"] == "valid"
            for record in self.paired(bundle)
        ))

    def test_decoy_annotations_are_not_promoted_to_targets(self):
        bundle = self.read_bundle()
        paired = self.paired(bundle)
        self.assertEqual(
            Counter(record["decoy_status"] for record in paired),
            {
                "decoy_only": 31,
                "target_labels_only": 44,
                "contains_target_and_decoy_labels": 74,
            },
        )
        for actual, expected in zip(
            bundle["records"], self.expected["records"]
        ):
            if expected["record_kind"] != "paired_crosslink_candidate":
                continue
            for side in ("a", "b"):
                actual_candidates = actual["endpoints"][side][
                    "protein_candidates"
                ]
                expected_candidates = expected["endpoints"][side][
                    "protein_candidates"
                ]
                self.assertEqual(
                    [
                        (item["protein_id"], item["protein_position"])
                        for item in actual_candidates
                    ],
                    [
                        (item["protein_id"], item["protein_position"])
                        for item in expected_candidates
                    ],
                )
                for item in actual_candidates:
                    self.assertEqual(
                        item["decoy_label_present"],
                        item["protein_id"].startswith("decoy_"),
                    )

    def test_modifications_coordinates_and_non_lysine_site_are_preserved(self):
        bundle = self.read_bundle()
        residues = Counter()
        modification_tokens = Counter()
        for actual, expected in zip(
            bundle["records"], self.expected["records"]
        ):
            if expected["record_kind"] != "paired_crosslink_candidate":
                continue
            for side in ("a", "b"):
                endpoint = actual["endpoints"][side]
                wanted = expected["endpoints"][side]
                for field in (
                    "raw_peptide", "peptide", "link_position",
                    "modifications", "protein_candidates",
                ):
                    self.assertEqual(endpoint[field], wanted[field])
                residues.update(
                    item["residue"]
                    for item in endpoint["protein_candidates"]
                )
                modification_tokens.update(
                    item["raw_token"] for item in endpoint["modifications"]
                )
                for item in endpoint["protein_candidates"]:
                    self.assertEqual(
                        item["peptide_start"],
                        item["protein_position"] - endpoint["link_position"] + 1,
                    )
            self.assertEqual(actual["spectrum_file"], None)
            self.assertEqual(actual["fdr_status"], "not_provided")
        self.assertEqual(residues, {"K": 297, "L": 1})
        self.assertEqual(
            dict(modification_tokens),
            self.expected["summary"]["paired_modification_tokens"],
        )

    def test_nonpaired_rows_do_not_acquire_beta_peptides(self):
        bundle = self.read_bundle()
        for record in bundle["records"]:
            if record["record_kind"] == "paired_crosslink_candidate":
                continue
            self.assertEqual(record["source_fields"]["Peptide #2"], "-")
            self.assertNotIn("b", record.get("endpoints", {}))

    def test_unresolved_protein_annotation_is_not_discarded_by_first_match(self):
        baseline = next(
            record for record in self.expected["records"]
            if record["record_kind"] == "paired_crosslink_candidate"
            and record["decoy_status"] == "target_labels_only"
        )
        raw_rows = list(csv.reader(
            io.StringIO(self.results.read_text(encoding="utf-8-sig")),
            delimiter="\t",
        ))
        columns = raw_rows[1]
        number = baseline["source_record_number"]
        row = raw_rows[number + 1]
        self.assertEqual(row, baseline["source_values"])
        column = columns.index("Protein #1")
        row[column] = row[column].rstrip(";") + ";missing_protein(1);"

        with tempfile.TemporaryDirectory() as temporary:
            changed = Path(temporary) / "changed.kojak.txt"
            with changed.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, delimiter="\t").writerows(raw_rows)
            bundle = self.read_bundle(results=changed)

        self.assertEqual(len(bundle["records"]), 344)
        record = bundle["records"][number - 1]
        self.assertEqual(record["source_fields"]["Protein #1"], row[column])
        self.assertNotEqual(record["record_status"], "valid")
        self.assertEqual(
            len(record["endpoints"]["a"]["protein_candidates"]), 2
        )
        self.assertEqual(sum(
            item["record_status"] == "valid"
            for item in self.paired(bundle)
        ), 148)

    def test_missing_decoy_fasta_entries_are_not_resolved_as_targets(self):
        baseline = next(
            record for record in self.expected["records"]
            if record["record_kind"] == "paired_crosslink_candidate"
            and record["decoy_status"] == "decoy_only"
        )
        lines = []
        keep = False
        for line in self.fasta.read_text(encoding="utf-8-sig").splitlines():
            if line.startswith(">"):
                keep = not line[1:].split()[0].startswith("decoy_")
            if keep:
                lines.append(line)

        with tempfile.TemporaryDirectory() as temporary:
            target_only = Path(temporary) / "target_only.fasta"
            target_only.write_text(
                "\n".join(lines) + "\n", encoding="utf-8"
            )
            bundle = self.read_bundle(fasta=target_only)

        self.assertEqual(len(bundle["records"]), 344)
        record = bundle["records"][baseline["source_record_number"] - 1]
        self.assertEqual(record["decoy_status"], "decoy_only")
        self.assertNotEqual(record["record_status"], "valid")
        for endpoint in record["endpoints"].values():
            for item in endpoint["protein_candidates"]:
                self.assertTrue(item["protein_id"].startswith("decoy_"))
                self.assertEqual(
                    item["reason"], "protein_not_in_search_fasta"
                )
        self.assertEqual(sum(
            item["record_status"] == "valid"
            for item in self.paired(bundle)
        ), 44)

class TestKojakAuditCLI(unittest.TestCase):
    def setUp(self):
        self.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures/kojak_pxd015044"
        )
        expected = json.loads(
            (self.fixture / "expected.json").read_text(encoding="utf-8")
        )
        self.results = self.fixture / expected["result_filename"]
        self.fasta = self.fixture / "search.fasta"
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def arguments(self, destination, *, include_decoy=True):
        arguments = [
            "standardize",
            "--results", str(self.results),
            "--fasta", str(self.fasta),
            "--dataset-accession", "PXD015044",
            "--output", str(destination),
        ]
        if include_decoy:
            arguments.extend(["--decoy-prefix", "decoy_"])
        return arguments

    def test_cli_audit_matches_the_native_import(self):
        import hashlib

        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main
        from claudio.dataset_catalog.kojak import read_kojak_csms

        destination = self.root / "audit"
        result = CliRunner().invoke(
            main, self.arguments(destination) + ["--audit-only"]
        )
        self.assertEqual(
            result.exit_code, 0,
            result.output + repr(result.exception),
        )
        report = json.loads(
            (destination / "import_report.json").read_text(encoding="utf-8")
        )
        expected = read_kojak_csms(
            self.results, self.fasta,
            dataset_accession="PXD015044",
            decoy_prefix="decoy_",
        )
        self.assertEqual(report, expected)
        manifest = json.loads(
            (destination / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["counts"]["source_records"], 344)
        self.assertEqual(
            manifest["counts"]["validated_source_records"], 149
        )
        self.assertIs(manifest["fdr_filter_applied"], False)
        self.assertIs(manifest["identity_mapping_applied"], False)
        self.assertIs(manifest["input_exported"], False)
        self.assertIs(manifest["assessment_performed"], False)
        self.assertFalse((destination / "claudio_input.csv").exists())
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(
                hashlib.sha256((destination / name).read_bytes()).hexdigest(),
                digest,
            )
        self.assertIn("Kojak", result.output)

    def test_cli_requires_audit_mode_and_rejects_selection_options(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main

        cases = [
            ([], "--audit-only"),
            (["--audit-only", "--residue", "K"], "--residue"),
            (
                ["--audit-only", "--crosslinker", "CBDPS"],
                "--crosslinker",
            ),
            (
                ["--audit-only", "--fdr-level", "psm",
                 "--fdr-cutoff", "0.01"],
                "--fdr-level",
            ),
            (
                ["--audit-only", "--export-unverified-candidates"],
                "--export-unverified-candidates",
            ),
        ]
        for number, (options, expected_option) in enumerate(cases):
            with self.subTest(options=options):
                destination = self.root / f"rejected_{number}"
                result = CliRunner().invoke(
                    main, self.arguments(destination) + options
                )
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(expected_option, result.output)
                self.assertFalse(destination.exists())

    def test_cli_requires_an_explicit_decoy_prefix(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main

        destination = self.root / "missing_prefix"
        result = CliRunner().invoke(
            main,
            self.arguments(destination, include_decoy=False)
            + ["--audit-only"],
        )
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("--decoy-prefix", result.output)
        self.assertFalse(destination.exists())

    def test_changed_sources_are_rejected_before_audit_output(self):
        from claudio.dataset_catalog.kojak import read_kojak_csms
        from claudio.dataset_catalog.xifdr import write_import_report

        results = self.root / "source.kojak.txt"
        fasta = self.root / "search.fasta"
        results.write_bytes(self.results.read_bytes())
        fasta.write_bytes(self.fasta.read_bytes())
        bundle = read_kojak_csms(
            results, fasta,
            dataset_accession="PXD015044",
            decoy_prefix="decoy_",
        )
        for name, path in (("psm", results), ("fasta", fasta)):
            with self.subTest(source=name):
                original = path.read_bytes()
                destination = self.root / ("changed_" + name)
                try:
                    path.write_bytes(original + b"\n")
                    with self.assertRaisesRegex(
                        ValueError, "(?i)changed|checksum"
                    ):
                        write_import_report(bundle, destination)
                    self.assertFalse(destination.exists())
                finally:
                    path.write_bytes(original)
