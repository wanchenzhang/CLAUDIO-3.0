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

class TestKojakInputExport(unittest.TestCase):
    def setUp(self):
        from claudio.dataset_catalog.cli import _read_reference_projection
        from claudio.dataset_catalog.kojak import read_kojak_csms

        self.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures/kojak_pxd015044"
        )
        native = json.loads(
            (self.fixture / "expected.json").read_text(encoding="utf-8")
        )
        self.expected = json.loads(
            (self.fixture / "evidence/expected_input.json").read_text(
                encoding="utf-8"
            )
        )
        self.bundle = read_kojak_csms(
            self.fixture / native["result_filename"],
            self.fixture / "search.fasta",
            dataset_accession="PXD015044",
            decoy_prefix="decoy_",
        )
        self.mapping = _read_reference_projection(
            self.fixture / "evidence/reference_mapping.json"
        )
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def export(self, *, bundle=None, mapping=None, name="input", **options):
        from claudio.dataset_catalog import kojak

        writer = getattr(kojak, "write_kojak_claudio_input", None)
        self.assertTrue(
            callable(writer), "Missing Kojak CLAUDIO input exporter"
        )
        arguments = {
            "allowed_linker_masses": {"509.0974", "517.1476"},
            "allowed_residues": {"K"},
            "export_unverified_candidates": True,
        }
        arguments.update(options)
        return Path(writer(
            self.bundle if bundle is None else bundle,
            self.mapping if mapping is None else mapping,
            self.root / name,
            **arguments,
        ))

    def rows(self, destination):
        with (destination / "claudio_input.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            reader = csv.DictReader(handle)
            return reader.fieldnames, list(reader)

    def test_exports_target_pairs_and_preserves_all_source_rows(self):
        import copy

        original = copy.deepcopy(self.bundle)
        destination = self.export()
        columns, rows = self.rows(destination)
        self.assertEqual(columns[:8], [
            "peptide1", "peptide2", "position1", "position2",
            "k_pos1", "k_pos2", "entry1", "entry2",
        ])
        self.assertEqual(len(rows), 44)
        for actual, expected in zip(rows, self.expected["rows"]):
            for key, value in expected.items():
                self.assertEqual(actual[key], value)
            self.assertEqual(actual["input_scope"], "unverified_candidate")
            self.assertEqual(actual["fdr_verification_status"], "unverified")
            self.assertEqual(
                actual["source_decoy_status"], "target_labels_only"
            )
            for key in (
                "selected_fdr_level", "selected_fdr_value",
                "selected_fdr_cutoff",
            ):
                self.assertEqual(actual[key], "")

        audit = json.loads(
            (destination / "psm_audit.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(audit), 344)
        self.assertEqual(
            Counter(item["input_status"] for item in audit),
            {"candidate_ready": 44, "not_evaluable": 300},
        )
        for actual, source in zip(audit, self.bundle["records"]):
            for key in (
                "source_record_id", "source_record_number",
                "source_fields", "source_values", "record_kind",
                "decoy_status", "fdr",
            ):
                self.assertEqual(actual[key], source[key])
            if (
                source["record_kind"] == "paired_crosslink_candidate"
                and source["decoy_status"] != "target_labels_only"
            ):
                self.assertIn("decoy_record", actual["input_reasons"])

        manifest = json.loads(
            (destination / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["kind"], "kojak_claudio_input")
        self.assertIs(
            manifest["fdr_policy"]["per_record_fdr_applied"], False
        )
        self.assertIs(manifest["assessment_performed"], False)
        self.assertEqual(self.bundle, original)

    def test_export_requires_explicit_unverified_candidate_mode(self):
        with self.assertRaisesRegex(ValueError, "(?i)unverified"):
            self.export(export_unverified_candidates=False)
        self.assertFalse((self.root / "input").exists())

    def test_linker_mass_selection_preserves_native_isotope_values(self):
        destination = self.export(allowed_linker_masses={"509.0974"})
        _, rows = self.rows(destination)
        self.assertEqual(len(rows), 21)
        self.assertEqual(
            {row["source_linker_mass"] for row in rows}, {"509.0974"}
        )

    def test_missing_construct_review_is_rejected_before_output(self):
        import copy

        mapping = copy.deepcopy(self.mapping)
        del mapping["references"]["TAU-441-His"]["sequence_projection_review"]
        with self.assertRaisesRegex(ValueError, "(?i)review"):
            self.export(mapping=mapping)
        self.assertFalse((self.root / "input").exists())

    def test_changed_reference_evidence_is_rejected_before_output(self):
        receipt = next(
            source
            for source in self.mapping["references"]["TAU-441-His"]["sources"]
            if source["request_url"].endswith("/P10636-8.fasta")
        )
        path = Path(receipt["absolute_path"])
        original = path.read_bytes()
        try:
            path.write_bytes(original + b"\n")
            with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
                self.export()
            self.assertFalse((self.root / "input").exists())
        finally:
            path.write_bytes(original)

    def test_changed_decoy_evidence_is_rejected_before_output(self):
        import copy

        bundle = copy.deepcopy(self.bundle)
        record = next(
            item for item in bundle["records"]
            if item["record_kind"] == "paired_crosslink_candidate"
            and item["decoy_status"] == "decoy_only"
        )
        record["decoy_status"] = "target_labels_only"
        with self.assertRaisesRegex(ValueError, "(?i)native|source|changed"):
            self.export(bundle=bundle)
        self.assertFalse((self.root / "input").exists())

    def test_packaged_reference_evidence_can_be_loaded_after_copying(self):
        import shutil
        from claudio.dataset_catalog.cli import _read_reference_projection

        destination = self.export()
        moved = self.root / "moved"
        shutil.copytree(destination, moved)
        mapping = _read_reference_projection(moved / "reference_projection.json")

        def check_receipts(value):
            if isinstance(value, dict):
                if (
                    "kind" in value and "sha256" in value
                    and "absolute_path" in value
                ):
                    self.assertTrue(
                        Path(value["absolute_path"]).is_relative_to(moved)
                    )
                for child in value.values():
                    check_receipts(child)
            elif isinstance(value, list):
                for child in value:
                    check_receipts(child)

        check_receipts(mapping)
        repeated = self.export(mapping=mapping, name="repeated")
        self.assertEqual(
            (destination / "claudio_input.csv").read_bytes(),
            (repeated / "claudio_input.csv").read_bytes(),
        )

class TestKojakExportCLI(unittest.TestCase):
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
        self.mapping = self.fixture / "evidence/reference_mapping.json"
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def arguments(self, destination):
        return [
            "standardize",
            "--results", str(self.results),
            "--fasta", str(self.fasta),
            "--dataset-accession", "PXD015044",
            "--decoy-prefix", "decoy_",
            "--reference-projection", str(self.mapping),
            "--kojak-linker-mass", "509.0974",
            "--kojak-linker-mass", "517.1476",
            "--residue", "K",
            "--export-unverified-candidates",
            "--output", str(destination),
        ]

    def test_cli_export_matches_the_verified_api_behavior(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import (
            main, _read_reference_projection,
        )
        from claudio.dataset_catalog.kojak import (
            read_kojak_csms, write_kojak_claudio_input,
        )

        bundle = read_kojak_csms(
            self.results, self.fasta,
            dataset_accession="PXD015044", decoy_prefix="decoy_",
        )
        api = Path(write_kojak_claudio_input(
            bundle, _read_reference_projection(self.mapping),
            self.root / "api",
            allowed_linker_masses={"509.0974", "517.1476"},
            allowed_residues={"K"},
            export_unverified_candidates=True,
        ))
        destination = self.root / "cli"
        result = CliRunner().invoke(main, self.arguments(destination))
        self.assertEqual(
            result.exit_code, 0, result.output + repr(result.exception)
        )
        self.assertEqual(
            (destination / "claudio_input.csv").read_bytes(),
            (api / "claudio_input.csv").read_bytes(),
        )
        self.assertEqual(
            json.loads(
                (destination / "psm_audit.json").read_text(encoding="utf-8")
            ),
            json.loads(
                (api / "psm_audit.json").read_text(encoding="utf-8")
            ),
        )
        self.assertIn("unverified", result.output)

    def test_export_requires_all_explicit_candidate_options(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main

        required = (
            "--reference-projection",
            "--kojak-linker-mass",
            "--residue",
            "--export-unverified-candidates",
            "--decoy-prefix",
        )
        for number, missing in enumerate(required):
            with self.subTest(missing=missing):
                destination = self.root / f"missing_{number}"
                arguments = self.arguments(destination)
                filtered = []
                index = 0
                while index < len(arguments):
                    if arguments[index] == missing:
                        index += (
                            1 if missing == "--export-unverified-candidates"
                            else 2
                        )
                    else:
                        filtered.append(arguments[index])
                        index += 1
                result = CliRunner().invoke(main, filtered)
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(missing, result.output)
                self.assertFalse(destination.exists())

    def test_changed_reference_evidence_is_rejected_before_output(self):
        import shutil
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main

        evidence = self.root / "evidence"
        shutil.copytree(self.fixture / "evidence", evidence)
        mapping = json.loads(
            (evidence / "reference_mapping.json").read_text(encoding="utf-8")
        )
        receipt = next(
            item
            for item in mapping["references"]["TAU-441-His"]["sources"]
            if item["request_url"].endswith("/P10636-8.fasta")
        )
        path = evidence / receipt["path"]
        path.write_bytes(path.read_bytes() + b"\n")
        destination = self.root / "rejected"
        arguments = self.arguments(destination)
        index = arguments.index("--reference-projection") + 1
        arguments[index] = str(evidence / "reference_mapping.json")
        result = CliRunner().invoke(main, arguments)
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("checksum", result.output.lower())
        self.assertFalse(destination.exists())

    def test_linker_mass_option_is_restricted_to_kojak_export(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main

        other_results = (
            self.fixture.parent / "merox_pxd010895/velos.zhrm"
        )
        self.assertTrue(other_results.is_file())
        for number, (results, audit_only) in enumerate((
            (self.results, True),
            (other_results, False),
        )):
            with self.subTest(results=results, audit_only=audit_only):
                destination = self.root / f"scope_{number}"
                arguments = [
                    "standardize",
                    "--results", str(results),
                    "--fasta", str(self.fasta),
                    "--dataset-accession", "PXD015044",
                    "--decoy-prefix", "decoy_",
                    "--kojak-linker-mass", "509.0974",
                    "--output", str(destination),
                ]
                if audit_only:
                    arguments.append("--audit-only")
                result = CliRunner().invoke(main, arguments)
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn("--kojak-linker-mass", result.output)
                self.assertFalse(destination.exists())
