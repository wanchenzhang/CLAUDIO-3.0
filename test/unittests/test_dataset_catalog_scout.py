import csv
import hashlib
import json
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

from claudio.dataset_catalog import xifdr
from click.testing import CliRunner
from claudio.dataset_catalog.cli import main
from claudio.dataset_catalog.xifdr import read_scout_csms
FIXTURE = Path(__file__).resolve().parents[1] / "fixtures" / "scout_pxd072720"


class TestScoutCSMs(unittest.TestCase):
    def reader(self):
        reader = getattr(xifdr, "read_scout_csms", None)
        self.assertTrue(callable(reader), "Missing Scout CSM importer")
        return reader

    def read_bundle(self, fasta=None):
        return self.reader()(
            FIXTURE / "csms.csv",
            fasta or FIXTURE / "search.fasta",
            dataset_accession="PXD072720",
            crosslinker="DSBU",
        )

    def test_import_restores_source_endpoints_and_retains_every_row(self):
        manifest = json.loads((FIXTURE / "manifest.json").read_text())
        for name in ("csms.csv", "search.fasta"):
            self.assertEqual(
                hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest(),
                manifest["files"][name]["sha256"],
            )

        with (FIXTURE / "csms.csv").open(
            newline="", encoding="utf-8-sig"
        ) as handle:
            original = list(csv.DictReader(handle))

        bundle = self.read_bundle()
        records = bundle["records"]

        self.assertEqual(bundle["format"], "scout_csm")
        self.assertEqual(len(records), 72)
        self.assertEqual(len({r["source_record_id"] for r in records}), 72)
        self.assertEqual(
            [r["source_record_number"] for r in records],
            list(range(1, 73)),
        )
        self.assertEqual(
            [r["source_fields"] for r in records], original
        )
        self.assertEqual(bundle["parser"]["parsed_records"], 63)

        paired = [
            r for r in records
            if r["record_kind"] == "paired_crosslink_candidate"
        ]
        self.assertEqual(len(paired), 63)
        self.assertEqual(
            Counter(r["parser_orientation"] for r in paired),
            {"same": 22, "swapped": 41},
        )
        self.assertTrue(all(r["record_status"] == "valid" for r in paired))

        first = paired[0]
        self.assertEqual(first["endpoints"]["a"]["peptide"], "VMGDKSTAR")
        self.assertEqual(first["endpoints"]["a"]["link_position"], 5)
        self.assertEqual(first["endpoints"]["a"]["protein_position"], 129)
        self.assertEqual(first["endpoints"]["b"]["peptide"], "VAMGEKLR")
        self.assertEqual(first["endpoints"]["b"]["protein_position"], 333)

        for record in paired:
            for endpoint in record["endpoints"].values():
                self.assertEqual(
                    endpoint["source_fasta_id"], "CAH2072050.1|BC|"
                )
                self.assertEqual(
                    endpoint["source_protein_id"],
                    "sp|CAH2072050_1_BC_|bc",
                )
                self.assertEqual(endpoint["status"], "valid")

        for record, row in zip(records, original):
            self.assertEqual(
                record["spectrum"]["scan_number"], int(row["Scan"])
            )
            self.assertEqual(
                record["spectrum"]["original_file"], row["File"]
            )
            self.assertEqual(record["fdr"]["psm"]["status"], "missing")

    def test_loop_links_keep_source_fields_without_fabricating_beta_peptides(self):
        bundle = self.read_bundle()
        loops = [
            r for r in bundle["records"]
            if r["record_kind"] == "loop_link"
        ]

        self.assertEqual(len(loops), 9)
        self.assertEqual(
            [r["source_record_number"] for r in loops],
            list(range(64, 73)),
        )
        for record in loops:
            self.assertEqual(record["record_status"], "not_evaluable")
            self.assertEqual(
                record["record_reason"],
                "loop_link_requires_single_peptide_support",
            )
            self.assertEqual(record["source_fields"]["Beta peptide"], "NA")
            self.assertEqual(
                record["source_fields"]["Beta peptide position"], "NA"
            )
            self.assertNotIn("b", record["endpoints"])
            self.assertNotIn("parser_error", record)
            self.assertNotIn("parser_orientation", record)
            self.assertIsNotNone(record["spectrum"])

    def test_identical_fasta_candidates_are_not_resolved_by_first_match(self):
        self.reader()
        proteins = xifdr.read_search_fasta(FIXTURE / "search.fasta")
        sequence = proteins["CAH2072050.1|BC|"]["sequence"]

        with TemporaryDirectory() as directory:
            fasta = Path(directory) / "ambiguous.fasta"
            fasta.write_bytes(
                (FIXTURE / "search.fasta").read_bytes()
                + b"\n>duplicate_BC\n"
                + sequence.encode("ascii")
                + b"\n"
            )
            bundle = self.read_bundle(fasta)

        first = bundle["records"][0]
        endpoint = first["endpoints"]["a"]
        self.assertEqual(first["record_status"], "not_evaluable")
        self.assertEqual(endpoint["status"], "ambiguous")
        self.assertIsNone(endpoint["source_fasta_id"])
        self.assertEqual(
            set(endpoint["source_fasta_candidates"]),
            {"CAH2072050.1|BC|", "duplicate_BC"},
        )
class TestScoutStandardizeCLI(unittest.TestCase):
    def arguments(self, output):
        return [
            "standardize",
            "--results", str(FIXTURE / "csms.csv"),
            "--fasta", str(FIXTURE / "search.fasta"),
            "--dataset-accession", "PXD072720",
            "--output", str(output),
        ]

    def test_cli_saves_all_source_records_and_loop_links(self):
        with TemporaryDirectory() as directory:
            output = Path(directory) / "audit"
            result = CliRunner().invoke(
                main,
                self.arguments(output)
                + ["--audit-only", "--crosslinker", "DSBU"],
            )
            self.assertEqual(
                result.exit_code, 0,
                result.output + repr(result.exception),
            )

            report = json.loads(
                (output / "import_report.json").read_text()
            )
            manifest = json.loads(
                (output / "manifest.json").read_text()
            )

            self.assertEqual(report["format"], "scout_csm")
            self.assertEqual(len(report["records"]), 72)
            self.assertEqual(report["parser"]["parsed_records"], 63)
            self.assertEqual(
                Counter(r["record_kind"] for r in report["records"]),
                {"paired_crosslink_candidate": 63, "loop_link": 9},
            )
            self.assertEqual(
                manifest["counts"],
                {
                    "source_records": 72,
                    "validated_source_records": 63,
                    "unresolved_source_records": 9,
                    "parser_errors": 0,
                },
            )
            for key in (
                "fdr_filter_applied",
                "identity_mapping_applied",
                "input_exported",
                "assessment_performed",
            ):
                self.assertFalse(manifest[key])

            for name, expected in manifest["output_sha256"].items():
                self.assertEqual(
                    hashlib.sha256(
                        (output / name).read_bytes()
                    ).hexdigest(),
                    expected,
                )
            self.assertFalse((output / "claudio_input.csv").exists())

    def test_cli_requires_one_explicit_crosslinker_for_scout(self):
        cases = [
            ["--audit-only"],
            [
                "--audit-only",
                "--crosslinker", "DSBU",
                "--crosslinker", "DSSO",
            ],
        ]
        for options in cases:
            with self.subTest(options=options), TemporaryDirectory() as directory:
                output = Path(directory) / "audit"
                result = CliRunner().invoke(
                    main, self.arguments(output) + options
                )
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn("--crosslinker", result.output)
                self.assertFalse(output.exists())

    def test_cli_keeps_export_and_selection_options_out_of_scout_audit(self):
        cases = [
            ["--crosslinker", "DSBU"],
            [
                "--audit-only", "--crosslinker", "DSBU",
                "--fdr-level", "psm", "--fdr-cutoff", "0.01",
            ],
            [
                "--audit-only", "--crosslinker", "DSBU",
                "--residue", "K",
            ],
        ]
        for options in cases:
            with self.subTest(options=options), TemporaryDirectory() as directory:
                output = Path(directory) / "audit"
                result = CliRunner().invoke(
                    main, self.arguments(output) + options
                )
                self.assertNotEqual(result.exit_code, 0)
                self.assertFalse(output.exists())
class TestScoutInputExport(unittest.TestCase):
    def setUp(self):
        import shutil

        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fixture = self.root / "fixture"
        shutil.copytree(FIXTURE, self.fixture)
        self.evidence = self.fixture / "evidence"
        self.output = self.root / "input"

        self.bundle = xifdr.read_scout_csms(
            self.fixture / "csms.csv",
            self.fixture / "search.fasta",
            dataset_accession="PXD072720",
            crosslinker="DSBU",
        )
        self.reports = json.loads(
            (self.evidence / "candidate_report.json").read_bytes()
        )
        for report in self.reports.values():
            for source in report["sources"]:
                source["absolute_path"] = str(
                    (self.evidence / source["path"]).resolve()
                )

        self.projection = xifdr.project_xifdr_candidates(
            self.bundle, {}, self.reports
        )

    def export(self):
        writer = getattr(xifdr, "write_scout_claudio_input", None)
        self.assertTrue(
            callable(writer), "Missing Scout CLAUDIO input exporter"
        )
        return writer(
            self.bundle,
            self.projection,
            self.output,
            author_filter_path=self.evidence / "manifest.json",
            allowed_crosslinkers={"DSBU"},
            allowed_residues={"K"},
        )

    def test_exports_projected_pairs_and_retains_every_source_record(self):
        destination = Path(self.export())

        with (destination / "claudio_input.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(
                reader.fieldnames[:8],
                [
                    "peptide1", "peptide2",
                    "position1", "position2",
                    "k_pos1", "k_pos2",
                    "entry1", "entry2",
                ],
            )
            self.assertNotIn("source_record_id", reader.fieldnames)
            rows = list(reader)

        paired = [
            record for record in self.bundle["records"]
            if record["record_kind"] == "paired_crosslink_candidate"
        ]
        self.assertEqual(len(paired), 63)
        expected = [
            record for record in paired
            if record["source_record_number"] != 47
        ]
        self.assertEqual(len(expected), 62)
        self.assertEqual(len(rows), 62)

        for row, original in zip(rows, expected):
            self.assertEqual(
                row["xlms_source_record_id"],
                original["source_record_id"],
            )
            self.assertEqual(
                row["identity_scope"], "candidate_hypothesis"
            )
            self.assertEqual(row.get("selected_fdr_value", ""), "")

            for site, suffix in (("a", "1"), ("b", "2")):
                endpoint = original["endpoints"][site]
                self.assertEqual(row["entry" + suffix], "A0AAU9SW86")
                self.assertEqual(
                    row["peptide" + suffix], endpoint["peptide"]
                )
                self.assertEqual(
                    int(row["position" + suffix]),
                    endpoint["protein_position"] + 57,
                )
                self.assertEqual(
                    int(row["k_pos" + suffix]),
                    endpoint["link_position"],
                )
                self.assertEqual(
                    int(row["source_protein_position_" + site]),
                    endpoint["protein_position"],
                )

        audit = json.loads(
            (destination / "psm_audit.json").read_bytes()
        )
        self.assertEqual(len(audit), 72)
        self.assertEqual(
            [record["source_record_id"] for record in audit],
            [
                record["source_record_id"]
                for record in self.bundle["records"]
            ],
        )
        self.assertEqual(
            [record["source_fields"] for record in audit],
            [
                record["source_fields"]
                for record in self.bundle["records"]
            ],
        )
        self.assertEqual(
            Counter(record["input_status"] for record in audit),
            {"ready": 62, "not_evaluable": 10},
        )
        excluded = next(
            record for record in audit
            if record["source_record_number"] == 47
        )
        self.assertEqual(excluded["input_status"], "not_evaluable")
        self.assertEqual(
            excluded["input_reasons"], ["residue_not_selected_a"]
        )
        self.assertEqual(
            excluded["endpoints"]["a"]["peptide"][
                excluded["endpoints"]["a"]["link_position"] - 1
            ],
            "T",
        )
        loops = [
            record for record in audit
            if record["record_kind"] == "loop_link"
        ]
        self.assertEqual(len(loops), 9)
        self.assertTrue(all(
            record["input_status"] == "not_evaluable"
            for record in loops
        ))
        self.assertTrue(all(
            record["fdr"]["psm"]["status"] == "missing"
            and record["fdr"]["psm"]["value"] is None
            for record in audit
        ))

        manifest = json.loads(
            (destination / "manifest.json").read_bytes()
        )
        self.assertEqual(manifest["kind"], "scout_claudio_input")
        self.assertEqual(manifest["counts"]["source_records"], 72)
        self.assertEqual(manifest["counts"]["pipeline_input_rows"], 62)
        self.assertEqual(
            manifest["author_filter_policy"]["reported_fdr_cutoff"],
            0.05,
        )
        self.assertEqual(
            manifest["author_filter_policy"]["reported_fdr_level"],
            "unspecified",
        )
        self.assertFalse(
            manifest["fdr_policy"]["per_record_fdr_applied"]
        )
        self.assertFalse(manifest["assessment_performed"])

        for name, expected in manifest["output_sha256"].items():
            self.assertEqual(
                hashlib.sha256(
                    (destination / name).read_bytes()
                ).hexdigest(),
                expected,
                name,
            )

    def test_changed_reference_evidence_is_rejected_before_output(self):
        report = next(iter(self.reports.values()))
        source = Path(report["sources"][0]["absolute_path"])
        source.write_bytes(source.read_bytes() + b"\n")

        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            self.export()

        self.assertFalse(self.output.exists())

    def test_changed_projection_coordinates_are_rejected_before_output(self):
        projected = next(
            record for record in self.projection["records"]
            if record["pairs"]
        )
        projected["pairs"][0]["reference_position_a"] += 1

        with self.assertRaisesRegex(
            ValueError, "(?i)projection|coordinate"
        ):
            self.export()

        self.assertFalse(self.output.exists())
    def cli_arguments(self):
        projection_path = self.root / "projection.json"
        projection_path.write_text(
            json.dumps(self.projection), encoding="utf-8"
        )
        return [
            "standardize",
            "--results", str(self.fixture / "csms.csv"),
            "--fasta", str(self.fixture / "search.fasta"),
            "--dataset-accession", "PXD072720",
            "--reference-projection", str(projection_path),
            "--author-filter-evidence",
            str(self.evidence / "manifest.json"),
            "--crosslinker", "DSBU",
            "--residue", "K",
            "--output", str(self.output),
        ]

    def test_cli_exports_the_same_input_as_the_verified_writer(self):
        result = CliRunner().invoke(main, self.cli_arguments())
        self.assertEqual(
            result.exit_code,
            0,
            result.output + repr(result.exception),
        )
        self.assertIn("Standardized input rows: 62", result.output)

        expected = xifdr.write_scout_claudio_input(
            self.bundle,
            self.projection,
            self.root / "api_input",
            author_filter_path=self.evidence / "manifest.json",
            allowed_crosslinkers={"DSBU"},
            allowed_residues={"K"},
        )
        self.assertEqual(
            (self.output / "claudio_input.csv").read_bytes(),
            (Path(expected) / "claudio_input.csv").read_bytes(),
        )

        manifest = json.loads(
            (self.output / "manifest.json").read_bytes()
        )
        self.assertEqual(manifest["kind"], "scout_claudio_input")
        self.assertEqual(manifest["counts"]["source_records"], 72)
        self.assertEqual(
            manifest["counts"]["pipeline_input_rows"], 62
        )
        self.assertFalse(
            manifest["fdr_policy"]["per_record_fdr_applied"]
        )
        self.assertFalse(manifest["assessment_performed"])

    def test_cli_rejects_changed_scout_reference_evidence(self):
        arguments = self.cli_arguments()
        report = next(iter(self.reports.values()))
        source = Path(report["sources"][0]["absolute_path"])
        source.write_bytes(source.read_bytes() + b"\n")

        result = CliRunner().invoke(main, arguments)

        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("checksum", result.output.lower())
        self.assertFalse(self.output.exists())

    def test_cli_requires_explicit_scout_export_evidence_and_residues(self):
        for option in (
            "--author-filter-evidence",
            "--reference-projection",
            "--residue",
        ):
            with self.subTest(missing_option=option):
                arguments = self.cli_arguments()
                index = arguments.index(option)
                del arguments[index:index + 2]

                result = CliRunner().invoke(main, arguments)

                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(option, result.output)
                self.assertFalse(self.output.exists())