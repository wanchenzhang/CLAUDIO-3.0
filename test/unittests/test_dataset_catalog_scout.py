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