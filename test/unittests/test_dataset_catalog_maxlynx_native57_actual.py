import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from click.testing import CliRunner
from claudio.dataset_catalog import cli, maxlynx


FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures" / "maxlynx_pxd053456"
)


class TestMaxLynxNative57Actual(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.bundle = maxlynx.read_maxlynx_csms(
            FIXTURE / "DSBU_crosslinkMsms.txt",
            FIXTURE / "search.fasta",
            dataset_accession="PXD053456",
            parameters_path=FIXTURE / "parameters.xml",
            backend="pyxlms",
        )

    def test_actual_records_and_source_values_are_preserved(self):
        with (FIXTURE / "DSBU_crosslinkMsms.txt").open(
            newline="", encoding="utf-8-sig"
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames
            rows = list(reader)

        records = self.bundle["records"]
        self.assertEqual(len(records), 3228)
        self.assertEqual(self.bundle["parser"]["parsed_records"], 42)
        self.assertEqual(
            [record["source_fields"] for record in records], rows
        )
        self.assertEqual(
            [record["source_values"] for record in records],
            [[row[column] for column in columns] for row in rows],
        )
        self.assertEqual(
            len({record["source_record_id"] for record in records}), 3228
        )
        self.assertEqual(self.bundle["native_settings"]["crosslinker"], "DSBU")

        policy = self.bundle["fdr_policy"]
        self.assertEqual(policy["verification_status"], "configured_only")
        self.assertEqual(policy["configured_cutoff"], 0.01)
        self.assertFalse(policy["per_record_fdr_applied"])
        self.assertFalse(policy["application_verified"])

        examples = [
            record for record in records
            if record["source_fields"]["Scan number"] == "1549"
            and record["source_fields"]["Sequence1"] == "GKSQAK"
            and record["source_fields"]["Sequence2"] == "LIKSR"
        ]
        self.assertEqual(len(examples), 1)
        record = examples[0]
        self.assertEqual(record["record_status"], "valid")
        self.assertEqual(record["endpoints"]["a"]["protein_position"], 368)
        self.assertEqual(record["endpoints"]["a"]["link_position"], 2)
        self.assertEqual(record["endpoints"]["a"]["peptide_start"], 367)
        self.assertEqual(record["endpoints"]["b"]["protein_position"], 364)
        self.assertEqual(record["endpoints"]["b"]["link_position"], 3)
        self.assertEqual(record["endpoints"]["b"]["peptide_start"], 362)

    def test_cli_reports_native_layout_and_matches_api_audit(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "audit"
            result = CliRunner().invoke(cli.main, [
                "standardize",
                "--results", str(FIXTURE / "DSBU_crosslinkMsms.txt"),
                "--fasta", str(FIXTURE / "search.fasta"),
                "--dataset-accession", "PXD053456",
                "--maxlynx-parameters", str(FIXTURE / "parameters.xml"),
                "--audit-only",
                "--output", str(output),
            ])
            self.assertEqual(
                result.exit_code, 0,
                result.output + repr(result.exception),
            )
            self.assertIn(
                "Detected format: MaxLynx CSM (native 57-column TSV)",
                result.output,
            )
            actual = json.loads(
                (output / "import_report.json").read_text(encoding="utf-8")
            )
            self.assertEqual(actual, self.bundle)

    def test_changed_parameter_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            parameters = Path(temporary) / "parameters.xml"
            original = (FIXTURE / "parameters.xml").read_bytes()
            parameters.write_bytes(original)
            bundle = maxlynx.read_maxlynx_csms(
                FIXTURE / "DSBU_crosslinkMsms.txt",
                FIXTURE / "search.fasta",
                dataset_accession="PXD053456",
                parameters_path=parameters,
                backend="native",
            )
            parameters.write_bytes(original + b"\n")
            with self.assertRaisesRegex(ValueError, "(?i)checksum"):
                maxlynx.verify_maxlynx_sources(bundle)


if __name__ == "__main__":
    unittest.main()
