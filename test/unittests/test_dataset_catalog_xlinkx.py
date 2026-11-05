import hashlib
import json
import shutil
import sqlite3
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

from claudio.dataset_catalog import xifdr
from click.testing import CliRunner
from claudio.dataset_catalog.cli import main

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "xlinkx_pxd021648"
)


class TestXlinkXCSMs(unittest.TestCase):
    def reader(self):
        reader = getattr(xifdr, "read_xlinkx_csms", None)
        self.assertTrue(callable(reader), "Missing XlinkX CSM importer")
        return reader

    def read_bundle(self, native=None):
        return self.reader()(
            native if native is not None else FIXTURE / "native.pdResult",
            FIXTURE / "search.fasta",
            dataset_accession="PXD021648",
        )

    @staticmethod
    def source_fields(row):
        return {
            name: (
                {"sqlite_blob_hex": value.hex()}
                if isinstance(value, bytes)
                else value
            )
            for name, value in dict(row).items()
        }

    def test_retains_native_rows_and_restores_source_endpoint_order(self):
        manifest = json.loads((FIXTURE / "manifest.json").read_bytes())
        for name, metadata in manifest["files"].items():
            self.assertEqual(
                hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest(),
                metadata["sha256"],
            )

        connection = sqlite3.connect(
            (FIXTURE / "native.pdResult").as_uri() + "?mode=ro",
            uri=True,
        )
        connection.row_factory = sqlite3.Row
        try:
            source = [
                self.source_fields(row)
                for row in connection.execute(
                    "SELECT * FROM CSMs ORDER BY WorkflowId, ID"
                )
            ]
            crosslinks = [
                self.source_fields(row)
                for row in connection.execute(
                    "SELECT * FROM Crosslinks ORDER BY WorkflowId, ID"
                )
            ]
        finally:
            connection.close()

        bundle = self.read_bundle()
        records = bundle["records"]
        self.assertEqual(bundle["format"], "xlinkx_csm")
        self.assertEqual(len(records), 156)
        self.assertEqual(bundle["parser"]["parsed_records"], 156)
        self.assertFalse(
            bundle["parser"]["modification_masses_assessed"]
        )
        self.assertEqual(
            bundle["native_table_counts"],
            {"CSMs": 156, "DecoyCSMs": 0, "Crosslinks": 89},
        )
        self.assertEqual(bundle["native_crosslinks"], crosslinks)
        self.assertTrue(all(
            "IsDecoy" not in row for row in bundle["native_crosslinks"]
        ))

        self.assertEqual(
            [record["source_fields"] for record in records], source
        )
        self.assertEqual(
            len({record["source_record_id"] for record in records}), 156
        )
        self.assertEqual(
            Counter(record["parser_orientation"] for record in records),
            {"same": 69, "swapped": 87},
        )

        for record, native in zip(records, source):
            self.assertEqual(record["source_table"], "CSMs")
            self.assertEqual(
                record["source_native_key"],
                {"WorkflowId": native["WorkflowId"], "ID": native["ID"]},
            )
            self.assertEqual(record["record_status"], "valid")
            self.assertEqual(record["decoy_status"], "target_target")
            self.assertIs(record["source_decoy_flag"], False)
            self.assertEqual(record["fdr"]["psm"]["status"], "missing")
            self.assertIsNone(record["fdr"]["psm"]["value"])
            self.assertEqual(
                record["spectrum"]["scan_number"], native["FirstScan"]
            )
            self.assertEqual(
                record["spectrum"]["original_file"],
                native["SpectrumFileName"],
            )

    def test_all_source_endpoints_match_the_search_fasta(self):
        bundle = self.read_bundle()
        checked = 0
        for record in bundle["records"]:
            native = record["source_fields"]
            self.assertEqual(
                record["record_kind"], "paired_crosslink_candidate"
            )
            for site, suffix in (("a", "A"), ("b", "B")):
                endpoint = record["endpoints"][site]
                self.assertEqual(endpoint["status"], "valid")
                self.assertEqual(
                    endpoint["source_fasta_id"], "sp|P02769|ALBU_BOVIN"
                )
                self.assertEqual(
                    endpoint["source_protein_id"],
                    native["ProteinAccession" + suffix],
                )
                self.assertEqual(
                    endpoint["peptide"], native["Sequence" + suffix]
                )
                self.assertEqual(
                    endpoint["link_position"],
                    int(native["CrosslinkerPosition" + suffix]),
                )
                self.assertEqual(
                    endpoint["protein_position"],
                    int(native["LeadingProteinPosition" + suffix]),
                )
                protein = bundle["proteins"][endpoint["source_fasta_id"]]
                start = endpoint["peptide_start"] - 1
                peptide = endpoint["peptide"]
                self.assertEqual(
                    protein["sequence"][start:start + len(peptide)],
                    peptide,
                )
                self.assertEqual(
                    peptide[endpoint["link_position"] - 1], "K"
                )
                checked += 1
        self.assertEqual(checked, 312)

    def test_decoy_table_rows_are_retained_without_target_promotion(self):
        with TemporaryDirectory() as temporary:
            native = Path(temporary) / "with_decoy.pdResult"
            shutil.copyfile(FIXTURE / "native.pdResult", native)

            connection = sqlite3.connect(native)
            try:
                columns = [
                    row[1] for row in connection.execute(
                        'PRAGMA table_info("DecoyCSMs")'
                    )
                ]

                def quoted(name):
                    return '"' + name.replace('"', '""') + '"'

                names = ", ".join(quoted(name) for name in columns)
                values = ", ".join(
                    "1" if name == "IsDecoy" else quoted(name)
                    for name in columns
                )
                connection.execute(
                    f'INSERT INTO "DecoyCSMs" ({names}) '
                    f'SELECT {values} FROM "CSMs" '
                    'ORDER BY WorkflowId, ID LIMIT 1'
                )
                connection.commit()
            finally:
                connection.close()

            bundle = self.read_bundle(native)
            records = bundle["records"]
            self.assertEqual(len(records), 157)
            self.assertEqual(
                len({record["source_record_id"] for record in records}), 157
            )
            decoys = [
                record for record in records
                if record["source_table"] == "DecoyCSMs"
            ]
            self.assertEqual(len(decoys), 1)
            self.assertIs(decoys[0]["source_decoy_flag"], True)
            self.assertEqual(decoys[0]["source_fields"]["IsDecoy"], 1)
            self.assertEqual(decoys[0]["decoy_status"], "unknown")
class TestXlinkXStandardizeCLI(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = Path(self.temporary.name) / "audit"

    def arguments(self):
        return [
            "standardize",
            "--results", str(FIXTURE / "native.pdResult"),
            "--fasta", str(FIXTURE / "search.fasta"),
            "--dataset-accession", "PXD021648",
            "--output", str(self.output),
        ]

    def test_cli_preserves_native_csms_and_separate_crosslinks(self):
        result = CliRunner().invoke(
            main, self.arguments() + ["--audit-only"]
        )
        self.assertEqual(
            result.exit_code,
            0,
            result.output + repr(result.exception),
        )
        self.assertIn("Source CSM records: 156", result.output)
        self.assertIn("Parsed CSM records: 156", result.output)

        report = json.loads(
            (self.output / "import_report.json").read_bytes()
        )
        self.assertEqual(report["format"], "xlinkx_csm")
        self.assertEqual(len(report["records"]), 156)
        self.assertEqual(
            report["native_table_counts"],
            {"CSMs": 156, "DecoyCSMs": 0, "Crosslinks": 89},
        )
        self.assertEqual(len(report["native_crosslinks"]), 89)
        self.assertEqual(
            Counter(
                record["parser_orientation"]
                for record in report["records"]
            ),
            {"same": 69, "swapped": 87},
        )
        self.assertEqual(
            len({
                record["source_record_id"]
                for record in report["records"]
            }),
            156,
        )
        self.assertTrue(all(
            record["fdr"]["psm"]["status"] == "missing"
            and record["fdr"]["psm"]["value"] is None
            for record in report["records"]
        ))

        manifest = json.loads(
            (self.output / "manifest.json").read_bytes()
        )
        self.assertEqual(
            manifest["counts"],
            {
                "source_records": 156,
                "validated_source_records": 156,
                "unresolved_source_records": 0,
                "parser_errors": 0,
            },
        )
        for flag in (
            "fdr_filter_applied",
            "identity_mapping_applied",
            "input_exported",
            "assessment_performed",
        ):
            self.assertIs(manifest[flag], False)

        for name, expected in manifest["output_sha256"].items():
            self.assertEqual(
                hashlib.sha256(
                    (self.output / name).read_bytes()
                ).hexdigest(),
                expected,
                name,
            )
        self.assertFalse(
            (self.output / "claudio_input.csv").exists()
        )

    def test_cli_requires_audit_mode_and_rejects_selection_options(self):
        cases = [
            ([], "--audit-only"),
            (["--audit-only", "--residue", "K"], "--residue"),
            (
                ["--audit-only", "--crosslinker", "DSBU"],
                "--crosslinker",
            ),
            (
                [
                    "--audit-only",
                    "--fdr-level", "psm",
                    "--fdr-cutoff", "0.05",
                ],
                "--fdr-level",
            ),
        ]
        for options, expected_option in cases:
            with self.subTest(options=options):
                result = CliRunner().invoke(
                    main, self.arguments() + options
                )
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(expected_option, result.output)
                self.assertFalse(self.output.exists())