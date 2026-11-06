import hashlib
import json
import shutil
import sqlite3
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory
import csv
from claudio.dataset_catalog import xifdr
from click.testing import CliRunner
from claudio.dataset_catalog.cli import main, _read_protein_mapping

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
    def export_arguments(self):
        return self.arguments() + [
            "--mapping",
            str(FIXTURE / "evidence/protein_mapping.json"),
            "--crosslinker", "DSBU",
            "--residue", "K",
        ]

    def test_cli_exports_the_same_input_as_the_verified_api(self):
        from claudio.dataset_catalog.cli import _read_protein_mapping
        from claudio.dataset_catalog.xifdr import (
            read_xlinkx_csms,
            write_xlinkx_claudio_input,
        )

        result = CliRunner().invoke(main, self.export_arguments())
        self.assertEqual(
            result.exit_code,
            0,
            result.output + repr(result.exception),
        )

        bundle = read_xlinkx_csms(
            FIXTURE / "native.pdResult",
            FIXTURE / "search.fasta",
            dataset_accession="PXD021648",
        )
        mappings = _read_protein_mapping(
            FIXTURE / "evidence/protein_mapping.json"
        )
        expected = Path(write_xlinkx_claudio_input(
            bundle,
            mappings,
            Path(self.temporary.name) / "api_input",
            allowed_crosslinkers={"DSBU"},
            allowed_residues={"K"},
        ))

        self.assertEqual(
            (self.output / "claudio_input.csv").read_bytes(),
            (expected / "claudio_input.csv").read_bytes(),
        )

        with (self.output / "claudio_input.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 156)
        self.assertTrue(all(
            row["fdr_verification_status"] == "unverified"
            and row["source_validation_strategy"] == "simple_validation"
            for row in rows
        ))
        self.assertTrue(all(
            not row.get(field)
            for row in rows
            for field in (
                "selected_fdr_level",
                "selected_fdr_value",
                "selected_fdr_cutoff",
            )
        ))

        manifest = json.loads(
            (self.output / "manifest.json").read_bytes()
        )
        self.assertEqual(manifest["kind"], "xlinkx_claudio_input")
        self.assertIs(manifest["assessment_performed"], False)
        for name, expected_hash in manifest["output_sha256"].items():
            self.assertEqual(
                hashlib.sha256(
                    (self.output / name).read_bytes()
                ).hexdigest(),
                expected_hash,
                name,
            )

        self.assertIn("Standardized input rows: 156", result.output)
        self.assertIn("Native FDR verification: unverified", result.output)
        self.assertIn(
            "Per-record FDR filtering: not performed",
            result.output,
        )

    def test_cli_export_requires_mapping_crosslinker_and_residue(self):
        for missing in ("--mapping", "--crosslinker", "--residue"):
            with self.subTest(missing=missing):
                arguments = self.export_arguments()
                index = arguments.index(missing)
                del arguments[index:index + 2]

                result = CliRunner().invoke(main, arguments)
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(missing, result.output)
                self.assertFalse(self.output.exists())

    def test_cli_rejects_changed_identity_evidence_before_output(self):
        import shutil

        from claudio.dataset_catalog.cli import _read_protein_mapping

        copied = Path(self.temporary.name) / "fixture"
        shutil.copytree(FIXTURE, copied)
        mapping_path = copied / "evidence/protein_mapping.json"
        mappings = _read_protein_mapping(mapping_path)

        receipts = [
            source
            for mapping in mappings.values()
            for source in mapping.get("sources", [])
        ]
        self.assertTrue(receipts)

        receipt_path = Path(receipts[0]["absolute_path"])
        self.assertTrue(receipt_path.is_relative_to(copied.resolve()))
        with receipt_path.open("ab") as handle:
            handle.write(b"\n")

        arguments = self.export_arguments()
        index = arguments.index("--mapping")
        arguments[index + 1] = str(mapping_path)

        result = CliRunner().invoke(main, arguments)
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("checksum", result.output.lower())
        self.assertFalse(self.output.exists())
class TestXlinkXInputExport(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fixture = self.root / "fixture"
        shutil.copytree(FIXTURE, self.fixture)
        self.native = self.fixture / "native.pdResult"
        self.output = self.root / "input"

        self.bundle = xifdr.read_xlinkx_csms(
            self.native,
            self.fixture / "search.fasta",
            dataset_accession="PXD021648",
        )
        self.mappings = _read_protein_mapping(
            self.fixture / "evidence/protein_mapping.json"
        )

    def export(self):
        writer = getattr(xifdr, "write_xlinkx_claudio_input", None)
        self.assertTrue(
            callable(writer), "Missing XlinkX CLAUDIO input exporter"
        )
        return Path(writer(
            self.bundle,
            self.mappings,
            self.output,
            allowed_crosslinkers={"DSBU"},
            allowed_residues={"K"},
        ))

    def read_rows(self, destination):
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
            return list(reader)

    def test_exports_exact_coordinates_with_unverified_fdr_labels(self):
        destination = self.export()
        rows = self.read_rows(destination)
        self.assertEqual(len(rows), 156)

        for row, original in zip(rows, self.bundle["records"]):
            self.assertEqual(
                row["xlms_source_record_id"],
                original["source_record_id"],
            )
            self.assertEqual(row["identity_scope"], "exact_sequence")
            self.assertEqual(
                row["selection_basis"], "native_reported_candidate"
            )
            self.assertEqual(
                row["fdr_verification_status"], "unverified"
            )
            self.assertEqual(
                row["source_validation_strategy"], "simple_validation"
            )
            self.assertEqual(
                row["source_native_csm_confidence"], "Unknown"
            )
            self.assertEqual(
                row["source_native_parent_confidence"], "High"
            )
            self.assertEqual(
                float(row["source_configured_fdr_cutoff"]), 0.05
            )
            for field in (
                "selected_fdr_level",
                "selected_fdr_value",
                "selected_fdr_cutoff",
            ):
                self.assertEqual(row.get(field, ""), "")

            for site, suffix in (("a", "1"), ("b", "2")):
                endpoint = original["endpoints"][site]
                self.assertEqual(row["entry" + suffix], "P02769")
                self.assertEqual(
                    row["peptide" + suffix], endpoint["peptide"]
                )
                self.assertEqual(
                    int(row["position" + suffix]),
                    endpoint["protein_position"],
                )
                self.assertEqual(
                    int(row["k_pos" + suffix]),
                    endpoint["link_position"],
                )

        audit = json.loads(
            (destination / "psm_audit.json").read_bytes()
        )
        self.assertEqual(len(audit), 156)
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
        self.assertTrue(all(
            record["input_status"] == "ready"
            and record["fdr"]["psm"]["status"] == "missing"
            and record["fdr"]["psm"]["value"] is None
            for record in audit
        ))

        manifest = json.loads(
            (destination / "manifest.json").read_bytes()
        )
        self.assertEqual(manifest["kind"], "xlinkx_claudio_input")
        self.assertEqual(manifest["counts"]["source_records"], 156)
        self.assertEqual(
            manifest["counts"]["pipeline_input_rows"], 156
        )
        policy = manifest["native_filter_policy"]
        self.assertEqual(policy["configured_fdr_cutoff"], 0.05)
        self.assertEqual(policy["configured_strategy"], "Percolator")
        self.assertEqual(policy["observed_strategy"], "simple_validation")
        self.assertEqual(policy["fdr_verification_status"], "unverified")
        self.assertTrue(any(
            "reverting to simple validation" in warning["Message"].lower()
            for warning in policy["runtime_warnings"]
        ))
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

    def test_changed_identity_evidence_is_rejected_before_output(self):
        identity = next(iter(self.mappings.values()))
        source = Path(identity["sources"][0]["absolute_path"])
        source.write_bytes(source.read_bytes() + b"\n")

        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            self.export()

        self.assertFalse(self.output.exists())

    def test_csm_without_parent_crosslink_remains_in_the_audit(self):
        connection = sqlite3.connect(self.native)
        try:
            missing_key = connection.execute(
                "SELECT CSMsWorkflowId, CSMsID FROM CrosslinksCSMs "
                "ORDER BY CSMsWorkflowId, CSMsID LIMIT 1"
            ).fetchone()
            connection.execute(
                "DELETE FROM CrosslinksCSMs "
                "WHERE CSMsWorkflowId = ? AND CSMsID = ?",
                missing_key,
            )
            connection.commit()
        finally:
            connection.close()

        self.bundle = xifdr.read_xlinkx_csms(
            self.native,
            self.fixture / "search.fasta",
            dataset_accession="PXD021648",
        )
        destination = self.export()
        rows = self.read_rows(destination)
        audit = json.loads(
            (destination / "psm_audit.json").read_bytes()
        )

        self.assertEqual(len(rows), 155)
        self.assertEqual(len(audit), 156)
        unresolved = [
            record for record in audit
            if record["input_status"] != "ready"
        ]
        self.assertEqual(len(unresolved), 1)
        self.assertEqual(
            unresolved[0]["source_native_key"],
            {"WorkflowId": missing_key[0], "ID": missing_key[1]},
        )
        self.assertIn(
            "native_parent_crosslink_missing",
            unresolved[0]["input_reasons"],
        )