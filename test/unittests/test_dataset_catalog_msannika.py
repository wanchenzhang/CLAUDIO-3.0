import copy
import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from claudio.dataset_catalog import xifdr


class TestMSAnnikaCSMs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures" / "msannika_pxd041955"
        )
        reader = getattr(xifdr, "read_msannika_csms", None)
        if not callable(reader):
            raise AssertionError("Missing MS Annika CSM importer")
        original_reader = pd.read_excel
        cls.bundle = cls.read_fixture()
        if pd.read_excel is not original_reader:
            raise AssertionError("Importer did not restore pandas.read_excel")

    @classmethod
    def read_fixture(cls, root=None, **options):
        root = cls.fixture if root is None else root
        return xifdr.read_msannika_csms(
            root / "csms.xlsx",
            root / "search.fasta",
            dataset_accession="PXD041955",
            crosslinks_path=root / "crosslinks.xlsx",
            contaminant_fasta_path=root / "contaminants.fasta",
            **options,
        )

    def test_source_levels_and_every_record_are_retained(self):
        bundle = self.bundle
        self.assertEqual(bundle["format"], "msannika_csm")
        self.assertEqual(len(bundle["records"]), 765)
        self.assertEqual(len(bundle["crosslink_records"]), 369)
        self.assertEqual(
            sum(record["record_status"] == "valid"
                for record in bundle["records"]),
            761,
        )
        self.assertEqual(bundle["parser"]["parsed_records"], 761)
        self.assertFalse(bundle["parser"]["modification_masses_assessed"])
        self.assertEqual(
            len({record["source_record_id"] for record in bundle["records"]}),
            765,
        )
        self.assertEqual(
            bundle["crosslink_records"][0]["source_fields"]["Confidence"],
            "High",
        )
        for record in bundle["records"]:
            self.assertTrue(all(
                value["status"] == "missing"
                for value in record["fdr"].values()
            ))

    def test_coordinates_and_original_source_values_are_preserved(self):
        first = self.bundle["records"][0]
        self.assertEqual(first["source_worksheet_row"], 2)
        self.assertEqual(first["source_fields"]["A in protein"], "12")
        self.assertIs(first["source_fields"]["Checked"], False)
        self.assertEqual(first["endpoints"]["a"]["peptide_start"], 13)
        self.assertEqual(first["endpoints"]["a"]["protein_position"], 13)
        self.assertEqual(first["endpoints"]["a"]["link_position"], 1)
        for record in self.bundle["records"]:
            if record["record_status"] != "valid":
                continue
            for endpoint in record["endpoints"].values():
                sequence = self.bundle["proteins"][
                    endpoint["source_fasta_id"]
                ]["sequence"]
                start = endpoint["peptide_start"] - 1
                peptide = endpoint["peptide"]
                self.assertEqual(sequence[start:start + len(peptide)], peptide)
                self.assertEqual(
                    endpoint["protein_position"],
                    endpoint["peptide_start"] + endpoint["link_position"] - 1,
                )

    def test_incomplete_identifiers_are_not_promoted(self):
        unresolved = [
            record for record in self.bundle["records"]
            if record["record_status"] != "valid"
        ]
        self.assertEqual(
            [record["source_worksheet_row"] for record in unresolved],
            [445, 617, 672, 689],
        )
        for record in unresolved:
            for side, column in (("a", "Accession A"), ("b", "Accession B")):
                if record["source_fields"][column] == "sp":
                    endpoint = record["endpoints"][side]
                    self.assertIsNone(endpoint["source_fasta_id"])
                    self.assertEqual(endpoint["status"], "not_evaluable")

    def test_missing_voltage_and_source_aliases_remain_explicit(self):
        self.assertNotIn("Compensation Voltage", self.bundle["source_columns"])
        for record in self.bundle["records"]:
            self.assertNotIn("Compensation Voltage", record["source_fields"])
            self.assertEqual(
                record["compensation_voltage_status"], "not_present_in_source"
            )
        protein = self.bundle["proteins"]["sg|P0A7X3|RS9_ECOLI"]
        self.assertIsNone(protein["uniprot_accession"])
        self.assertEqual(
            sum(bool(record.get("parser_error"))
                for record in self.bundle["records"]),
            0,
        )

    def test_record_identity_is_independent_of_file_location(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in (
                "csms.xlsx", "crosslinks.xlsx",
                "search.fasta", "contaminants.fasta",
            ):
                shutil.copyfile(self.fixture / name, root / name)
            moved = self.read_fixture(root, backend="native")
            self.assertEqual(
                [record["source_record_id"] for record in moved["records"]],
                [record["source_record_id"]
                 for record in self.bundle["records"]],
            )

    def test_changed_secondary_sources_are_rejected_before_output(self):
        for name in ("crosslinks.xlsx", "contaminants.fasta"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                for source in (
                    "csms.xlsx", "crosslinks.xlsx",
                    "search.fasta", "contaminants.fasta",
                ):
                    shutil.copyfile(self.fixture / source, root / source)
                bundle = self.read_fixture(root, backend="native")
                path = root / name
                path.write_bytes(path.read_bytes() + b"\n")
                output = root / "audit"
                with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
                    xifdr.write_import_report(bundle, output)
                self.assertFalse(output.exists())

    def test_audit_writer_retains_both_result_levels(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "audit"
            xifdr.write_import_report(copy.deepcopy(self.bundle), output)
            report = json.loads((output / "import_report.json").read_text())
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(len(report["records"]), 765)
            self.assertEqual(len(report["crosslink_records"]), 369)
            self.assertEqual(manifest["counts"]["parser_errors"], 0)
            self.assertFalse(manifest["input_exported"])
            for name, expected in manifest["output_sha256"].items():
                self.assertEqual(
                    hashlib.sha256((output / name).read_bytes()).hexdigest(),
                    expected,
                )


from click.testing import CliRunner
from unittest.mock import patch

import openpyxl
from claudio.dataset_catalog import cli as catalog_cli


class TestMSAnnikaCLI(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures" / "msannika_pxd041955"
        )
        self.output = self.root / "audit"

    def arguments(self):
        return [
            "standardize",
            "--results", str(self.fixture / "csms.xlsx"),
            "--fasta", str(self.fixture / "search.fasta"),
            "--dataset-accession", "PXD041955",
            "--msannika-crosslinks", str(self.fixture / "crosslinks.xlsx"),
            "--msannika-contaminants", str(self.fixture / "contaminants.fasta"),
            "--output", str(self.output),
        ]

    def test_cli_retains_source_csms_and_separate_crosslinks(self):
        result = CliRunner().invoke(
            catalog_cli.main, self.arguments() + ["--audit-only"]
        )
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        bundle = json.loads((self.output / "import_report.json").read_text())
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertEqual(bundle["format"], "msannika_csm")
        self.assertEqual(len(bundle["records"]), 765)
        self.assertEqual(len(bundle["crosslink_records"]), 369)
        self.assertEqual(bundle["parser"]["parsed_records"], 761)
        self.assertEqual(manifest["counts"]["parser_errors"], 0)
        self.assertFalse(manifest["input_exported"])
        self.assertFalse(manifest["fdr_filter_applied"])

    def test_cli_requires_audit_mode(self):
        result = CliRunner().invoke(catalog_cli.main, self.arguments())
        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("--audit-only", result.output)
        self.assertFalse(self.output.exists())

    def test_cli_rejects_selection_options_during_audit(self):
        for extra in (
            ["--mapping", str(self.fixture / "manifest.json")],
            ["--fdr-level", "link", "--fdr-cutoff", "0.01"],
            ["--crosslinker", "DSSO"],
            ["--residue", "K"],
            ["--reference-projection", str(self.fixture / "manifest.json")],
            ["--allow-modification", "cm"],
            ["--export-unverified-candidates"],
        ):
            with self.subTest(extra=extra):
                result = CliRunner().invoke(
                    catalog_cli.main,
                    self.arguments() + ["--audit-only"] + extra,
                )
                self.assertEqual(result.exit_code, 2, result.output)
                self.assertIn(extra[0], result.output)
                self.assertFalse(self.output.exists())

    def test_secondary_options_are_rejected_for_other_formats(self):
        with patch.object(
            catalog_cli, "_detect_psm_format", return_value="xifdr_psm"
        ):
            result = CliRunner().invoke(
                catalog_cli.main, self.arguments() + ["--audit-only"]
            )
        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("MS Annika", result.output)
        self.assertFalse(self.output.exists())

    def test_mixed_result_workbook_is_rejected(self):
        workbook = openpyxl.Workbook()
        sheet = workbook.active
        sheet.append([
            "Sequence A", "Sequence B", "Accession A", "Accession B",
            "A in protein", "B in protein",
            "Crosslinker Position A", "Crosslinker Position B",
            "Alpha T/D", "Beta T/D", "Spectrum File", "First Scan",
            "Crosslink Strategy",
        ])
        workbook.create_sheet("pLink").append([
            "Title", "Peptide", "Proteins", "Peptide_Type", "Linker",
        ])
        path = self.root / "mixed.xlsx"
        workbook.save(path)
        workbook.close()
        with self.assertRaisesRegex(ValueError, "(?i)multiple|ambiguous|unique"):
            catalog_cli._detect_psm_format(path)

    def test_reader_retains_other_engine_rows_as_unresolved(self):
        path = self.root / "other_engine.xlsx"
        workbook = openpyxl.load_workbook(self.fixture / "csms.xlsx")
        sheet = workbook["CSMs"]
        columns = [cell.value for cell in sheet[1]]
        sheet.cell(2, columns.index("Crosslink Strategy") + 1, "XlinkX")
        workbook.save(path)
        workbook.close()
        bundle = xifdr.read_msannika_csms(
            path, self.fixture / "search.fasta",
            dataset_accession="PXD041955", backend="native",
        )
        self.assertEqual(
            bundle["records"][0]["record_reason"], "source_engine_unresolved"
        )
        self.assertEqual(
            bundle["records"][0]["source_fields"]["Crosslink Strategy"], "XlinkX"
        )
