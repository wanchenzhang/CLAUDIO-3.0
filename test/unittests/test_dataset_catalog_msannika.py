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

    def test_cli_export_requires_evidence_mapping_and_selection(self):
        complete = self.arguments() + [
            "--mapping", str(self.fixture / "evidence" / "protein_mapping.json"),
            "--author-filter-evidence",
            str(self.fixture / "evidence" / "author_filter.json"),
            "--crosslinker", "DSSO",
            "--residue", "K",
        ]
        for missing in (
            "--mapping", "--author-filter-evidence", "--crosslinker", "--residue"
        ):
            with self.subTest(missing=missing):
                arguments = list(complete)
                index = arguments.index(missing)
                del arguments[index:index + 2]
                result = CliRunner().invoke(catalog_cli.main, arguments)
                self.assertEqual(result.exit_code, 2, result.output)
                self.assertIn(missing, result.output)
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


from claudio.dataset_catalog.cli import _read_protein_mapping


class TestMSAnnikaInputExport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures" / "msannika_pxd041955"
        )
        cls.evidence = cls.fixture / "evidence"
        cls.writer = staticmethod(
            getattr(xifdr, "write_msannika_claudio_input", None)
        )
        if not callable(cls.writer):
            raise AssertionError("Missing MS Annika CLAUDIO input exporter")
        cls.bundle = xifdr.read_msannika_csms(
            cls.fixture / "csms.xlsx",
            cls.fixture / "search.fasta",
            dataset_accession="PXD041955",
            crosslinks_path=cls.fixture / "crosslinks.xlsx",
            contaminant_fasta_path=cls.fixture / "contaminants.fasta",
        )

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "input"
        self.mappings = _read_protein_mapping(
            self.evidence / "protein_mapping.json"
        )

    def export(self, bundle=None, policy=None):
        return Path(self.writer(
            copy.deepcopy(self.bundle) if bundle is None else bundle,
            self.mappings,
            self.output,
            author_filter_path=(
                self.evidence / "author_filter.json"
                if policy is None else policy
            ),
            allowed_crosslinkers={"DSSO"},
            allowed_residues={"K"},
        ))

    def load(self):
        import csv
        with (self.output / "claudio_input.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames
            rows = list(reader)
        audit = json.loads((self.output / "psm_audit.json").read_text())
        manifest = json.loads((self.output / "manifest.json").read_text())
        return columns, rows, audit, manifest

    def expected_records(self):
        return [
            record for record in self.bundle["records"]
            if record["record_status"] == "valid"
            and record["decoy_status"] == "target_target"
            and record["source_fields"]["Crosslinker"] == "DSSO"
            and all(
                endpoint["peptide"][endpoint["link_position"] - 1] == "K"
                for endpoint in record["endpoints"].values()
            )
        ]

    def test_exports_default_m1_columns_and_retains_all_source_rows(self):
        original = copy.deepcopy(self.bundle)
        self.export()
        columns, rows, audit, manifest = self.load()
        self.assertEqual(columns[:8], [
            "peptide1", "peptide2", "position1", "position2",
            "k_pos1", "k_pos2", "entry1", "entry2",
        ])
        self.assertEqual(
            [row["xlms_source_record_id"] for row in rows],
            [record["source_record_id"] for record in self.expected_records()],
        )
        self.assertTrue(rows)
        self.assertEqual(len(audit), 765)
        crosslinks = json.loads(
            (self.output / "crosslink_audit.json").read_text()
        )
        self.assertEqual(len(crosslinks), 369)
        for before, after in zip(self.bundle["records"], audit):
            for field in (
                "source_record_id", "source_fields", "source_values",
                "record_status", "decoy_status", "fdr",
            ):
                self.assertEqual(before[field], after[field])
        self.assertEqual(self.bundle, original)
        self.assertFalse(manifest["assessment_performed"])

        from claudio.module01.src.io.read_in import read_inputfile
        from claudio.module01.src.main import main as m1_main
        headers = next(
            parameter.default for parameter in m1_main.params
            if parameter.name == "projections"
        ).split(",")
        internal = [
            "pep_a", "pep_b", "pos_a", "pos_b",
            "res_pos_a", "res_pos_b", "unip_id_a", "unip_id_b",
        ]
        data = read_inputfile(
            str(self.output / "claudio_input.csv"),
            dict(zip(headers, internal)),
        )
        self.assertEqual(len(data), len(rows))

    def test_author_reported_confidence_does_not_become_numeric_csm_fdr(self):
        self.export()
        _, rows, _, manifest = self.load()
        for row in rows:
            self.assertEqual(row["identity_scope"], "exact_sequence")
            self.assertEqual(row["fdr_verification_status"], "author_reported_only")
            self.assertEqual(row["source_reported_fdr_level"], "unspecified")
            self.assertEqual(float(row["source_reported_fdr_cutoff"]), 0.01)
            for field in (
                "selected_fdr_level", "selected_fdr_value", "selected_fdr_cutoff"
            ):
                self.assertEqual(row[field], "")
        policy = manifest["fdr_policy"]
        self.assertEqual(policy["level"], "unspecified")
        self.assertEqual(policy["reported_cutoff"], 0.01)
        self.assertIs(policy["per_record_fdr_applied"], False)
        self.assertIs(policy["target_decoy_validation_performed"], False)
        self.assertTrue((self.output / "author_filter_evidence.json").is_file())

    def test_missing_parent_does_not_override_source_csm_evidence(self):
        membership = json.loads(
            (self.evidence / "membership_review.json").read_text()
        )
        missing = {
            record["source_record_id"] for record in membership["records"]
            if not record["parent_candidates"]["residue_pair"]
        }
        expected = {
            record["source_record_id"] for record in self.expected_records()
        }
        self.assertTrue(missing & expected)
        self.export()
        _, rows, _, _ = self.load()
        exported = {row["xlms_source_record_id"] for row in rows}
        self.assertTrue((missing & expected) <= exported)

    def test_changed_author_evidence_is_rejected_before_output(self):
        copied = self.root / "evidence"
        shutil.copytree(self.evidence, copied)
        policy = copied / "author_filter.json"
        policy.write_bytes(policy.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            self.export(policy=policy)
        self.assertFalse(self.output.exists())

    def test_changed_identity_response_is_rejected_before_output(self):
        receipt = next(
            source for item in self.mappings.values()
            for source in item["sources"]
        )
        changed = self.root / "changed_source"
        changed.write_bytes(Path(receipt["absolute_path"]).read_bytes() + b"\n")
        receipt["absolute_path"] = str(changed)
        with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
            self.export()
        self.assertFalse(self.output.exists())

    def test_changed_crosslink_source_is_rejected_before_output(self):
        copied = self.root / "sources"
        copied.mkdir()
        for name in (
            "csms.xlsx", "crosslinks.xlsx",
            "search.fasta", "contaminants.fasta",
        ):
            shutil.copyfile(self.fixture / name, copied / name)
        bundle = xifdr.read_msannika_csms(
            copied / "csms.xlsx", copied / "search.fasta",
            dataset_accession="PXD041955",
            crosslinks_path=copied / "crosslinks.xlsx",
            contaminant_fasta_path=copied / "contaminants.fasta",
            backend="native",
        )
        path = copied / "crosslinks.xlsx"
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
            self.export(bundle=bundle)
        self.assertFalse(self.output.exists())

    def test_missing_author_evidence_is_rejected_before_output(self):
        with self.assertRaisesRegex(ValueError, "(?i)author|evidence"):
            self.writer(
                copy.deepcopy(self.bundle), self.mappings, self.output,
                author_filter_path=None,
                allowed_crosslinkers={"DSSO"},
                allowed_residues={"K"},
            )
        self.assertFalse(self.output.exists())
class TestMSAnnikaExportCLI(unittest.TestCase):
    setUp = TestMSAnnikaCLI.setUp
    arguments = TestMSAnnikaCLI.arguments

    def export_options(self, mapping=None, evidence=None):
        return [
            "--mapping", str(
                self.fixture / "evidence" / "protein_mapping.json"
                if mapping is None else mapping
            ),
            "--author-filter-evidence", str(
                self.fixture / "evidence" / "author_filter.json"
                if evidence is None else evidence
            ),
            "--crosslinker", "DSSO",
            "--residue", "K",
        ]

    def test_cli_reproduces_the_verified_api_input(self):
        bundle = xifdr.read_msannika_csms(
            self.fixture / "csms.xlsx",
            self.fixture / "search.fasta",
            dataset_accession="PXD041955",
            crosslinks_path=self.fixture / "crosslinks.xlsx",
            contaminant_fasta_path=self.fixture / "contaminants.fasta",
        )
        mappings = _read_protein_mapping(
            self.fixture / "evidence" / "protein_mapping.json"
        )
        reference = self.root / "api_input"
        xifdr.write_msannika_claudio_input(
            bundle, mappings, reference,
            author_filter_path=self.fixture / "evidence" / "author_filter.json",
            allowed_crosslinkers={"DSSO"},
            allowed_residues={"K"},
        )
        result = CliRunner().invoke(
            catalog_cli.main, self.arguments() + self.export_options()
        )
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertEqual(
            (self.output / "claudio_input.csv").read_bytes(),
            (reference / "claudio_input.csv").read_bytes(),
        )
        self.assertIn("Standardized input rows: 759", result.output)
        manifest = json.loads((self.output / "manifest.json").read_text())
        self.assertEqual(manifest["fdr_policy"]["level"], "unspecified")
        self.assertFalse(manifest["fdr_policy"]["per_record_fdr_applied"])

    def test_export_rejects_per_record_fdr_options(self):
        for extra in (
            ["--fdr-level", "psm"],
            ["--fdr-cutoff", "0.01"],
            ["--fdr-column", "Score"],
        ):
            with self.subTest(extra=extra):
                result = CliRunner().invoke(
                    catalog_cli.main,
                    self.arguments() + self.export_options() + extra,
                )
                self.assertEqual(result.exit_code, 2, result.output)
                self.assertIn(extra[0], result.output)
                self.assertFalse(self.output.exists())

    def test_cli_reuses_a_moved_evidence_package(self):
        first = CliRunner().invoke(
            catalog_cli.main, self.arguments() + self.export_options()
        )
        self.assertEqual(first.exit_code, 0, first.output + repr(first.exception))
        moved = self.root / "moved_input"
        self.output.rename(moved)
        mappings = _read_protein_mapping(moved / "protein_mapping.json")
        for identity in mappings.values():
            for receipt in identity["sources"]:
                self.assertTrue(
                    Path(receipt["absolute_path"]).resolve().is_relative_to(
                        moved.resolve()
                    )
                )
        second = CliRunner().invoke(
            catalog_cli.main,
            self.arguments() + self.export_options(
                mapping=moved / "protein_mapping.json",
                evidence=moved / "author_filter_evidence.json",
            ),
        )
        self.assertEqual(second.exit_code, 0, second.output + repr(second.exception))
        self.assertEqual(
            (self.output / "claudio_input.csv").read_bytes(),
            (moved / "claudio_input.csv").read_bytes(),
        )