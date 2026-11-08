import copy
import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from claudio.dataset_catalog import openms_mzidentml as adapter


FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures/openms_mzidentml_pxd014523"
)


class TestOpenMSLabeledInputExport(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.evidence = self.root / "evidence"
        shutil.copytree(FIXTURE / "evidence", self.evidence)
        self.expected = json.loads(
            (self.evidence / "expected_input.json").read_text(encoding="utf-8")
        )
        self.bundle = adapter.read_openms_labeled_mzidentml(
            FIXTURE / "OpenPepXL_BSA_PDH_CID.mzid",
            FIXTURE / "search.fasta",
            dataset_accession="PXD014523",
            decoy_prefix="decoy_",
        )
        self.mapping = json.loads(
            (self.evidence / "reference_mapping.json").read_text(
                encoding="utf-8"
            )
        )
        for report in self.mapping.values():
            for receipt in report["sources"]:
                receipt["absolute_path"] = str(
                    self.evidence / receipt["path"]
                )
            review = report["sequence_projection_review"]
            review["absolute_path"] = str(self.evidence / review["path"])

    def writer(self):
        writer = getattr(
            adapter, "write_openms_labeled_claudio_input", None
        )
        self.assertTrue(
            callable(writer), "Missing OpenMS labeled CLAUDIO input exporter"
        )
        return writer

    def export(self, name="input", **options):
        return Path(self.writer()(
            self.bundle,
            self.mapping,
            self.root / name,
            allowed_residues={"D", "E"},
            export_unverified_candidates=True,
            **options,
        ))

    def load_output(self, destination):
        with (destination / "claudio_input.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames
            rows = list(reader)
        audit = json.loads(
            (destination / "psm_audit.json").read_text(encoding="utf-8")
        )
        manifest = json.loads(
            (destination / "manifest.json").read_text(encoding="utf-8")
        )
        return columns, rows, audit, manifest

    def test_exports_one_candidate_per_source_crosslink_group(self):
        original = copy.deepcopy(self.bundle)
        destination = self.export()
        columns, rows, audit, manifest = self.load_output(destination)

        self.assertEqual(columns[:8], [
            "peptide1", "peptide2", "position1", "position2",
            "k_pos1", "k_pos2", "entry1", "entry2",
        ])
        self.assertEqual(len(rows), 32)
        self.assertEqual(len(audit), 70)
        self.assertEqual(len({r["xlms_source_record_id"] for r in rows}), 32)

        for row, expected in zip(rows, self.expected):
            self.assertEqual(
                row["xlms_source_record_id"], expected["source_record_id"]
            )
            for name in columns[:8]:
                self.assertEqual(row[name], str(expected[name]))
            self.assertEqual(row["input_scope"], "unverified_candidate")
            self.assertEqual(row["identity_scope"], "candidate_hypothesis")
            self.assertEqual(row["fdr_verification_status"], "unverified")
            self.assertEqual(
                json.loads(row["source_native_ids_json"]),
                expected["source_native_ids"],
            )
            self.assertEqual(
                json.loads(row["source_rank_values_json"]),
                expected["source_rank_values"],
            )
            self.assertEqual(
                json.loads(row["native_fdr_raw_values_json"]),
                expected["native_fdr"]["raw_values"],
            )
            for name in (
                "selected_fdr_level", "selected_fdr_value",
                "selected_fdr_cutoff",
            ):
                self.assertEqual(row[name], "")

        for source, saved in zip(original["records"], audit):
            for key, value in source.items():
                self.assertEqual(saved[key], value, key)
            if source["record_kind"] == "loop_link_group":
                self.assertEqual(saved["input_status"], "not_evaluable")
                self.assertIn(
                    "loop_link_requires_single_peptide_support",
                    saved["input_reasons"],
                )
            else:
                self.assertEqual(saved["input_status"], "candidate_ready")

        self.assertEqual(
            sum(
                record["native_fdr"]["status"] == "out_of_range"
                for record in audit
            ),
            7,
        )
        policy = manifest["fdr_policy"]
        self.assertEqual(policy["kind"], "unverified_candidates")
        self.assertIsNone(policy["level"])
        self.assertIsNone(policy["cutoff"])
        self.assertIs(policy["per_record_fdr_applied"], False)
        self.assertIs(manifest["assessment_performed"], False)
        self.assertEqual(self.bundle, original)

        for name, expected in manifest["output_sha256"].items():
            self.assertEqual(
                hashlib.sha256((destination / name).read_bytes()).hexdigest(),
                expected,
            )

    def test_default_export_requires_explicit_candidate_mode(self):
        destination = self.root / "strict"
        with self.assertRaisesRegex(ValueError, "(?i)candidate|unverified"):
            self.writer()(
                self.bundle, self.mapping, destination,
                allowed_residues={"D", "E"},
            )
        self.assertFalse(destination.exists())

    def test_changed_reference_response_is_rejected_before_output(self):
        receipt = next(iter(self.mapping.values()))["sources"][0]
        path = Path(receipt["absolute_path"])
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
            self.export()
        self.assertFalse((self.root / "input").exists())

    def test_changed_projection_review_is_rejected_before_output(self):
        review = next(iter(self.mapping.values()))["sequence_projection_review"]
        path = Path(review["absolute_path"])
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
            self.export()
        self.assertFalse((self.root / "input").exists())

    def test_changed_projection_offset_is_rejected_before_output(self):
        next(iter(self.mapping.values()))[
            "source_to_reference_position_offset"
        ] = 25
        with self.assertRaisesRegex(
            ValueError, "(?i)projection|offset|coordinate|reference"
        ):
            self.export()
        self.assertFalse((self.root / "input").exists())

    def test_changed_native_sources_are_rejected_before_output(self):
        for kind, name in (
            ("psm", "OpenPepXL_BSA_PDH_CID.mzid"),
            ("fasta", "search.fasta"),
        ):
            with self.subTest(source=kind):
                source_dir = self.root / kind
                source_dir.mkdir()
                results = source_dir / "native.mzid"
                fasta = source_dir / "search.fasta"
                shutil.copy2(
                    FIXTURE / "OpenPepXL_BSA_PDH_CID.mzid", results
                )
                shutil.copy2(FIXTURE / "search.fasta", fasta)
                self.bundle = adapter.read_openms_labeled_mzidentml(
                    results, fasta,
                    dataset_accession="PXD014523",
                    decoy_prefix="decoy_",
                )
                changed = results if kind == "psm" else fasta
                changed.write_bytes(changed.read_bytes() + b"\n")
                with self.assertRaisesRegex(
                    ValueError, "(?i)checksum|changed"
                ):
                    self.export(name="changed_" + kind)
                self.assertFalse(
                    (self.root / ("changed_" + kind)).exists()
                )



from click.testing import CliRunner as OpenMSExportRunner
from claudio.dataset_catalog.cli import main as openms_export_main


class TestOpenMSLabeledExportCLI(unittest.TestCase):
    def setUp(self):
        TestOpenMSLabeledInputExport.setUp(self)

    def arguments(self, destination):
        return [
            "standardize",
            "--results", str(FIXTURE / "OpenPepXL_BSA_PDH_CID.mzid"),
            "--fasta", str(FIXTURE / "search.fasta"),
            "--dataset-accession", "PXD014523",
            "--mzidentml-profile", "openms-labeled",
            "--output", str(destination),
        ]

    def selection(self):
        return [
            "--reference-projection",
            str(self.evidence / "reference_mapping.json"),
            "--residue", "D", "--residue", "E",
            "--export-unverified-candidates",
        ]

    def test_cli_exports_the_verified_native_candidates(self):
        destination = self.root / "cli"
        result = OpenMSExportRunner().invoke(
            openms_export_main,
            self.arguments(destination) + self.selection(),
        )
        self.assertEqual(result.exit_code, 0, result.output)
        api = Path(adapter.write_openms_labeled_claudio_input(
            self.bundle, self.mapping, self.root / "api",
            allowed_residues={"D", "E"},
            export_unverified_candidates=True,
        ))
        self.assertEqual(
            (destination / "claudio_input.csv").read_bytes(),
            (api / "claudio_input.csv").read_bytes(),
        )
        audit = json.loads(
            (destination / "psm_audit.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(audit), 70)
        self.assertIn("unverified", result.output.lower())

    def test_cli_requires_explicit_candidate_inputs(self):
        cases = [
            (
                "--reference-projection",
                ["--residue", "D", "--residue", "E",
                 "--export-unverified-candidates"],
            ),
            (
                "--residue",
                ["--reference-projection",
                 str(self.evidence / "reference_mapping.json"),
                 "--export-unverified-candidates"],
            ),
            (
                "--export-unverified-candidates",
                ["--reference-projection",
                 str(self.evidence / "reference_mapping.json"),
                 "--residue", "D", "--residue", "E"],
            ),
        ]
        for index, (missing, options) in enumerate(cases):
            with self.subTest(missing=missing):
                destination = self.root / f"missing_{index}"
                result = OpenMSExportRunner().invoke(
                    openms_export_main,
                    self.arguments(destination) + options,
                )
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(missing, result.output)
                self.assertFalse(destination.exists())

    def test_cli_rejects_confidence_and_linker_overrides(self):
        for index, (option, values) in enumerate([
            ("--fdr-level", ["psm"]),
            ("--crosslinker", ["PDH"]),
            ("--mapping", [str(self.evidence / "reference_mapping.json")]),
        ]):
            with self.subTest(option=option):
                destination = self.root / f"rejected_{index}"
                result = OpenMSExportRunner().invoke(
                    openms_export_main,
                    self.arguments(destination) + self.selection()
                    + [option] + values,
                )
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(option, result.output)
                self.assertFalse(destination.exists())

    def test_cli_rejects_changed_reference_evidence_before_output(self):
        receipt = next(iter(self.mapping.values()))["sources"][0]
        path = Path(receipt["absolute_path"])
        path.write_bytes(path.read_bytes() + b"\n")
        destination = self.root / "changed"
        result = OpenMSExportRunner().invoke(
            openms_export_main,
            self.arguments(destination) + self.selection(),
        )
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("checksum", result.output.lower())
        self.assertFalse(destination.exists())
