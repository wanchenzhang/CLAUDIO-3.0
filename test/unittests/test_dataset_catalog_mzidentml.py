import copy
import csv
import gzip
import importlib
import importlib.util
import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

from click.testing import CliRunner
from claudio.dataset_catalog import cli
from claudio.dataset_catalog.cli import _read_reference_projection


def tag(element):
    return element.tag.rsplit("}", 1)[-1]


class TestMzIdentMLStandardize(unittest.TestCase):
    def setUp(self):
        self.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures" / "mzidentml_pxd046392"
        )
        self.expected = json.loads(
            (self.fixture / "expected.json").read_text(encoding="utf-8")
        )
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def api(self):
        self.assertIsNotNone(
            importlib.util.find_spec("claudio.dataset_catalog.mzidentml"),
            "Missing mzIdentML adapter",
        )
        return importlib.import_module("claudio.dataset_catalog.mzidentml")

    def read_bundle(self, results=None, backend="pyxlms"):
        return self.api().read_mzidentml_csms(
            self.fixture / self.expected["native_filename"]
            if results is None else results,
            self.fixture / "search.fasta",
            dataset_accession="PXD046392",
            backend=backend,
        )

    def mapping(self):
        return _read_reference_projection(
            self.fixture / "evidence" / "reference_mapping.json"
        )

    def export(self, bundle=None, mapping=None, name="input", **options):
        arguments = {
            "allowed_crosslinkers": {"DSS"},
            "allowed_residues": {"K"},
            "export_unverified_candidates": True,
        }
        arguments.update(options)
        return Path(self.api().write_mzidentml_claudio_input(
            self.read_bundle() if bundle is None else bundle,
            self.mapping() if mapping is None else mapping,
            self.root / name,
            **arguments,
        ))

    def load_output(self, directory):
        with (directory / "claudio_input.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames
            rows = list(reader)
        audit = json.loads(
            (directory / "psm_audit.json").read_text(encoding="utf-8")
        )
        manifest = json.loads(
            (directory / "manifest.json").read_text(encoding="utf-8")
        )
        return columns, rows, audit, manifest

    def mutated_xml(self, transform):
        root = ET.parse(self.fixture / "native.mzid").getroot()
        first_sir = next(
            element for element in root.iter()
            if tag(element) == "SpectrumIdentificationResult"
            and any(
                child.get("accession") == "MS:1002511"
                for item in element
                if tag(item) == "SpectrumIdentificationItem"
                for child in item
            )
        )
        items = [
            item for item in first_sir
            if tag(item) == "SpectrumIdentificationItem"
        ]
        transform(root, items)
        destination = self.root / "changed.mzid.gz"
        destination.write_bytes(gzip.compress(
            ET.tostring(root, encoding="utf-8"), mtime=0
        ))
        return destination

    def test_import_preserves_all_identifications_and_unknown_decoys(self):
        bundle = self.read_bundle()
        self.assertEqual(bundle["format"], "mzidentml_csm")
        self.assertEqual(len(bundle["records"]), 17)
        self.assertEqual(bundle["parser"]["parsed_records"], 17)
        self.assertEqual(len(bundle["source_identifications"]), 1212)
        self.assertEqual(
            sum(len(item["items"]) for item in bundle["source_identifications"]),
            1249,
        )
        self.assertEqual(
            Counter(record["parser_orientation"] for record in bundle["records"]),
            {"same": 8, "swapped": 9},
        )
        for record, expected in zip(
            bundle["records"], self.expected["groups"]
        ):
            self.assertEqual(record["source_record_id"], expected["source_record_id"])
            self.assertEqual(
                record["source_sir_attributes"], expected["sir_attributes"]
            )
            self.assertEqual(record["source_items"], expected["items"])
            self.assertEqual(record["decoy_status"], "unknown")
            self.assertEqual(record["record_status"], "valid")
            for endpoint in record["endpoints"].values():
                self.assertIsNone(endpoint["decoy"])

    def test_broken_donor_receiver_pair_remains_in_audit(self):
        def change(root, items):
            peptide_id = items[0].get("peptide_ref")
            peptide = next(
                element for element in root.iter()
                if tag(element) == "Peptide"
                and element.get("id") == peptide_id
            )
            marker = next(
                element for element in peptide.iter()
                if element.get("accession") in ("MS:1002509", "MS:1002510")
            )
            marker.set("value", "9999999")
        bundle = self.read_bundle(
            self.mutated_xml(change), backend="native"
        )
        self.assertEqual(len(bundle["records"]), 17)
        self.assertEqual(bundle["records"][0]["record_status"], "not_evaluable")
        self.assertTrue(bundle["records"][0].get("source_error"))

    def test_export_matches_verified_coordinates_and_retains_unknown_confidence(self):
        bundle = self.read_bundle()
        original = copy.deepcopy(bundle)
        directory = self.export(bundle=bundle)
        columns, rows, audit, manifest = self.load_output(directory)
        self.assertEqual(columns[:8], [
            "peptide1", "peptide2", "position1", "position2",
            "k_pos1", "k_pos2", "entry1", "entry2",
        ])
        self.assertEqual(len(rows), 16)
        self.assertEqual(len(audit), 17)
        for row, expected in zip(rows, self.expected["rows"]):
            for name, value in expected.items():
                self.assertEqual(row[name], value)
            self.assertEqual(row["source_crosslinker"], "DSS")
            self.assertEqual(row["source_decoy_status"], "unknown")
            self.assertEqual(row["input_scope"], "unverified_candidate")
            self.assertEqual(row["fdr_verification_status"], "unverified")
            for name in (
                "selected_fdr_level", "selected_fdr_value", "selected_fdr_cutoff"
            ):
                self.assertEqual(row[name], "")
        for record, exported in zip(bundle["records"], audit):
            self.assertEqual(record["source_items"], exported["source_items"])
            self.assertEqual(record["decoy_status"], exported["decoy_status"])
        identifications = json.loads(
            (directory / "identification_audit.json").read_text(encoding="utf-8")
        )
        self.assertEqual(len(identifications), 1212)
        self.assertEqual(
            Counter(row["identity_scope"] for row in rows),
            {"exact_sequence": 8, "candidate_hypothesis": 8},
        )
        self.assertEqual(manifest["kind"], "mzidentml_claudio_input")
        self.assertFalse(manifest["fdr_policy"]["per_record_fdr_applied"])
        self.assertFalse(manifest["assessment_performed"])
        self.assertEqual(bundle, original)

    def test_default_export_requires_explicit_candidate_mode(self):
        with self.assertRaisesRegex(ValueError, "(?i)unverified"):
            self.export(export_unverified_candidates=False)
        self.assertFalse((self.root / "input").exists())

    def test_explicit_native_decoy_is_excluded_without_target_promotion(self):
        def change(root, items):
            reference = next(
                child for child in items[0]
                if tag(child) == "PeptideEvidenceRef"
            ).get("peptideEvidence_ref")
            evidence = next(
                element for element in root.iter()
                if tag(element) == "PeptideEvidence"
                and element.get("id") == reference
            )
            evidence.set("isDecoy", "true")
        bundle = self.read_bundle(self.mutated_xml(change), backend="native")
        first = bundle["records"][0]
        directory = self.export(bundle=bundle)
        _, rows, audit, _ = self.load_output(directory)
        self.assertEqual(len(rows), 15)
        self.assertNotIn(
            first["source_record_id"],
            {row["xlms_source_record_id"] for row in rows},
        )
        self.assertIn("decoy_record", audit[0]["input_reasons"])

    def test_explicit_threshold_failure_is_retained_and_excluded(self):
        def change(root, items):
            items[0].set("passThreshold", "false")
        bundle = self.read_bundle(self.mutated_xml(change), backend="native")
        directory = self.export(bundle=bundle)
        _, rows, audit, _ = self.load_output(directory)
        self.assertEqual(len(rows), 15)
        self.assertEqual(
            audit[0]["source_items"][0]["attributes"]["passThreshold"], "false"
        )
        self.assertIn("native_pass_threshold_failed", audit[0]["input_reasons"])

    def test_missing_candidate_reviews_are_rejected_before_output(self):
        cases = [
            ("sp|P0AEQ3|GLNH_ECOLI", "sequence_projection_review"),
            ("tr|A0A024R5I4|A0A024R5I4_HUMAN", "reference_accession_review"),
        ]
        for number, (identifier, key) in enumerate(cases):
            with self.subTest(identifier=identifier):
                mapping = self.mapping()
                del mapping["references"][identifier][key]
                name = f"missing_review_{number}"
                with self.assertRaisesRegex(ValueError, "(?i)review"):
                    self.export(mapping=mapping, name=name)
                self.assertFalse((self.root / name).exists())

    def test_changed_review_evidence_is_rejected_before_output(self):
        mapping = self.mapping()
        review = mapping["references"][
            "tr|A0A024R5I4|A0A024R5I4_HUMAN"
        ]["reference_accession_review"]
        receipt = review["sources"][0]
        changed = self.root / "changed_review.json"
        changed.write_bytes(Path(receipt["absolute_path"]).read_bytes() + b"\n")
        receipt["absolute_path"] = str(changed)
        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            self.export(mapping=mapping)
        self.assertFalse((self.root / "input").exists())

    def cli_arguments(self, output, projection=None):
        return [
            "standardize",
            "--results", str(self.fixture / self.expected["native_filename"]),
            "--fasta", str(self.fixture / "search.fasta"),
            "--dataset-accession", "PXD046392",
            "--reference-projection", str(
                self.fixture / "evidence" / "reference_mapping.json"
                if projection is None else projection
            ),
            "--crosslinker", "DSS",
            "--residue", "K",
            "--export-unverified-candidates",
            "--output", str(output),
        ]

    def test_cli_reproduces_verified_api_input(self):
        api_output = self.export(name="api_input")
        output = self.root / "cli_input"
        result = CliRunner().invoke(cli.main, self.cli_arguments(output))
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertEqual(
            (output / "claudio_input.csv").read_bytes(),
            (api_output / "claudio_input.csv").read_bytes(),
        )

    def test_cli_reuses_moved_reference_evidence(self):
        directory = self.export(name="api_input")
        moved = self.root / "moved"
        directory.rename(moved)
        output = self.root / "cli_input"
        result = CliRunner().invoke(
            cli.main,
            self.cli_arguments(output, moved / "reference_projection.json"),
        )
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertEqual(
            (output / "claudio_input.csv").read_bytes(),
            (moved / "claudio_input.csv").read_bytes(),
        )

class TestMzIdentMLResidueSelection(unittest.TestCase):
    def test_q_endpoint_is_preserved_and_requires_explicit_residue_selection(self):
        import csv
        import json
        import tempfile
        from collections import Counter
        from pathlib import Path

        from claudio.dataset_catalog import mzidentml
        from claudio.dataset_catalog.cli import _read_reference_projection

        fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures/mzidentml_pxd046392"
        )
        expected = json.loads(
            (fixture / "expected.json").read_text(encoding="utf-8")
        )
        bundle = mzidentml.read_mzidentml_csms(
            fixture / expected["native_filename"],
            fixture / "search.fasta",
            dataset_accession="PXD046392",
            backend="native",
        )
        mapping = _read_reference_projection(
            fixture / "evidence/reference_mapping.json"
        )
        matches = [
            record for record in bundle["records"]
            if record["source_sir_attributes"]["id"] == "SIR_148871_XLX47"
        ]
        self.assertEqual(len(matches), 1)
        record = matches[0]
        identifier = record["source_record_id"]
        self.assertEqual(record["endpoints"]["b"]["residue"], "Q")
        self.assertEqual(record["endpoints"]["b"]["peptide"], "QCLEMTTK")
        self.assertEqual(record["endpoints"]["b"]["link_position"], 1)
        self.assertEqual(record["endpoints"]["b"]["kind"], "residue")

        with tempfile.TemporaryDirectory() as temporary:
            for name, residues, count, exact_count in (
                ("kk", {"K"}, 16, 8),
                ("explicit_kq", {"K", "Q"}, 17, 9),
            ):
                with self.subTest(selection=name):
                    destination = Path(
                        mzidentml.write_mzidentml_claudio_input(
                            bundle,
                            mapping,
                            Path(temporary) / name,
                            allowed_crosslinkers={"DSS"},
                            allowed_residues=residues,
                            export_unverified_candidates=True,
                        )
                    )
                    with (destination / "claudio_input.csv").open(
                        newline="", encoding="utf-8"
                    ) as handle:
                        rows = list(csv.DictReader(handle))
                    audit = json.loads(
                        (destination / "psm_audit.json").read_text(
                            encoding="utf-8"
                        )
                    )
                    self.assertEqual(len(rows), count)
                    self.assertEqual(len(audit), 17)
                    self.assertEqual(
                        Counter(row["identity_scope"] for row in rows),
                        {
                            "exact_sequence": exact_count,
                            "candidate_hypothesis": 8,
                        },
                    )
                    retained = next(
                        item for item in audit
                        if item["source_record_id"] == identifier
                    )
                    self.assertEqual(
                        retained["source_items"], record["source_items"]
                    )
                    self.assertEqual(
                        retained["endpoints"]["b"]["residue"], "Q"
                    )
                    matching_rows = [
                        row for row in rows
                        if row["xlms_source_record_id"] == identifier
                    ]
                    if name == "kk":
                        self.assertEqual(matching_rows, [])
                        self.assertEqual(
                            retained["input_reasons"],
                            ["residue_not_selected_b"],
                        )
                    else:
                        self.assertEqual(len(matching_rows), 1)
                        self.assertEqual(
                            matching_rows[0]["peptide2"], "QCLEMTTK"
                        )
                        self.assertEqual(matching_rows[0]["k_pos2"], "1")
                        self.assertEqual(
                            matching_rows[0]["fdr_verification_status"],
                            "unverified",
                        )
