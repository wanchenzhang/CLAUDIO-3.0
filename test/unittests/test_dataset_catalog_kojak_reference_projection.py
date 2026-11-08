import copy
import csv
import hashlib
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from claudio.dataset_catalog import kojak


FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/kojak_pxd015275"
EVIDENCE = FIXTURE / "evidence"


def load_mapping(directory):
    mapping = json.loads(
        (directory / "reference_mapping.json").read_text(encoding="utf-8")
    )

    def resolve(value):
        if isinstance(value, dict):
            if {"kind", "sha256", "path"} <= set(value):
                path = (directory / value["path"]).resolve()
                if not path.is_relative_to(directory.resolve()):
                    raise ValueError("Reference receipt leaves its evidence directory")
                value["absolute_path"] = str(path)
            for child in list(value.values()):
                resolve(child)
        elif isinstance(value, list):
            for child in value:
                resolve(child)

    resolve(mapping)
    return mapping


class TestKojakReviewedReferenceProjection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        manifest = json.loads(
            (EVIDENCE / "manifest.json").read_text(encoding="utf-8")
        )
        for name, expected in manifest["output_sha256"].items():
            observed = hashlib.sha256((EVIDENCE / name).read_bytes()).hexdigest()
            if observed != expected:
                raise AssertionError(f"Evidence checksum changed: {name}")

        cls.expected = json.loads(
            (EVIDENCE / "expected_input.json").read_text(encoding="utf-8")
        )
        cls.bundle = kojak.read_kojak_csms(
            FIXTURE / "native.kojak.txt",
            FIXTURE / "search.fasta",
            dataset_accession="PXD015275",
            decoy_prefix="decoy_",
        )

    def export(self, mapping, destination):
        return Path(kojak.write_kojak_claudio_input(
            self.bundle,
            mapping,
            destination,
            allowed_linker_masses={"110.0362", "116.0564"},
            allowed_residues={"K"},
            export_unverified_candidates=True,
        ))

    def test_exports_reviewed_intact_peptide_subset(self):
        original = copy.deepcopy(self.bundle)
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "input"
            try:
                self.export(load_mapping(EVIDENCE), destination)
            except ValueError as error:
                self.fail(f"Reviewed intact-peptide export is unsupported: {error}")

            with (destination / "claudio_input.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle))
            audit = json.loads(
                (destination / "psm_audit.json").read_text(encoding="utf-8")
            )
            manifest = json.loads(
                (destination / "manifest.json").read_text(encoding="utf-8")
            )

            self.assertEqual(len(rows), 56)
            self.assertEqual(len(audit), 1769)
            self.assertEqual(
                [row["xlms_source_record_id"] for row in rows],
                [
                    row["xlms_source_record_id"]
                    for row in self.expected["rows"]
                ],
            )
            for actual, expected in zip(rows, self.expected["rows"]):
                for field, value in expected.items():
                    self.assertEqual(actual[field], str(value))
                self.assertEqual(actual["input_scope"], "unverified_candidate")
                self.assertEqual(actual["identity_scope"], "candidate_hypothesis")
                self.assertEqual(actual["fdr_verification_status"], "unverified")
                for field in (
                    "selected_fdr_level", "selected_fdr_value",
                    "selected_fdr_cutoff",
                ):
                    self.assertEqual(actual[field], "")

            for original_record, audited in zip(self.bundle["records"], audit):
                for field, value in original_record.items():
                    self.assertEqual(audited[field], value)

            self.assertEqual(
                sum(record["input_status"] == "candidate_ready" for record in audit),
                56,
            )
            self.assertEqual(
                sum(record["input_status"] == "not_evaluable" for record in audit),
                1713,
            )
            by_id = {record["source_record_id"]: record for record in audit}
            for identifier in self.expected["sequence_mismatch_record_ids"]:
                self.assertEqual(by_id[identifier]["input_status"], "not_evaluable")
                self.assertTrue(any(
                    reason.startswith("reference_projection_unresolved_")
                    for reason in by_id[identifier]["input_reasons"]
                ))

            policy = manifest["reference_projection_policy"]
            self.assertEqual(policy["kind"], "reviewed_intact_peptide_projection")
            self.assertIs(policy["require_all_optimal_alignments"], True)
            self.assertIs(policy["identity_accepted"], False)
            self.assertEqual(manifest["fdr_policy"]["kind"], "unverified_candidates")
            self.assertIs(manifest["fdr_policy"]["per_record_fdr_applied"], False)
            self.assertIs(manifest["assessment_performed"], False)

        self.assertEqual(self.bundle, original)

    def test_changed_reference_body_is_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            evidence = root / "evidence"
            shutil.copytree(EVIDENCE, evidence)
            mapping = load_mapping(evidence)
            receipt = mapping["references"]["Im7_7-45"]["sources"][0]
            path = Path(receipt["absolute_path"])
            path.write_bytes(path.read_bytes() + b" ")
            destination = root / "input"
            with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
                self.export(mapping, destination)
            self.assertFalse(destination.exists())

    def test_changed_review_document_is_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            evidence = root / "evidence"
            shutil.copytree(EVIDENCE, evidence)
            mapping = load_mapping(evidence)
            review = mapping["references"]["Im7_7-45"]["sequence_projection_review"]
            path = Path(review["sources"][0]["absolute_path"])
            path.write_bytes(path.read_bytes() + b" ")
            destination = root / "input"
            with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
                self.export(mapping, destination)
            self.assertFalse(destination.exists())

    def test_review_cannot_override_recomputed_coordinate_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            evidence = root / "evidence"
            shutil.copytree(EVIDENCE, evidence)
            mapping = load_mapping(evidence)
            profile = mapping["references"]["SpyL32P"]
            profile["coordinate_evidence"]["mapping_option_count"] += 1
            review = profile["sequence_projection_review"]
            review["coordinate_evidence"] = copy.deepcopy(
                profile["coordinate_evidence"]
            )
            declaration = {
                key: value for key, value in review.items() if key != "sources"
            }
            receipt = review["sources"][0]
            path = Path(receipt["absolute_path"])
            path.write_text(
                json.dumps(declaration, ensure_ascii=False, indent=2) + "\n",
                encoding="utf-8",
            )
            receipt["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            destination = root / "input"
            with self.assertRaisesRegex(
                ValueError, "(?i)coordinate|projection|review"
            ):
                self.export(mapping, destination)
            self.assertFalse(destination.exists())

    def test_all_alignment_requirement_cannot_be_disabled(self):
        mapping = load_mapping(EVIDENCE)
        mapping["references"]["Im7_7-45"]["require_all_optimal_alignments"] = False
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "input"
            with self.assertRaises(ValueError):
                self.export(mapping, destination)
            self.assertFalse(destination.exists())

    def test_unknown_projection_rule_is_rejected(self):
        mapping = load_mapping(EVIDENCE)
        mapping["projection_rule_version"] = "unreviewed-rule"
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "input"
            with self.assertRaisesRegex(
                ValueError, "(?i)projection|unsupported|rule"
            ):
                self.export(mapping, destination)
            self.assertFalse(destination.exists())
