import csv
import hashlib
import io
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from click.testing import CliRunner

from claudio.dataset_catalog import kojak
from claudio.dataset_catalog.cli import main


FIXTURE = Path(__file__).resolve().parents[1] / "fixtures/kojak_pxd015275"


class TestKojak151Native(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        manifest = json.loads(
            (FIXTURE / "manifest.json").read_text(encoding="utf-8")
        )
        for name, expected in manifest["output_sha256"].items():
            observed = hashlib.sha256((FIXTURE / name).read_bytes()).hexdigest()
            if observed != expected:
                raise AssertionError(f"Fixture checksum changed: {name}")

        cls.expected = json.loads(
            (FIXTURE / "expected.json").read_text(encoding="utf-8")
        )
        cls.results = FIXTURE / "native.kojak.txt"
        cls.fasta = FIXTURE / "search.fasta"
        lines = cls.results.read_text(encoding="utf-8-sig").splitlines()
        cls.columns = next(csv.reader([lines[1]], delimiter="\t"))
        cls.source_rows = list(csv.DictReader(
            io.StringIO("\n".join(lines[1:])),
            delimiter="\t",
        ))
        cls._bundle = None

    def read_bundle(self):
        if self.__class__._bundle is None:
            try:
                self.__class__._bundle = kojak.read_kojak_csms(
                    self.results,
                    self.fasta,
                    dataset_accession="PXD015275",
                    decoy_prefix="decoy_",
                )
            except ValueError as error:
                self.fail(f"Cannot read reviewed Kojak 1.5.1 fixture: {error}")
        return self.__class__._bundle

    def test_import_preserves_every_source_record(self):
        bundle = self.read_bundle()
        records = bundle["records"]
        self.assertEqual(len(records), self.expected["source_records"])
        self.assertEqual(len(records), len(self.source_rows))
        self.assertEqual(bundle["native_version"], "1.5.1")
        self.assertEqual(bundle["source_preamble"], ["Kojak version 1.5.1"])
        self.assertEqual(
            bundle["parser"]["native_format_version"], "1.5.1"
        )
        self.assertEqual(bundle["parser"]["validated_paired_records"], 362)
        self.assertEqual(
            dict(Counter(record["record_kind"] for record in records)),
            self.expected["record_kinds"],
        )

        for number, (original, record) in enumerate(
            zip(self.source_rows, records), start=1
        ):
            self.assertEqual(record["source_record_number"], number)
            self.assertEqual(record["source_fields"], original)
            self.assertEqual(
                record["source_values"],
                [original[column] for column in self.columns],
            )
            self.assertIsNone(record["spectrum_file"])
            self.assertEqual(
                record["scan_number"], int(original["Scan Number"])
            )

        self.assertEqual(
            bundle["provenance"]["psm_sha256"],
            hashlib.sha256(self.results.read_bytes()).hexdigest(),
        )
        paired = [
            record for record in records
            if record["record_kind"] == "paired_crosslink_candidate"
        ]
        self.assertEqual(len(paired), 362)
        self.assertEqual(
            dict(Counter(record["linker_mass_raw"] for record in paired)),
            self.expected["paired_native_linker_masses"],
        )
        self.assertFalse(bundle["identity_mapping_applied"])
        self.assertFalse(bundle["input_exported"])
        self.assertFalse(bundle["assessment_performed"])
        self.assertEqual(
            bundle["fdr_policy"]["verification_status"], "unverified"
        )
        self.assertFalse(bundle["fdr_policy"]["per_record_fdr_applied"])

    def test_coordinates_match_independent_protein_associations(self):
        records = self.read_bundle()["records"]
        residues = Counter()
        modifications = Counter()
        candidate_fields = (
            "raw", "protein_id", "protein_position",
            "peptide_start", "peptide_end", "residue",
        )

        for reference in self.expected["paired_records"]:
            number = reference["source_record_number"]
            record = records[number - 1]
            self.assertEqual(record["record_status"], "valid")

            for side in ("a", "b"):
                endpoint = record["endpoints"][side]
                expected = reference["endpoints"][side]
                self.assertEqual(endpoint["status"], "valid")
                self.assertEqual(endpoint["peptide"], expected["peptide"])
                self.assertEqual(
                    endpoint["link_position"], expected["link_position"]
                )

                actual_candidates = []
                for candidate in endpoint["protein_candidates"]:
                    self.assertEqual(candidate["status"], "valid")
                    actual_candidates.append({
                        field: candidate[field]
                        for field in candidate_fields
                    })
                    self.assertEqual(
                        candidate["decoy_label_present"],
                        candidate["protein_id"].startswith("decoy_"),
                    )
                    self.assertEqual(
                        candidate["residue"],
                        endpoint["peptide"][endpoint["link_position"] - 1],
                    )
                    residues[candidate["residue"]] += 1

                self.assertEqual(
                    actual_candidates, expected["protein_candidates"]
                )
                for modification in endpoint["modifications"]:
                    modifications[modification["raw_token"]] += 1
                    self.assertFalse(
                        modification["chemical_identity_assigned"]
                    )

        self.assertEqual(
            dict(residues),
            self.expected["validated_association_residues"],
        )
        self.assertEqual(
            dict(modifications),
            self.expected["paired_modification_tokens"],
        )

    def test_source_decoy_classifications_are_preserved(self):
        classifications = Counter()

        for record in self.read_bundle()["records"]:
            if record["record_kind"] != "paired_crosslink_candidate":
                continue
            labels = [
                token.startswith("decoy_")
                for field in ("Protein #1", "Protein #2")
                for token in record["source_fields"][field].split(";")
                if token
            ]
            self.assertTrue(labels)
            if all(labels):
                expected = "decoy_only"
            elif any(labels):
                expected = "contains_target_and_decoy_labels"
            else:
                expected = "target_labels_only"
            self.assertEqual(record["decoy_status"], expected)
            classifications[expected] += 1

        self.assertGreater(
            classifications["decoy_only"]
            + classifications["contains_target_and_decoy_labels"],
            0,
        )
        self.assertFalse(
            self.read_bundle()["decoy_policy"][
                "target_decoy_validation_performed"
            ]
        )

    def test_unrelated_fasta_does_not_validate_paired_records(self):
        bundle = kojak.read_kojak_csms(
            self.results,
            FIXTURE / "unrelated.fasta",
            dataset_accession="PXD015275",
            decoy_prefix="decoy_",
        )
        self.assertEqual(len(bundle["records"]), 1769)
        self.assertEqual(bundle["parser"]["validated_paired_records"], 0)

        paired = [
            record for record in bundle["records"]
            if record["record_kind"] == "paired_crosslink_candidate"
        ]
        self.assertEqual(len(paired), 362)
        for record in paired:
            self.assertEqual(record["record_status"], "not_evaluable")
            for endpoint in record["endpoints"].values():
                self.assertEqual(endpoint["status"], "not_evaluable")
                self.assertTrue(endpoint["protein_candidates"])
                for candidate in endpoint["protein_candidates"]:
                    self.assertEqual(
                        candidate["reason"], "protein_not_in_search_fasta"
                    )

    def test_unreviewed_versions_remain_rejected(self):
        original = self.results.read_bytes()
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "unsupported.kojak.txt"
            for version in ("1.5.2", "9.9.9"):
                with self.subTest(version=version):
                    path.write_bytes(original.replace(
                        b"Kojak version 1.5.1",
                        f"Kojak version {version}".encode("ascii"),
                        1,
                    ))
                    with self.assertRaisesRegex(
                        ValueError, "not yet validated"
                    ):
                        kojak.read_kojak_csms(
                            path,
                            self.fasta,
                            dataset_accession="PXD015275",
                            decoy_prefix="decoy_",
                        )

    def test_changed_column_layout_remains_rejected(self):
        original = self.results.read_bytes()
        changed = original.replace(
            b"\tScore\t", b"\tUnexpectedScore\t", 1
        )
        self.assertNotEqual(changed, original)

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "changed.kojak.txt"
            path.write_bytes(changed)
            with self.assertRaisesRegex(ValueError, "column layout"):
                kojak.read_kojak_csms(
                    path,
                    self.fasta,
                    dataset_accession="PXD015275",
                    decoy_prefix="decoy_",
                )

    def test_cli_audit_matches_native_api(self):
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary) / "audit"
            result = CliRunner().invoke(main, [
                "standardize",
                "--results", str(self.results),
                "--fasta", str(self.fasta),
                "--dataset-accession", "PXD015275",
                "--decoy-prefix", "decoy_",
                "--audit-only",
                "--output", str(destination),
            ])
            self.assertEqual(
                result.exit_code, 0,
                result.output + repr(result.exception),
            )
            audit = json.loads(
                (destination / "import_report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(audit, self.read_bundle())
            self.assertFalse((destination / "claudio_input.csv").exists())
