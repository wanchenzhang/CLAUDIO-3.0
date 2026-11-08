import copy
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from claudio.dataset_catalog import xifdr


class TestConstructCandidateReferencePolicy(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fasta = self.root / "search.fasta"
        self.fasta.write_text(
            ">construct_a OX=9913\nAKAAK\n", encoding="utf-8"
        )
        self.results = self.root / "results.csv"
        columns = [
            "Title", "Charge", "Precursor_Mass", "Peptide",
            "Peptide_Type", "Linker", "Modifications",
            "Evalue", "Score", "Proteins",
        ]
        with self.results.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerow({
                "Title": "run.100.100.3",
                "Charge": "3",
                "Precursor_Mass": "1000",
                "Peptide": "AKAAK(2)-AKAAK(5)",
                "Peptide_Type": "Cross-Linked",
                "Linker": "DSS",
                "Modifications": "null",
                "Evalue": "0.01",
                "Score": "1",
                "Proteins": "construct_a(2)-construct_a(5)/",
            })
        self.bundle = xifdr.read_plink_psms(
            self.results, self.fasta,
            dataset_accession="PXD000001",
        )
        self.reference_path = self.root / "synthetic_reference.json"
        self.reference_body = {
            "primaryAccession": "P02769",
            "organism": {"taxonId": 9913},
            "sequence": {"value": "MAKAAK"},
        }
        self.receipt = {
            "kind": "uniprot_api",
            "request_url": "https://rest.uniprot.org/uniprotkb/P02769.json",
            "response_url": "https://rest.uniprot.org/uniprotkb/P02769.json",
            "status": 200,
            "path": self.reference_path.name,
            "absolute_path": str(self.reference_path),
        }
        self.write_reference(self.reference_body)
        self.report = {
            "status": "candidates_for_review",
            "reason": "",
            "gene_hint": "POLICY_TEST",
            "taxon_id": 9913,
            "sources": [self.receipt],
            "identity_accepted": False,
            "candidates": [{
                "accession": "P02769",
                "sequence": "MAKAAK",
                "taxon_id": 9913,
            }],
        }
        self.mappings = {
            "construct_a": {
                "status": "unmapped",
                "reason": "source_construct_requires_reference_review",
                "sources": [],
            }
        }

    def write_reference(self, body):
        self.reference_path.write_text(
            json.dumps(body) + "\n", encoding="utf-8"
        )
        self.receipt["sha256"] = hashlib.sha256(
            self.reference_path.read_bytes()
        ).hexdigest()

    def export(self, name="input"):
        projection = xifdr.project_xifdr_candidates(
            self.bundle, self.mappings,
            {"construct_a": self.report},
            export_unverified_candidates=True,
        )
        return Path(xifdr.write_claudio_input(
            self.bundle, self.mappings, self.root / name,
            allowed_crosslinkers={"DSS"},
            allowed_residues={"K"},
            reference_projection=projection,
            export_unverified_candidates=True,
        ))

    def profile(self):
        return {
            "accession": "P02769",
            "taxon_id": 9913,
            "reference_sequence": "MAKAAK",
        }

    def test_explicit_candidate_mode_supports_undeclared_construct_accession(self):
        original = copy.deepcopy(self.bundle)
        try:
            destination = self.export()
        except ValueError as error:
            self.fail(f"Construct candidate export was rejected: {error}")

        with (destination / "claudio_input.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["position1"], "3")
        self.assertEqual(rows[0]["position2"], "6")
        self.assertEqual(rows[0]["entry1"], "P02769")
        self.assertEqual(rows[0]["entry2"], "P02769")
        self.assertEqual(rows[0]["input_scope"], "unverified_candidate")
        self.assertEqual(rows[0]["identity_scope"], "candidate_hypothesis")
        self.assertEqual(rows[0]["fdr_verification_status"], "unverified")

        manifest = json.loads(
            (destination / "manifest.json").read_text(encoding="utf-8")
        )
        scope = manifest["reference_projection_policy"].get("scope")
        self.assertTrue(scope)
        self.assertNotEqual(scope, "declared_source_accession")
        self.assertEqual(self.bundle, original)
        self.assertIsNone(
            self.bundle["proteins"]["construct_a"]["uniprot_accession"]
        )

    def test_default_evidence_policy_still_requires_declared_accession(self):
        with self.assertRaisesRegex(ValueError, "(?i)declared.*accession"):
            xifdr._reference_evidence(
                self.profile(), self.report,
                self.bundle["proteins"]["construct_a"],
            )

    def test_changed_reference_checksum_is_rejected_before_output(self):
        self.reference_path.write_bytes(
            self.reference_path.read_bytes() + b"\n"
        )
        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            self.export()
        self.assertFalse((self.root / "input").exists())

    def test_candidate_mode_still_checks_reference_sequence_and_organism(self):
        for index, field in enumerate(("sequence", "organism")):
            with self.subTest(field=field):
                body = copy.deepcopy(self.reference_body)
                if field == "sequence":
                    body["sequence"]["value"] = "MATAAK"
                else:
                    body["organism"]["taxonId"] = 9606
                self.write_reference(body)
                with self.assertRaisesRegex(
                    ValueError, "(?i)sequence|organism"
                ):
                    self.export(name=f"wrong_{index}")
                self.assertFalse((self.root / f"wrong_{index}").exists())
        self.write_reference(self.reference_body)

    def test_optional_declaration_does_not_bypass_accession_change_review(self):
        source = copy.deepcopy(self.bundle["proteins"]["construct_a"])
        source["uniprot_accession"] = "P02768"
        with self.assertRaisesRegex(ValueError, "(?i)review|accession"):
            xifdr._reference_evidence(
                self.profile(), self.report, source,
                require_declared_accession=False,
            )
