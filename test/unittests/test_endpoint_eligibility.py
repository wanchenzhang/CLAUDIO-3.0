"""Check endpoint prerequisites before interpreting structural distances."""
import unittest

import numpy as np
import pandas as pd

from claudio.module04.src.algorithm.distance_assessment import (
    check_endpoint_eligibility,
)


class TestEndpointEligibility(unittest.TestCase):
    def setUp(self):
        self.row = pd.Series({
            "source_record_id": "record_000001", "experiment": "condition_A",
            "source_crosslink_type": "Homodimer", "XL_type": "inter",
            "position_status_a": "valid", "position_status_b": "corrected",
            "chain_a": "AA", "chain_b": "AB",
            "pdb_pos_a": 104, "pdb_pos_b": 108,
            "pdb_hetflag_a": " ", "pdb_icode_a": "A",
            "pdb_hetflag_b": " ", "pdb_icode_b": " ",
            "mapping_status_a": "mapped", "mapping_status_b": "mapped",
            "mapping_reason_a": None, "mapping_reason_b": None,
            "res_crit_a": True, "res_crit_b": True,
            "pLDDT_a": 90.0, "pLDDT_b": 85.0,
            "eucl_dist": 15.0, "topo_dist": 20.0,
        })

    def check(self, row=None, source="af_complex", cutoff=70.0):
        return check_endpoint_eligibility(
            self.row if row is None else row,
            structure_source=source, plddt_cutoff=cutoff,
        )

    def test_all_structure_types_accept_valid_mapped_endpoints(self):
        before = self.row.copy(deep=True)
        for source in ("pdb", "af_monomer", "af_complex"):
            with self.subTest(source=source):
                result = self.check(source=source)
                self.assertTrue(result["eligible"])
                self.assertEqual(result["reason"], "endpoints_evaluable")
                self.assertIsNone(result["site"])
        pd.testing.assert_series_equal(self.row, before)

    def test_unresolved_m1_positions_block_stale_numeric_distances(self):
        for site in ("a", "b"):
            for status in ("invalid", "ambiguous"):
                with self.subTest(site=site, status=status):
                    row = self.row.copy()
                    row[f"position_status_{site}"] = status
                    result = self.check(row)
                    self.assertFalse(result["eligible"])
                    self.assertEqual(result["reason"], "unresolved_position")
                    self.assertEqual(result["site"], site)
                    self.assertEqual(row["topo_dist"], 20.0)

    def test_missing_or_unrecognized_m1_status_is_not_assumed_valid(self):
        for site in ("a", "b"):
            for status, reason in (
                (None, "missing_position_status"),
                (pd.NA, "missing_position_status"),
                ("-", "missing_position_status"),
                ("valid?", "unrecognized_position_status"),
            ):
                with self.subTest(site=site, status=status):
                    row = self.row.copy()
                    row[f"position_status_{site}"] = status
                    result = self.check(row)
                    self.assertFalse(result["eligible"])
                    self.assertEqual(result["reason"], reason)
                    self.assertEqual(result["site"], site)

    def test_af_mapping_failures_keep_the_upstream_reason(self):
        row = self.row.copy()
        row["mapping_status_b"] = "not_evaluable"
        row["mapping_reason_b"] = "missing_atom"
        result = self.check(row)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "unmapped_endpoint")
        self.assertEqual(result["site"], "b")
        self.assertEqual(result["upstream_reason"], "missing_atom")

    def test_af_complex_requires_an_explicit_mapping_status(self):
        row = self.row.drop(["mapping_status_b", "mapping_reason_b"])
        result = self.check(row)
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "missing_mapping_status")
        self.assertEqual(result["site"], "b")

    def test_ordinary_mapping_uses_legacy_endpoint_confirmation(self):
        row = self.row.drop([
            "mapping_status_a", "mapping_status_b",
            "mapping_reason_a", "mapping_reason_b",
        ])
        for source in ("pdb", "af_monomer"):
            with self.subTest(source=source):
                self.assertTrue(self.check(row, source)["eligible"])
        row["res_crit_b"] = False
        result = self.check(row, "pdb")
        self.assertFalse(result["eligible"])
        self.assertEqual(result["reason"], "unmapped_endpoint")
        self.assertEqual(result["site"], "b")

    def test_mapping_confirmation_is_not_inferred_from_truthy_values(self):
        for confirmation in (None, pd.NA, "False", "True", 1):
            with self.subTest(confirmation=confirmation):
                row = self.row.copy()
                row["res_crit_a"] = confirmation
                result = self.check(row)
                self.assertFalse(result["eligible"])
                self.assertEqual(result["reason"], "unmapped_endpoint")
        row = self.row.copy()
        row["res_crit_a"] = np.bool_(True)
        self.assertTrue(self.check(row)["eligible"])

    def test_missing_chain_or_residue_identity_blocks_evaluation(self):
        cases = (
            ("chain_a", "-", "missing_chain_assignment", "a"),
            ("chain_b", None, "missing_chain_assignment", "b"),
            ("pdb_pos_a", None, "missing_residue_identity", "a"),
            ("pdb_pos_b", 108.5, "invalid_residue_identity", "b"),
            ("pdb_pos_b", True, "invalid_residue_identity", "b"),
        )
        for field, value, reason, site in cases:
            with self.subTest(field=field, value=value):
                row = self.row.copy()
                row[field] = value
                result = self.check(row)
                self.assertFalse(result["eligible"])
                self.assertEqual((result["reason"], result["site"]), (reason, site))

    def test_author_residue_numbers_are_not_treated_as_protein_positions(self):
        for number in (0, -3, 104.0):
            with self.subTest(number=number):
                row = self.row.copy()
                row["pdb_pos_a"] = number
                self.assertTrue(self.check(row)["eligible"])

    def test_af_confidence_at_the_cutoff_is_accepted_and_lower_values_are_blocked(self):
        for source in ("af_monomer", "af_complex"):
            for site in ("a", "b"):
                with self.subTest(source=source, site=site):
                    row = self.row.copy()
                    row[f"pLDDT_{site}"] = 70.0
                    self.assertTrue(self.check(row, source)["eligible"])
                    row[f"pLDDT_{site}"] = 69.999
                    result = self.check(row, source)
                    self.assertFalse(result["eligible"])
                    self.assertEqual((result["reason"], result["site"]), ("low_plddt", site))

    def test_af_confidence_accepts_numeric_text_from_mixed_csv_columns(self):
        row = self.row.copy()
        row["pLDDT_a"] = "90.0"
        row["pLDDT_b"] = "70"
        self.assertTrue(self.check(row)["eligible"])

    def test_missing_and_invalid_af_confidence_are_reported_separately(self):
        for value, reason in (
            (None, "missing_plddt"), (pd.NA, "missing_plddt"),
            ("-", "missing_plddt"), (float("nan"), "missing_plddt"),
            ("unknown", "invalid_plddt"), (float("inf"), "invalid_plddt"),
            (-1.0, "invalid_plddt"), (100.1, "invalid_plddt"),
            (True, "invalid_plddt"), (np.bool_(True), "invalid_plddt"),
        ):
            with self.subTest(value=value):
                row = self.row.copy()
                row["pLDDT_b"] = value
                result = self.check(row)
                self.assertFalse(result["eligible"])
                self.assertEqual((result["reason"], result["site"]), (reason, "b"))

    def test_pdb_confidence_fields_are_not_interpreted_as_af_plddt(self):
        for value in ("-", None, 0.0):
            with self.subTest(value=value):
                row = self.row.copy()
                row["pLDDT_a"] = row["pLDDT_b"] = value
                self.assertTrue(self.check(row, "pdb")["eligible"])

    def test_records_scope_cutoff_and_rule_version_without_reclassifying_the_source(self):
        result = self.check()
        self.assertEqual(result["structure_source"], "af_complex")
        self.assertEqual(result["plddt_cutoff"], 70.0)
        self.assertEqual(result["rule_version"], "endpoint-eligibility-v1")
        self.assertIsNone(result["upstream_reason"])
        self.assertEqual(self.row["source_crosslink_type"], "Homodimer")
        self.assertEqual(self.row["XL_type"], "inter")
        self.assertEqual(self.row["experiment"], "condition_A")

    def test_invalid_source_and_cutoff_configuration_are_rejected(self):
        for source in ("unknown", "AF", None):
            with self.subTest(source=source):
                with self.assertRaisesRegex(ValueError, "structure_source"):
                    self.check(source=source)
        for cutoff in (-1, 101, float("nan"), float("inf"), True, "70"):
            with self.subTest(cutoff=cutoff):
                with self.assertRaisesRegex(ValueError, "cutoff"):
                    self.check(cutoff=cutoff)


if __name__ == "__main__":
    unittest.main()
