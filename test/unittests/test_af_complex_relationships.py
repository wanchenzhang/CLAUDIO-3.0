"""Describe AF candidate chains without changing the reported cross-link type."""
import unittest

import pandas as pd

from claudio.module04.src.algorithm.af_complex_assessment import (
    annotate_complex_relationships,
)


class TestAFComplexRelationships(unittest.TestCase):
    def setUp(self):
        self.row = {
            "inventory_record_id": "dataset:sheet:2",
            "source_record_id": "record_000001",
            "candidate_id": "candidate_000001",
            "source_crosslink_type": "Intra",
            "biological_replicate": "Replicate 1",
            "organelle_fraction": "Nucleus",
            "unip_id_a": "P1", "unip_id_b": "P1",
            "pos_a": 10, "pos_b": 20,
            "original_pos_a": 10, "original_pos_b": 20,
            "position_status_a": "valid", "position_status_b": "valid",
            "af_match_status": "matched",
            "af_model_id": "AF-0000000000000001",
            "af_model_version": 1,
            "path": "/archived/model.cif",
            "chain_a": "A", "chain_b": "A",
            "distance_status": "computed", "distance_reason": None,
            "eucl_dist": 15.0, "topo_dist": 20.0,
        }

    def test_chain_relation_is_separate_from_reported_classification(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "source_crosslink_type": "Homodimer", "chain_b": "B"},
            {**self.row, "source_crosslink_type": "Inter",
             "unip_id_b": "P2", "chain_b": "B"},
        ], index=[7, "7_1", 12])
        before = data.copy(deep=True)
        result = annotate_complex_relationships(data)

        pd.testing.assert_frame_equal(result[data.columns], before)
        pd.testing.assert_frame_equal(data, before)
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["same_chain", "cross_chain", "cross_chain"])
        self.assertEqual(result["model_chain_relation_reason"].tolist(),
                         ["same_model_chain", "different_model_chains",
                          "different_model_chains"])
        self.assertEqual(result["chain_relation_rule_version"].tolist(),
                         ["af-chain-v1"] * 3)

    def test_multicharacter_chain_ids_are_compared_in_full(self):
        data = pd.DataFrame([
            {**self.row, "chain_a": "AA", "chain_b": "AB"},
            {**self.row, "chain_a": "AA", "chain_b": "AA"},
        ])
        result = annotate_complex_relationships(data)
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["cross_chain", "same_chain"])
        self.assertEqual(result["chain_a"].tolist(), ["AA", "AA"])
        self.assertEqual(result["chain_b"].tolist(), ["AB", "AA"])

    def test_unusable_models_do_not_turn_stale_chain_labels_into_a_relation(self):
        cases = [
            ({"af_match_status": "no_candidate_model", "chain_b": "B"},
             "no_candidate_model"),
            ({"af_match_status": "no_matching_chain_pair"},
             "no_matching_chain_pair"),
            ({"af_match_status": pd.NA}, "missing_af_match_status"),
            ({"af_match_status": "unexpected"}, "unusable_af_match_status"),
            ({"af_model_id": None}, "missing_model_identity"),
            ({"af_model_version": pd.NA}, "missing_model_identity"),
            ({"path": "-"}, "no_usable_structure"),
        ]
        for changes, reason in cases:
            with self.subTest(changes=changes):
                result = annotate_complex_relationships(
                    pd.DataFrame([{**self.row, **changes}])
                )
                self.assertEqual(result.iloc[0]["model_chain_relation"], "unknown")
                self.assertEqual(result.iloc[0]["model_chain_relation_reason"], reason)
                self.assertEqual(result.iloc[0]["source_crosslink_type"], "Intra")

    def test_missing_or_placeholder_chains_are_unknown(self):
        for value in (None, pd.NA, float("nan"), "", "-"):
            for site in ("a", "b"):
                with self.subTest(chain=value, site=site):
                    data = pd.DataFrame([{**self.row, f"chain_{site}": value}])
                    result = annotate_complex_relationships(data)
                    self.assertEqual(result.iloc[0]["model_chain_relation"], "unknown")
                    self.assertEqual(result.iloc[0]["model_chain_relation_reason"],
                                     "missing_chain_assignment")

    def test_distance_exclusion_does_not_erase_a_known_chain_assignment(self):
        data = pd.DataFrame([
            {**self.row, "position_status_a": "ambiguous",
             "distance_status": "not_evaluable", "distance_reason": "unresolved_position",
             "eucl_dist": float("nan"), "topo_dist": float("nan")},
            {**self.row, "chain_b": "B", "distance_status": "not_evaluable",
             "distance_reason": "low_plddt", "eucl_dist": float("nan"),
             "topo_dist": float("nan")},
        ])
        result = annotate_complex_relationships(data)
        pd.testing.assert_frame_equal(result[data.columns], data)
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["same_chain", "cross_chain"])
        self.assertEqual(result["distance_status"].tolist(),
                         ["not_evaluable", "not_evaluable"])

    def test_missing_source_classification_is_not_invented(self):
        data = pd.DataFrame([self.row]).drop(columns=["source_crosslink_type"])
        result = annotate_complex_relationships(data)
        self.assertNotIn("source_crosslink_type", result.columns)
        self.assertNotIn("XL_type", result.columns)
        self.assertEqual(result.iloc[0]["model_chain_relation"], "same_chain")

    def test_reannotation_preserves_rows_and_refreshes_only_derived_fields(self):
        data = pd.DataFrame([
            {**self.row, "chain_b": "B"},
            self.row,
        ], index=[7, 7])
        data["model_chain_relation"] = "stale"
        data["model_chain_relation_reason"] = "stale"
        data["chain_relation_rule_version"] = "old"
        before = data.copy(deep=True)
        result = annotate_complex_relationships(data)
        self.assertEqual(result.index.tolist(), [7, 7])
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["cross_chain", "same_chain"])
        columns = list(self.row)
        pd.testing.assert_frame_equal(result[columns], before[columns])
        pd.testing.assert_frame_equal(data, before)

    def test_empty_input_keeps_headers(self):
        data = pd.DataFrame([self.row]).iloc[:0]
        result = annotate_complex_relationships(data)
        self.assertTrue(result.empty)
        for column in ("model_chain_relation", "model_chain_relation_reason",
                       "chain_relation_rule_version"):
            self.assertIn(column, result.columns)
        pd.testing.assert_frame_equal(result[data.columns], data)


if __name__ == "__main__":
    unittest.main()
