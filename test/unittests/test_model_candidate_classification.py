"""Classify resolved model candidates without turning violations into inter links."""
import unittest

import numpy as np
import pandas as pd

from claudio.module04.src.algorithm.combine_reevals import combine_inter_reevaluations, def_xl_type
from claudio.module04.src.algorithm.relationship_classification import classify_model_candidates
from claudio.module04.src.algorithm.distance_assessment import assess_candidate_distance
from claudio.utils.utils import clean_dataset, minimize_dataset


class TestModelCandidateClassification(unittest.TestCase):
    def setUp(self):
        self.row = {
            "source_record_id": "record_000001", "candidate_id": "candidate_000001",
            "experiment": "condition_A", "source_crosslink_type": "Homodimer",
            "unip_id_a": "P1", "unip_id_b": "P1",
            "seq_a": "KKK", "seq_b": "KKK", "pep_a": "K", "pep_b": "K",
            "pos_a": 1, "pos_b": 2, "res_pos_a": 1, "res_pos_b": 1,
            "position_status_a": "valid", "position_status_b": "corrected",
            "pdb_id": "1ABC", "pdb_method": "X-RAY DIFFRACTION",
            "pdb_resolution": 2.0, "path": "/archive/1ABC.cif",
            "chain_a": "A", "chain_b": "A", "pdb_pos_a": 1, "pdb_pos_b": 2,
            "pdb_hetflag_a": " ", "pdb_hetflag_b": " ",
            "pdb_icode_a": " ", "pdb_icode_b": " ",
            "res_crit_a": True, "res_crit_b": True,
            "pLDDT_a": "-", "pLDDT_b": "-",
            "eucl_dist": 15.0, "topo_dist": 20.0,
            "homo_adjacency": 0.1, "homo_int_overl": 0.0, "homo_pep_overl": False,
        }

    def combine(self, data, scoring=True):
        return combine_inter_reevaluations(
            data, 70.0, 5.0, 35.0, None, 50.0, 0.0, scoring,
        )

    def test_same_chain_distance_violation_does_not_become_inter(self):
        row = {**self.row, "eucl_dist": 55.0, "topo_dist": 60.0}
        result = self.combine(pd.DataFrame([row])).iloc[0]
        self.assertEqual(result["XL_type"], "intra")
        self.assertEqual(result["XL_type_reason"], "mapped_same_chain")
        self.assertEqual(result["legacy_XL_type"], "inter")
        self.assertIn("distance", result["legacy_evidence"])
        self.assertEqual(result["source_crosslink_type"], "Homodimer")

    def test_same_accession_cross_chain_is_inter_even_if_legacy_score_says_intra(self):
        result = self.combine(pd.DataFrame([
            {**self.row, "chain_b": "B"},
        ])).iloc[0]
        self.assertEqual(result["XL_type"], "inter")
        self.assertEqual(result["protein_relationship"], "same_accession")
        self.assertEqual(result["legacy_score_XL_type"], "intra")
        self.assertEqual(result["legacy_inter_score"], 0.0)
        self.assertNotIn("score_XL_type", result.index)
        self.assertNotIn("inter_score", result.index)
        self.assertNotIn("XL_confirmed", result.index)

    def test_scalar_type_helper_uses_the_same_candidate_relationship_rule(self):
        row = pd.Series({**self.row, "chain_b": "B", "evidence": ""})
        self.assertEqual(def_xl_type(row), "inter")

    def test_unresolved_positions_stay_undetermined_despite_numeric_distances(self):
        for site in ("a", "b"):
            for status in ("ambiguous", "invalid"):
                with self.subTest(site=site, status=status):
                    data = pd.DataFrame([{**self.row, f"position_status_{site}": status}])
                    result = self.combine(data).iloc[0]
                    self.assertEqual(result["XL_type"], "undetermined")
                    self.assertEqual(result["XL_type_reason"], "unresolved_position")
                    self.assertEqual(result["legacy_XL_type"], "intra")

    def test_missing_model_chain_or_position_status_does_not_default_to_intra(self):
        cases = (
            ("path", "-", "no_usable_structure"),
            ("chain_b", "-", "missing_chain_assignment"),
            ("position_status_a", None, "missing_position_status"),
            ("position_status_a", "valid?", "unrecognized_position_status"),
        )
        for field, value, reason in cases:
            with self.subTest(field=field):
                result = classify_model_candidates(pd.DataFrame([
                    {**self.row, field: value},
                ])).iloc[0]
                self.assertEqual(result["XL_type"], "undetermined")
                self.assertEqual(result["XL_type_reason"], reason)

    def test_failed_or_missing_mapping_confirmation_is_not_classified(self):
        cases = (
            ("res_crit_b", False, "unmapped_endpoint"),
            ("res_crit_b", "True", "unmapped_endpoint"),
            ("res_crit_b", 1, "unmapped_endpoint"),
            ("pdb_pos_b", None, "missing_residue_identity"),
            ("pdb_pos_b", 2.5, "invalid_residue_identity"),
        )
        for field, value, reason in cases:
            with self.subTest(field=field):
                result = classify_model_candidates(pd.DataFrame([
                    {**self.row, field: value},
                ])).iloc[0]
                self.assertEqual(result["XL_type"], "undetermined")
                self.assertEqual(result["XL_type_reason"], reason)

    def test_af_complex_requires_a_matched_model_and_explicit_mapping(self):
        base = {
            **self.row, "af_match_status": "matched", "af_model_id": "AF-model",
            "af_model_version": 1, "pdb_method": "ALPHAFOLD",
            "mapping_status_a": "mapped", "mapping_status_b": "mapped",
        }
        for changes, reason in (
            ({"af_match_status": "no_candidate_model"}, "no_candidate_model"),
            ({"af_model_version": None}, "missing_model_identity"),
            ({"mapping_status_b": None}, "missing_mapping_status"),
            ({"mapping_status_b": "not_evaluable"}, "unmapped_endpoint"),
        ):
            with self.subTest(changes=changes):
                result = classify_model_candidates(pd.DataFrame([{**base, **changes}])).iloc[0]
                self.assertEqual(result["XL_type"], "undetermined")
                self.assertEqual(result["XL_type_reason"], reason)

    def test_different_accessions_on_one_chain_are_marked_as_a_conflict(self):
        result = classify_model_candidates(pd.DataFrame([
            {**self.row, "unip_id_b": "P2"},
        ])).iloc[0]
        self.assertEqual(result["XL_type"], "undetermined")
        self.assertEqual(result["XL_type_reason"], "protein_chain_conflict")
        self.assertEqual(result["protein_relationship"], "different_accessions")

    def test_confidence_affects_assessment_without_rewriting_chain_topology(self):
        for confidence in (20.0, None, "unknown"):
            with self.subTest(confidence=confidence):
                data = pd.DataFrame([{
                    **self.row, "pdb_method": "ALPHAFOLD", "pdb_id": "afP1",
                    "pLDDT_a": confidence, "pLDDT_b": confidence,
                }])
                result = self.combine(data).iloc[0]
                self.assertEqual(result["XL_type"], "intra")
                assessed = assess_candidate_distance(
                    result, metric="topological", minimum=0.0, maximum=35.0,
                    boundary="inclusive", plddt_cutoff=70.0,
                )
                self.assertEqual(assessed["status"], "not_evaluable")

    def test_identical_atom_or_coincident_residue_endpoints_are_undetermined(self):
        for changes, reason in (
            ({"distance_reason": "identical_atom_endpoints"}, "identical_atom_endpoints"),
            ({"pdb_pos_b": 1}, "coincident_residue_endpoints"),
        ):
            with self.subTest(changes=changes):
                result = self.combine(pd.DataFrame([{**self.row, **changes}])).iloc[0]
                self.assertEqual(result["XL_type"], "undetermined")
                self.assertEqual(result["XL_type_reason"], reason)
                assessed = assess_candidate_distance(
                    result, metric="topological", minimum=0.0, maximum=35.0,
                    boundary="inclusive", plddt_cutoff=70.0,
                )
                self.assertEqual(assessed["status"], "not_evaluable")

    def test_author_residue_numbers_and_insertion_codes_are_preserved(self):
        for number in (0, -3, 104.0):
            with self.subTest(number=number):
                row = {**self.row, "pdb_pos_a": number, "pdb_pos_b": number,
                       "pdb_icode_a": "A", "pdb_icode_b": "B"}
                result = classify_model_candidates(pd.DataFrame([row])).iloc[0]
                self.assertEqual(result["XL_type"], "intra")
                self.assertEqual(result["pdb_pos_a"], number)
                self.assertEqual(result["pdb_icode_b"], "B")

    def test_preserves_reported_input_classification_and_all_source_rows(self):
        data = pd.DataFrame([
            {**self.row, "XL_type": "reported_inter"},
            {**self.row, "source_record_id": "record_000002",
             "experiment": "condition_B", "XL_type": "reported_intra"},
        ], index=["7", "12_1"])
        before = data.copy(deep=True)
        result = self.combine(data)
        self.assertEqual(result["input_XL_type"].tolist(), ["reported_inter", "reported_intra"])
        self.assertEqual(result["source_crosslink_type"].tolist(), ["Homodimer"] * 2)
        self.assertEqual(result["experiment"].tolist(), ["condition_A", "condition_B"])
        self.assertTrue(result["XL_type_scope"].eq("model_candidate").all())
        self.assertTrue(result["XL_type_rule_version"].eq("model-candidate-type-v1").all())
        pd.testing.assert_series_equal(result["source_record_id"], before["source_record_id"])
        pd.testing.assert_frame_equal(data, before)

    def test_modern_candidate_tables_are_not_reduced_by_legacy_cleanup(self):
        data = self.combine(pd.DataFrame([
            self.row, {**self.row, "source_record_id": "record_000002",
                       "experiment": "condition_B"},
        ], index=["7", "12_1"]))
        for cleanup in (clean_dataset, minimize_dataset):
            with self.subTest(cleanup=cleanup.__name__):
                events = []
                result = cleanup(data, selection_log=events)
                pd.testing.assert_frame_equal(result, data)
                self.assertEqual(events, [])

    def test_refreshes_stale_classification_and_handles_empty_tables(self):
        data = pd.DataFrame([self.row])
        data["XL_type"] = "inter"
        data["XL_type_scope"] = "stale"
        result = classify_model_candidates(data)
        self.assertEqual(result.iloc[0]["XL_type"], "intra")
        self.assertEqual(data.iloc[0]["XL_type"], "inter")
        empty = self.combine(data.iloc[:0])
        self.assertTrue(empty.empty)
        for column in ("XL_type", "XL_type_reason", "XL_type_scope",
                       "XL_type_rule_version", "legacy_XL_type", "legacy_evidence"):
            self.assertIn(column, empty)

    def test_legacy_inference_failure_does_not_discard_a_mapped_candidate(self):
        data = pd.DataFrame([{
            **self.row, "pdb_method": "ALPHAFOLD", "pdb_id": "afP1",
            "pLDDT_a": "unknown", "pLDDT_b": "unknown",
        }])
        result = self.combine(data).iloc[0]
        self.assertEqual(result["XL_type"], "intra")
        self.assertEqual(result["legacy_inference_status"], "unavailable")
        self.assertEqual(result["source_record_id"], "record_000001")


if __name__ == "__main__":
    unittest.main()
