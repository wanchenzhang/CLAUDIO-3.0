"""Keep per-method assessments aligned with all candidate source rows."""
import copy
import json
import unittest

import pandas as pd

from claudio.module04.src.algorithm.distance_assessment import (
    annotate_distance_assessments,
)
from claudio.module04.src.algorithm.candidate_audit import build_candidate_audit


class TestDistanceAssessmentTable(unittest.TestCase):
    def setUp(self):
        self.row = {
            "source_record_id": "record_000001",
            "candidate_id": "candidate_000001",
            "experiment": "condition_A", "source_crosslink_type": "Homodimer",
            "XL_type": "inter", "score_XL_type": "intra",
            "af_match_status": "matched", "af_model_id": "AF-model",
            "af_model_version": 1, "path": "/archive/model.cif",
            "pdb_id": "AF-model", "pdb_method": "ALPHAFOLD",
            "chain_a": "AA", "chain_b": "AB",
            "position_status_a": "valid", "position_status_b": "corrected",
            "mapping_status_a": "mapped", "mapping_status_b": "mapped",
            "res_crit_a": True, "res_crit_b": True,
            "pdb_pos_a": 104, "pdb_pos_b": 108,
            "pLDDT_a": 90.0, "pLDDT_b": 90.0,
            "distance_status": "computed", "distance_reason": None,
            "eucl_dist": 15.0, "topo_dist": 20.0,
            "eucl_dist_lower_bound": None, "topo_dist_lower_bound": None,
            "eucl_dist_raw": "15.000", "topo_dist_raw": "20.000",
            "operator_note": "keep this source",
        }
        self.euclidean = {"minimum": 5.0, "maximum": 30.0, "boundary": "inclusive"}
        self.topological = {"minimum": 5.0, "maximum": 35.0, "boundary": "exclusive"}

    def annotate(self, data, **changes):
        options = {
            "euclidean_policy": self.euclidean,
            "topological_policy": self.topological,
            "plddt_cutoff": 70.0,
        }
        options.update(changes)
        return annotate_distance_assessments(data, **options)

    def test_adds_separate_results_without_pooling_candidates_or_sources(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "topo_dist": 45.0, "topo_dist_raw": "45.000"},
            {**self.row, "distance_status": "lower_bound",
             "topo_dist": None, "topo_dist_lower_bound": 75.328,
             "topo_dist_raw": ">75.328", "af_model_id": "AF-second"},
            {**self.row, "distance_status": "lower_bound",
             "topo_dist": None, "topo_dist_lower_bound": 20.0,
             "topo_dist_raw": ">20.000"},
            {**self.row, "position_status_b": "ambiguous"},
            {**self.row, "source_record_id": "record_000002",
             "experiment": "condition_B", "distance_status": "not_evaluable",
             "distance_reason": "identical_atom_endpoints"},
            {**self.row, "source_record_id": "record_000003",
             "experiment": "condition_C", "af_match_status": "no_candidate_model",
             "af_model_id": None, "path": "-", "chain_a": "-", "chain_b": "-"},
        ], index=["7", "7_1", "7_2", "7_3", "12", "19", "23"])
        before = data.copy(deep=True)
        result = self.annotate(data)
        self.assertEqual(result["topo_assessment_status"].tolist(), [
            "satisfied", "violated", "violated", "indeterminate",
            "not_evaluable", "not_evaluable", "not_evaluable",
        ])
        self.assertEqual(result["eucl_assessment_status"].tolist(), [
            "satisfied", "satisfied", "satisfied", "satisfied",
            "not_evaluable", "not_evaluable", "not_evaluable",
        ])
        self.assertEqual(result["topo_assessment_stage"].tolist(), [
            "interval", "interval", "interval", "interval", "endpoints", "distance",
            "structure",
        ])
        self.assertEqual(result.loc["12", "topo_assessment_site"], "b")
        self.assertEqual(result.loc["19", "eucl_assessment_upstream_reason"],
                         "identical_atom_endpoints")
        self.assertEqual(result.loc["23", "topo_assessment_reason"], "no_candidate_model")
        pd.testing.assert_frame_equal(result[before.columns], before)
        pd.testing.assert_frame_equal(data, before)

    def test_duplicate_and_non_contiguous_indices_keep_their_row_order(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "topo_dist": 45.0, "experiment": "condition_B"},
            {**self.row, "pLDDT_b": 60.0},
        ], index=["7", "7", "21_3"])
        before = data.copy(deep=True)
        result = self.annotate(data)
        pd.testing.assert_index_equal(result.index, before.index)
        pd.testing.assert_frame_equal(result[before.columns], before)
        self.assertEqual(result["topo_assessment_status"].tolist(),
                         ["satisfied", "violated", "not_evaluable"])

    def test_empty_input_has_the_same_assessment_headers(self):
        data = pd.DataFrame([self.row], index=["7"])
        expected_columns = self.annotate(data).columns
        result = self.annotate(data.iloc[:0])
        self.assertTrue(result.empty)
        self.assertEqual(result.columns.tolist(), expected_columns.tolist())
        self.assertIn("source_record_id", result.columns)
        self.assertEqual(result["eucl_assessment_minimum"].dtype, "float64")

    def test_policies_and_rule_versions_are_saved_for_each_method(self):
        result = self.annotate(pd.DataFrame([self.row], index=["7"]))
        row = result.iloc[0]
        self.assertEqual(row["eucl_assessment_minimum"], 5.0)
        self.assertEqual(row["eucl_assessment_maximum"], 30.0)
        self.assertEqual(row["eucl_assessment_boundary"], "inclusive")
        self.assertEqual(row["topo_assessment_maximum"], 35.0)
        self.assertEqual(row["topo_assessment_boundary"], "exclusive")
        for prefix, metric in (("eucl", "euclidean"), ("topo", "topological")):
            self.assertEqual(row[f"{prefix}_assessment_metric"], metric)
            self.assertEqual(row[f"{prefix}_assessment_plddt_cutoff"], 70.0)
            self.assertEqual(row[f"{prefix}_assessment_rule_version"], "candidate-distance-v1")
            self.assertEqual(row[f"{prefix}_assessment_interval_rule_version"], "distance-interval-v1")
            self.assertEqual(row[f"{prefix}_assessment_endpoint_rule_version"], "endpoint-eligibility-v1")
            self.assertEqual(row[f"{prefix}_assessment_chain_relation_rule_version"], "structure-chain-v1")
            self.assertEqual(row[f"{prefix}_assessment_upstream_distance_status"], "computed")

    def test_rerun_refreshes_old_assessments_and_structure_annotations(self):
        data = pd.DataFrame([self.row], index=["7"])
        first = self.annotate(data)
        self.assertEqual(first.iloc[0]["topo_assessment_status"], "satisfied")
        first["structure_source"] = "pdb"
        first["model_chain_relation"] = "same_chain"
        first["topo_assessment_reason"] = "stale"
        before = first.copy(deep=True)
        result = self.annotate(
            first,
            topological_policy={"minimum": 5.0, "maximum": 15.0, "boundary": "inclusive"},
        )
        self.assertEqual(result.iloc[0]["topo_assessment_status"], "violated")
        self.assertEqual(result.iloc[0]["topo_assessment_reason"], "exact_above_range")
        self.assertEqual(result.iloc[0]["topo_assessment_maximum"], 15.0)
        self.assertEqual(result.iloc[0]["structure_source"], "af_complex")
        self.assertEqual(result.iloc[0]["model_chain_relation"], "cross_chain")
        self.assertFalse(result.columns.duplicated().any())
        pd.testing.assert_frame_equal(first, before)

    def test_mixed_structure_sources_use_their_own_confidence_rules(self):
        pdb = {**self.row, "af_match_status": None, "af_model_id": None,
               "af_model_version": None, "pdb_id": "1ABC",
               "pdb_method": "X-RAY DIFFRACTION", "pLDDT_a": "-", "pLDDT_b": "-",
               "mapping_status_a": None, "mapping_status_b": None,
               "distance_status": None}
        monomer = {**pdb, "pdb_id": "afP1", "pdb_method": "ALPHAFOLD",
                   "pLDDT_a": "90.0", "pLDDT_b": "60.0"}
        result = self.annotate(pd.DataFrame([pdb, monomer, self.row],
                                           index=["7", "12", "20"]))
        self.assertEqual(result["structure_source"].tolist(),
                         ["pdb", "af_monomer", "af_complex"])
        self.assertEqual(result["topo_assessment_status"].tolist(),
                         ["satisfied", "not_evaluable", "satisfied"])
        self.assertEqual(result.loc["12", "topo_assessment_reason"], "low_plddt")
        self.assertIsNone(result.loc["7", "topo_assessment_upstream_distance_status"])

    def test_candidate_audit_keeps_assessments_for_excluded_rows(self):
        data = pd.DataFrame([
            self.row, {**self.row, "topo_dist": 45.0, "experiment": "condition_B"},
        ], index=["7", "7_1"])
        assessed = self.annotate(data)
        events = [{
            "removed_index": "7_1", "stage": "minimize_dataset",
            "rule": "lower_ranked_chain_candidate", "preferred_indices": ["7"],
            "rule_version": "legacy-m4-selection-v1",
        }]
        audit = build_candidate_audit(assessed, assessed.iloc[[0]], events)
        self.assertEqual(audit.loc["7_1", "m4_selection_status"], "excluded")
        self.assertEqual(audit.loc["7_1", "topo_assessment_status"], "violated")
        self.assertEqual(audit.loc["7_1", "experiment"], "condition_B")
        self.assertEqual(json.loads(audit.loc["7_1", "m4_exclusion_events"]), events)

    def test_invalid_configuration_is_rejected_even_when_there_are_no_candidates(self):
        empty = pd.DataFrame()
        for name in ("euclidean_policy", "topological_policy"):
            for policy in (
                None, {}, {"minimum": 5, "maximum": 35},
                {"minimum": 5, "maximum": 35, "boundary": "unknown"},
                {"minimum": 5, "maximum": 35, "boundary": "inclusive", "extra": 1},
                {"minimum": -1, "maximum": 35, "boundary": "inclusive"},
            ):
                with self.subTest(name=name, policy=policy):
                    with self.assertRaises(ValueError):
                        self.annotate(empty, **{name: policy})
        for cutoff in (float("nan"), 101.0, True, "70"):
            with self.subTest(cutoff=cutoff):
                with self.assertRaises(ValueError):
                    self.annotate(empty, plddt_cutoff=cutoff)

    def test_input_policy_dictionaries_are_not_modified(self):
        before_e = copy.deepcopy(self.euclidean)
        before_t = copy.deepcopy(self.topological)
        self.annotate(pd.DataFrame([self.row]))
        self.assertEqual(self.euclidean, before_e)
        self.assertEqual(self.topological, before_t)

    def test_duplicate_input_columns_are_rejected_as_ambiguous(self):
        data = pd.DataFrame([[15.0, 20.0]], columns=["topo_dist", "topo_dist"])
        with self.assertRaisesRegex(ValueError, "column"):
            self.annotate(data)

    def test_both_method_policies_are_required(self):
        with self.assertRaises(TypeError):
            annotate_distance_assessments(pd.DataFrame([self.row]), plddt_cutoff=70.0)


if __name__ == "__main__":
    unittest.main()
