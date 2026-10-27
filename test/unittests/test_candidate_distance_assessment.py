"""Join model context, endpoint gates and method-specific distance assessment."""
import unittest

import pandas as pd

from claudio.module04.src.algorithm.distance_assessment import (
    assess_candidate_distance,
)


class TestCandidateDistanceAssessment(unittest.TestCase):
    def setUp(self):
        self.row = pd.Series({
            "source_record_id": "record_000001", "candidate_id": "candidate_000001",
            "experiment": "condition_A", "source_crosslink_type": "Homodimer",
            "XL_type": "inter", "score_XL_type": "intra",
            "af_match_status": "matched", "af_model_id": "AF-model",
            "af_model_version": 1, "path": "/archive/model.cif",
            "pdb_id": "AF-model", "pdb_method": "ALPHAFOLD",
            "chain_a": "AA", "chain_b": "AB",
            "position_status_a": "valid", "position_status_b": "corrected",
            "mapping_status_a": "mapped", "mapping_status_b": "mapped",
            "mapping_reason_a": None, "mapping_reason_b": None,
            "res_crit_a": True, "res_crit_b": True,
            "pdb_pos_a": 104, "pdb_pos_b": 108,
            "pLDDT_a": 90.0, "pLDDT_b": 90.0,
            "distance_status": "computed", "distance_reason": None,
            "eucl_dist": 15.0, "topo_dist": 20.0,
            "eucl_dist_lower_bound": None, "topo_dist_lower_bound": None,
        })

    def assess(self, row=None, metric="topological", **changes):
        options = {
            "metric": metric, "minimum": 5.0, "maximum": 35.0,
            "boundary": "inclusive", "plddt_cutoff": 70.0,
        }
        options.update(changes)
        return assess_candidate_distance(self.row if row is None else row, **options)

    def test_each_distance_method_uses_its_own_value_and_interval(self):
        row = self.row.copy()
        row["topo_dist"] = 45.0
        euclidean = self.assess(row, "euclidean", maximum=30.0)
        topological = self.assess(row)
        self.assertEqual(euclidean["status"], "satisfied")
        self.assertEqual(euclidean["maximum"], 30.0)
        self.assertEqual(topological["status"], "violated")
        self.assertEqual(topological["reason"], "exact_above_range")
        self.assertEqual(topological["maximum"], 35.0)

    def test_a_bound_for_one_method_does_not_replace_the_other_exact_distance(self):
        row = self.row.copy()
        row["distance_status"] = "lower_bound"
        row["topo_dist"] = None
        row["topo_dist_lower_bound"] = 75.328
        self.assertEqual(self.assess(row, "euclidean")["status"], "satisfied")
        result = self.assess(row)
        self.assertEqual(result["status"], "violated")
        self.assertEqual(result["reason"], "lower_bound_above_range")
        row["topo_dist_lower_bound"] = 20.0
        self.assertEqual(self.assess(row)["status"], "indeterminate")

    def test_upstream_not_evaluable_blocks_stale_numbers_for_both_methods(self):
        row = self.row.copy()
        row["distance_status"] = "not_evaluable"
        row["distance_reason"] = "identical_atom_endpoints"
        for metric in ("euclidean", "topological"):
            with self.subTest(metric=metric):
                result = self.assess(row, metric)
                self.assertEqual(result["status"], "not_evaluable")
                self.assertEqual(result["reason"], "upstream_distance_not_evaluable")
                self.assertEqual(result["upstream_reason"], "identical_atom_endpoints")
                self.assertEqual(result["assessment_stage"], "distance")

    def test_af_complex_requires_a_distance_stage_status(self):
        for value in (None, pd.NA, "-"):
            with self.subTest(value=value):
                row = self.row.copy()
                row["distance_status"] = value
                result = self.assess(row)
                self.assertEqual(result["status"], "not_evaluable")
                self.assertEqual(result["reason"], "missing_distance_status")

    def test_unknown_distance_status_is_not_silently_ignored(self):
        row = self.row.copy()
        row["distance_status"] = "failed"
        result = self.assess(row)
        self.assertEqual(result["status"], "not_evaluable")
        self.assertEqual(result["reason"], "unrecognized_distance_status")

    def test_inconsistent_computed_and_lower_bound_states_are_blocked(self):
        for field, value in (
            ("topo_dist_lower_bound", 10.0),
            ("topo_dist", None),
        ):
            with self.subTest(field=field):
                row = self.row.copy()
                row[field] = value
                result = self.assess(row)
                self.assertEqual(result["status"], "not_evaluable")
                self.assertEqual(result["reason"], "inconsistent_distance_status")
        row = self.row.copy()
        row["distance_status"] = "lower_bound"
        result = self.assess(row)
        self.assertEqual(result["reason"], "inconsistent_distance_status")

    def test_missing_model_path_chain_and_af_identity_block_stale_distances(self):
        cases = (
            ("path", "-", "no_usable_structure"),
            ("chain_b", None, "missing_chain_assignment"),
            ("af_model_id", None, "missing_model_identity"),
            ("af_model_version", None, "missing_model_identity"),
            ("af_match_status", "no_candidate_model", "no_candidate_model"),
        )
        for field, value, reason in cases:
            with self.subTest(field=field):
                row = self.row.copy()
                row[field] = value
                result = self.assess(row)
                self.assertEqual(result["status"], "not_evaluable")
                self.assertEqual(result["reason"], reason)
                self.assertEqual(result["assessment_stage"], "structure")

    def test_refreshes_stale_structure_annotations(self):
        row = self.row.copy()
        row["structure_source"] = "pdb"
        row["model_chain_relation"] = "same_chain"
        result = self.assess(row)
        self.assertEqual(result["structure_source"], "af_complex")
        self.assertEqual(result["model_chain_relation"], "cross_chain")
        row["chain_b"] = "-"
        row["model_chain_relation"] = "cross_chain"
        self.assertEqual(self.assess(row)["status"], "not_evaluable")

    def test_unresolved_endpoint_and_low_confidence_stop_before_interval_assessment(self):
        cases = (
            ("position_status_b", "ambiguous", "unresolved_position"),
            ("mapping_status_b", "not_evaluable", "unmapped_endpoint"),
            ("pLDDT_b", 60.0, "low_plddt"),
        )
        for field, value, reason in cases:
            with self.subTest(field=field):
                row = self.row.copy()
                row[field] = value
                result = self.assess(row)
                self.assertEqual(result["status"], "not_evaluable")
                self.assertEqual(result["reason"], reason)
                self.assertEqual(result["site"], "b")
                self.assertEqual(result["assessment_stage"], "endpoints")

    def test_legacy_pdb_and_af_monomer_can_use_distances_without_a_stage_status(self):
        row = self.row.drop([
            "af_match_status", "af_model_id", "af_model_version",
            "mapping_status_a", "mapping_status_b", "distance_status",
        ])
        for method, source in (
            ("X-RAY DIFFRACTION", "pdb"), ("ALPHAFOLD", "af_monomer"),
        ):
            with self.subTest(source=source):
                row["pdb_method"] = method
                row["pdb_id"] = "1ABC" if source == "pdb" else "afP1"
                result = self.assess(row)
                self.assertEqual(result["status"], "satisfied")
                self.assertEqual(result["structure_source"], source)
                self.assertIsNone(result["upstream_distance_status"])
        row["distance_status"] = "not_evaluable"
        self.assertEqual(self.assess(row)["status"], "not_evaluable")

    def test_unknown_structure_context_is_not_evaluable(self):
        row = self.row.drop([
            "af_match_status", "af_model_id", "af_model_version",
            "pdb_id", "pdb_method",
        ])
        result = self.assess(row)
        self.assertEqual(result["status"], "not_evaluable")
        self.assertEqual(result["reason"], "unknown_structure_source")

    def test_preserves_sources_labels_conditions_and_raw_input(self):
        row = self.row.copy()
        row["topo_dist_raw"] = "20.000"
        before = row.copy(deep=True)
        result = self.assess(row)
        pd.testing.assert_series_equal(row, before)
        self.assertEqual(result["status"], "satisfied")
        self.assertEqual(row["source_crosslink_type"], "Homodimer")
        self.assertEqual(row["XL_type"], "inter")
        self.assertEqual(row["score_XL_type"], "intra")

    def test_records_all_applied_policy_versions_and_the_method(self):
        result = self.assess(boundary="exclusive")
        self.assertEqual(result["metric"], "topological")
        self.assertEqual(result["boundary"], "exclusive")
        self.assertEqual(result["plddt_cutoff"], 70.0)
        self.assertEqual(result["rule_version"], "candidate-distance-v1")
        self.assertEqual(result["interval_rule_version"], "distance-interval-v1")
        self.assertEqual(result["endpoint_rule_version"], "endpoint-eligibility-v1")
        self.assertEqual(result["chain_relation_rule_version"], "structure-chain-v1")
        self.assertEqual(result["upstream_distance_status"], "computed")

    def test_invalid_configuration_is_rejected_even_for_an_unusable_candidate(self):
        row = self.row.copy()
        row["path"] = "-"
        for options in (
            {"metric": "combined"}, {"minimum": -1.0},
            {"boundary": None}, {"plddt_cutoff": 101},
        ):
            with self.subTest(options=options):
                with self.assertRaises(ValueError):
                    self.assess(row, **options)


if __name__ == "__main__":
    unittest.main()
