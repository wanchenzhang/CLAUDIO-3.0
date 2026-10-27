"""Build a complete candidate table from M4 selection events."""
import json
import unittest

import pandas as pd

from claudio.module04.src.algorithm.candidate_audit import build_candidate_audit


class TestM4CandidateAudit(unittest.TestCase):
    def setUp(self):
        self.candidates = pd.DataFrame({
            "source_record_id": ["record_000001", "record_000002", "record_000002"],
            "candidate_id": ["candidate_000001"] * 3,
            "experiment": ["condition_A", "condition_B", "condition_B"],
            "source_crosslink_type": ["Inter", "Homodimer", "Homodimer"],
            "topo_dist_raw": ["18.180", ">75.328", "20.000"],
            "topo_dist": [18.18, float("nan"), 20.0],
            "mapping_status_a": ["mapped", "not_evaluable", "mapped"],
        }, index=["7", "12", "12_1"])
        self.events = [
            {"removed_index": "12", "stage": "clean_dataset",
             "rule": "duplicate_structural_candidate",
             "preferred_indices": ["12_1"], "rule_version": "legacy-m4-selection-v1"},
            {"removed_index": "12_1", "stage": "minimize_dataset",
             "rule": "repeated_peptide_source",
             "preferred_indices": ["7"], "rule_version": "legacy-m4-selection-v1"},
        ]

    def test_keeps_excluded_sources_conditions_and_raw_distance_text(self):
        before = self.candidates.copy(deep=True)
        result = build_candidate_audit(
            self.candidates, self.candidates.loc[["7"]], self.events,
        )
        pd.testing.assert_frame_equal(result[before.columns], before)
        pd.testing.assert_frame_equal(self.candidates, before)
        self.assertEqual(result["m4_selection_status"].tolist(),
                         ["retained", "excluded", "excluded"])
        self.assertEqual(json.loads(result.loc["7", "m4_exclusion_events"]), [])
        self.assertEqual(json.loads(result.loc["12", "m4_exclusion_events"]),
                         [self.events[0]])
        self.assertEqual(json.loads(result.loc["12_1", "m4_exclusion_events"]),
                         [self.events[1]])
        self.assertTrue(result["m4_audit_rule_version"].eq("m4-candidate-audit-v1").all())

    def test_keeps_multiple_exclusion_reasons_in_event_order(self):
        extra = {**self.events[0], "rule": "unconfirmed_alternative"}
        result = build_candidate_audit(
            self.candidates, self.candidates.loc[["7"]],
            [self.events[0], extra, self.events[1]],
        )
        self.assertEqual(json.loads(result.loc["12", "m4_exclusion_events"]),
                         [self.events[0], extra])

    def test_retained_does_not_mean_the_mapping_was_evaluable(self):
        data = self.candidates.loc[["12"]]
        result = build_candidate_audit(data, data, [])
        self.assertEqual(result.iloc[0]["m4_selection_status"], "retained")
        self.assertEqual(result.iloc[0]["mapping_status_a"], "not_evaluable")
        self.assertTrue(pd.isna(result.iloc[0]["topo_dist"]))

    def test_rejects_removed_rows_without_an_exclusion_event(self):
        with self.assertRaisesRegex(ValueError, "exclusion"):
            build_candidate_audit(self.candidates, self.candidates.loc[["7"]], [])

    def test_rejects_events_for_retained_or_unknown_candidates(self):
        for removed in ("7", "99"):
            with self.subTest(removed=removed):
                events = [*self.events, {**self.events[0], "removed_index": removed}]
                with self.assertRaisesRegex(ValueError, "exclusion"):
                    build_candidate_audit(
                        self.candidates, self.candidates.loc[["7"]], events,
                    )

    def test_rejects_unknown_preferred_candidate_and_unknown_selected_row(self):
        with self.assertRaisesRegex(ValueError, "preferred"):
            build_candidate_audit(
                self.candidates, self.candidates.loc[["7"]],
                [{**self.events[0], "preferred_indices": ["99"]}, self.events[1]],
            )
        unknown = self.candidates.loc[["7"]].rename(index={"7": "99"})
        with self.assertRaisesRegex(ValueError, "selected"):
            build_candidate_audit(self.candidates, unknown, [])

    def test_rejects_ambiguous_index_labels(self):
        for labels in (["7", "7", "12"], [7, "7", "12"]):
            with self.subTest(labels=labels):
                data = self.candidates.copy()
                data.index = labels
                with self.assertRaisesRegex(ValueError, "unique"):
                    build_candidate_audit(data, data, [])

    def test_empty_audit_still_has_headers(self):
        empty = self.candidates.iloc[:0]
        result = build_candidate_audit(empty, empty, [])
        self.assertTrue(result.empty)
        for column in (
            "source_record_id", "experiment", "m4_selection_status",
            "m4_exclusion_events", "m4_audit_rule_version",
        ):
            self.assertIn(column, result.columns)

    def test_does_not_overwrite_input_columns_with_audit_fields(self):
        for column in (
            "m4_selection_status", "m4_exclusion_events", "m4_audit_rule_version",
        ):
            with self.subTest(column=column):
                data = self.candidates.assign(**{column: "input value"})
                with self.assertRaisesRegex(ValueError, column):
                    build_candidate_audit(data, data, [])


if __name__ == "__main__":
    unittest.main()
