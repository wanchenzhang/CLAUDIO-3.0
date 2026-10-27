"""Keep a reason for candidates removed by M4's final selection."""
import unittest

import pandas as pd

from claudio.utils.utils import minimize_dataset


class TestM4MinimizeSelectionLog(unittest.TestCase):
    def setUp(self):
        self.row = {
            "source_record_id": "record_000001",
            "candidate_id": "candidate_000001",
            "experiment": "condition_A",
            "source_crosslink_type": "Inter",
            "unip_id_a": "P1", "unip_id_b": "P1",
            "seq_a": "MAAKTGGK", "seq_b": "MAAKTGGK",
            "pep_a": "AAKT", "pep_b": "GGK",
            "pos_a": 4, "pos_b": 8,
            "chain_a": "A", "chain_b": "B",
            "evidence": "", "XL_type": "inter", "topo_dist": 20.0,
        }

    def test_logs_candidates_with_a_different_type_from_the_root(self):
        for root_type, other_type in (("intra", "inter"), ("inter", "intra")):
            with self.subTest(root_type=root_type):
                data = pd.DataFrame([
                    {**self.row, "XL_type": root_type},
                    {**self.row, "XL_type": other_type, "topo_dist": 10.0},
                ], index=["7", "7_1"])
                events = []
                result = minimize_dataset(data, selection_log=events)
                self.assertEqual(result.index.tolist(), ["7"])
                self.assertEqual(len(events), 1)
                self.assertEqual(events[0]["removed_index"], "7_1")
                self.assertEqual(events[0]["stage"], "minimize_dataset")
                self.assertEqual(events[0]["rule"], "different_crosslink_type")
                self.assertEqual(events[0]["preferred_indices"], ["7"])
                self.assertTrue(events[0]["rule_version"])

    def test_records_evidence_priority_even_when_its_distance_is_longer(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "topo_dist": 5.0, "evidence": "'peptides overlap'"},
            {**self.row, "topo_dist": 25.0},
        ], index=["7", "7_1", "7_2"])
        events = []
        result = minimize_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(
            [(event["removed_index"], event["rule"], event["preferred_indices"])
             for event in events],
            [("7_1", "lower_ranked_chain_candidate", ["7"]),
             ("7_2", "lower_ranked_chain_candidate", ["7"])],
        )

    def test_prefers_an_available_distance_when_all_candidates_have_evidence(self):
        data = pd.DataFrame([
            {**self.row, "evidence": "'peptides overlap'", "topo_dist": float("nan")},
            {**self.row, "evidence": "'distance'", "topo_dist": 10.0},
        ], index=["7", "7_1"])
        events = []
        result = minimize_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7_1"])
        self.assertEqual(events[0]["removed_index"], "7")
        self.assertEqual(events[0]["preferred_indices"], ["7_1"])
        self.assertEqual(events[0]["rule"], "lower_ranked_chain_candidate")

    def test_keeps_the_root_fallback_out_of_the_exclusion_log(self):
        for evidence in ("", "'distance'"):
            with self.subTest(evidence=evidence):
                data = pd.DataFrame([
                    {**self.row, "evidence": evidence, "topo_dist": float("nan")},
                    {**self.row, "evidence": evidence, "topo_dist": float("nan")},
                ], index=["7", "7_1"])
                events = []
                result = minimize_dataset(data, selection_log=events)
                self.assertEqual(result.index.tolist(), ["7"])
                self.assertEqual([event["removed_index"] for event in events], ["7_1"])
                self.assertEqual(events[0]["preferred_indices"], ["7"])

    def test_records_the_shorter_intra_candidate(self):
        data = pd.DataFrame([
            {**self.row, "XL_type": "intra"},
            {**self.row, "XL_type": "intra", "topo_dist": 12.0},
        ], index=["7", "7_1"])
        events = []
        result = minimize_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7_1"])
        self.assertEqual(events[0]["removed_index"], "7")
        self.assertEqual(events[0]["rule"], "lower_ranked_chain_candidate")
        self.assertEqual(events[0]["preferred_indices"], ["7_1"])

    def test_equal_distances_keep_the_first_candidate(self):
        data = pd.DataFrame([self.row, self.row], index=["7", "7_1"])
        events = []
        result = minimize_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(events[0]["removed_index"], "7_1")
        self.assertEqual(events[0]["preferred_indices"], ["7"])

    def test_logs_repeated_peptide_sources_without_changing_their_conditions(self):
        data = pd.DataFrame([
            {**self.row, "seq_a": "MAAKTAAKTGGK"},
            {**self.row, "seq_a": "MAAKTAAKTGGK",
             "source_record_id": "record_000002", "experiment": "condition_B",
             "source_crosslink_type": "Homodimer"},
        ], index=["7", "12"])
        before = data.copy(deep=True)
        events = []
        result = minimize_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["removed_index"], "12")
        self.assertEqual(events[0]["rule"], "repeated_peptide_source")
        self.assertEqual(events[0]["preferred_indices"], ["7"])
        self.assertEqual(events[0]["stage"], "minimize_dataset")
        self.assertEqual(data.loc["12", "experiment"], "condition_B")
        pd.testing.assert_frame_equal(data, before)

    def test_records_chain_selection_before_repeated_source_removal(self):
        data = pd.DataFrame([
            {**self.row, "seq_a": "MAAKTAAKTGGK", "topo_dist": 15.0},
            {**self.row, "seq_a": "MAAKTAAKTGGK",
             "source_record_id": "record_000002", "experiment": "condition_B",
             "topo_dist": 30.0},
            {**self.row, "seq_a": "MAAKTAAKTGGK",
             "source_record_id": "record_000002", "experiment": "condition_B",
             "topo_dist": 10.0},
        ], index=["7", "12", "12_1"])
        before = data.copy(deep=True)
        legacy_result = minimize_dataset(data.copy(deep=True))
        events = []
        result = minimize_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(
            [(event["removed_index"], event["rule"], event["preferred_indices"])
             for event in events],
            [("12", "lower_ranked_chain_candidate", ["12_1"]),
             ("12_1", "repeated_peptide_source", ["7"])],
        )
        pd.testing.assert_frame_equal(result, legacy_result)
        pd.testing.assert_frame_equal(data, before)

    def test_keeps_existing_log_entries_when_nothing_is_removed(self):
        for data in (
            pd.DataFrame([self.row], index=["7"]),
            pd.DataFrame(columns=self.row),
        ):
            with self.subTest(rows=len(data)):
                events = [{"previous_stage": "clean_dataset"}]
                result = minimize_dataset(data, selection_log=events)
                pd.testing.assert_frame_equal(result, data)
                self.assertEqual(events, [{"previous_stage": "clean_dataset"}])


if __name__ == "__main__":
    unittest.main()
