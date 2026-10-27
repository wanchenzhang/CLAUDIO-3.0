"""Record the cleanup rule behind each removed structural candidate."""
import unittest

import pandas as pd

from claudio.utils.utils import clean_dataset


class TestM4CleanupSelectionLog(unittest.TestCase):
    def setUp(self):
        self.row = {
            "source_record_id": "record_000001",
            "candidate_id": "candidate_000001",
            "experiment": "condition_A",
            "source_crosslink_type": "Intra",
            "unip_id_a": "P1", "unip_id_b": "P1",
            "seq_a": "MAAKTGGK", "seq_b": "MAAKTGGK",
            "pep_a": "AAKT", "pep_b": "GGK",
            "pos_a": 4, "pos_b": 8, "res_pos_a": 3, "res_pos_b": 3,
            "pdb_id": "1ABC", "pdb_method": "X-RAY DIFFRACTION",
            "pdb_resolution": 2.0, "chain_a": "A", "chain_b": "B",
            "pdb_pos_a": 4, "pdb_pos_b": 8,
            "pLDDT_a": "-", "pLDDT_b": "-",
            "eucl_dist": 15.0, "topo_dist": 20.0,
            "homo_adjacency": 0.1, "homo_int_overl": 0.0,
            "homo_pep_overl": False,
            "evidence": "", "XL_type": "inter", "XL_confirmed": True,
            "swiss_model_homology": "",
        }

    def test_structural_duplicate_records_the_preferred_candidate(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "source_record_id": "record_000002",
             "experiment": "condition_B", "source_crosslink_type": "Homodimer"},
        ], index=["7", "12_1"])
        before = data.copy(deep=True)
        events = []
        result = clean_dataset(data, selection_log=events)

        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(result.loc["7", "source_record_id"], "record_000001")
        self.assertEqual(result.loc["7", "experiment"], "condition_A")
        self.assertEqual(len(events), 1)
        event = events[0]
        self.assertEqual(event["removed_index"], "12_1")
        self.assertEqual(event["stage"], "clean_dataset")
        self.assertEqual(event["rule"], "duplicate_structural_candidate")
        self.assertEqual(event["preferred_indices"], ["7"])
        self.assertTrue(event["rule_version"])
        pd.testing.assert_frame_equal(data, before)

    def test_unconfirmed_alternative_is_logged(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "chain_b": "C", "topo_dist": float("nan"),
             "XL_confirmed": False},
        ], index=["7", "7_1"])
        events = []
        result = clean_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["removed_index"], "7_1")
        self.assertEqual(events[0]["rule"], "unconfirmed_alternative")
        self.assertEqual(events[0]["preferred_indices"], ["7"])

    def test_evidence_bearing_alternative_is_logged_separately(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "chain_b": "C", "evidence": "'peptides overlap'"},
        ], index=["7", "7_1"])
        events = []
        result = clean_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["removed_index"], "7_1")
        self.assertEqual(events[0]["rule"], "alternative_with_evidence")
        self.assertEqual(events[0]["preferred_indices"], ["7"])

    def test_root_row_kept_by_legacy_cleanup_is_not_reported_as_removed(self):
        data = pd.DataFrame([
            {**self.row, "XL_confirmed": False, "topo_dist": float("nan")},
            {**self.row, "chain_b": "C"},
        ], index=["7", "7_1"])
        events = []
        result = clean_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7", "7_1"])
        self.assertEqual(events, [])

    def test_multiple_duplicate_candidates_each_get_a_record(self):
        data = pd.DataFrame([self.row, self.row, self.row],
                            index=["7", "7_1", "7_2"])
        events = []
        result = clean_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual([event["removed_index"] for event in events],
                         ["7_1", "7_2"])
        self.assertTrue(all(event["preferred_indices"] == ["7"] for event in events))

    def test_logging_preserves_legacy_selection_and_input(self):
        data = pd.DataFrame([
            self.row,
            self.row,
            {**self.row, "chain_b": "C", "topo_dist": float("nan"),
             "XL_confirmed": False},
        ], index=["7", "7_1", "7_2"])
        before = data.copy(deep=True)
        legacy_result = clean_dataset(data.copy())
        events = []
        result = clean_dataset(data, selection_log=events)
        pd.testing.assert_frame_equal(result, legacy_result)
        pd.testing.assert_frame_equal(data, before)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(
            [(event["removed_index"], event["rule"]) for event in events],
            [("7_1", "duplicate_structural_candidate"),
             ("7_2", "unconfirmed_alternative")],
        )

    def test_existing_log_entries_are_preserved_when_no_row_is_removed(self):
        events = [{"previous_stage": "kept"}]
        data = pd.DataFrame([self.row], index=["7"])
        result = clean_dataset(data, selection_log=events)
        self.assertEqual(result.index.tolist(), ["7"])
        self.assertEqual(events, [{"previous_stage": "kept"}])


if __name__ == "__main__":
    unittest.main()
