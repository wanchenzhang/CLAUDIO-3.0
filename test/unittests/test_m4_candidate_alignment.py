"""M3 evidence follows its parent source row for ordinary structures too."""
import unittest

import pandas as pd

from claudio.module04.src.io.read_ins import merge_datasets


class TestM4CandidateAlignment(unittest.TestCase):
    def setUp(self):
        self.ops = pd.DataFrame({
            "source_record_id": ["record_000001", "record_000002"],
            "unip_id_a": ["P1", "P2"], "unip_id_b": ["P1", "P2"],
            "pos_a": [3, 7], "pos_b": [9, 12],
            "homo_adjacency": [0.1, 0.8],
            "homo_int_overl": [0.0, 0.5], "homo_pep_overl": [False, True],
        }, index=[7, 12])
        self.structural = self.ops.loc[[12, 7]].drop(
            columns=["homo_adjacency", "homo_int_overl", "homo_pep_overl"],
        )
        self.structural["pdb_id"] = "1ABC"
        self.structural["chain_a"] = "A"
        self.structural["chain_b"] = "B"

    def test_reordered_ordinary_structure_rows_receive_their_own_m3_evidence(self):
        before = self.structural.copy(deep=True)
        result = merge_datasets(self.structural, self.ops)
        self.assertEqual(result.index.tolist(), [12, 7])
        self.assertEqual(result["homo_adjacency"].tolist(), [0.8, 0.1])
        self.assertEqual(result["homo_pep_overl"].tolist(), [True, False])
        pd.testing.assert_frame_equal(self.structural, before)

    def test_mismatched_source_or_position_is_rejected(self):
        for column, value in (("source_record_id", "wrong_source"), ("pos_a", 99)):
            with self.subTest(column=column):
                wrong = self.structural.copy()
                wrong.loc[12, column] = value
                with self.assertRaisesRegex(ValueError, column):
                    merge_datasets(wrong, self.ops)

    def test_empty_ordinary_table_retains_m3_headers(self):
        result = merge_datasets(self.structural.iloc[:0], self.ops.iloc[:0])
        self.assertTrue(result.empty)
        self.assertIn("homo_pep_overl", result)
        self.assertIn("source_record_id", result)


if __name__ == "__main__":
    unittest.main()
