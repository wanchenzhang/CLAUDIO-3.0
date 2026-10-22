"""Associate AF structural candidates with their original M3 rows."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from claudio.module04.src.io.read_ins import read_inputs


class TestAFComplexM4Alignment(unittest.TestCase):
    def setUp(self):
        self.ops = pd.DataFrame([
            {
                "source_record_id": "record_000001",
                "candidate_id": "candidate_000001",
                "unip_id_a": "P1", "unip_id_b": "P1",
                "pos_a": 3, "pos_b": 9, "experiment": "condition_A",
                "homo_adjacency": 0.1, "homo_int_overl": 0.2,
                "homo_pep_overl": True,
            },
            {
                "source_record_id": "record_000002",
                "candidate_id": "candidate_000002",
                "unip_id_a": "P1", "unip_id_b": "P1",
                "pos_a": 6, "pos_b": 9, "experiment": "condition_B",
                "homo_adjacency": 0.8, "homo_int_overl": 0.9,
                "homo_pep_overl": False,
            },
        ], index=[7, 12])
        self.structural = self.ops.drop(columns=[
            "homo_adjacency", "homo_int_overl", "homo_pep_overl",
        ]).assign(af_match_status="matched", chain_a="A", chain_b="A")

    def read_pair(self, structural, ops):
        with TemporaryDirectory() as folder:
            root = Path(folder)
            structural.to_csv(root / "structures.csv")
            ops.to_csv(root / "ops.csv")
            return read_inputs(str(root / "structures.csv"), str(root / "ops.csv"))

    def test_reordered_integer_indices_keep_their_own_m3_results(self):
        structural = self.structural.loc[[12, 7]]
        result = self.read_pair(structural, self.ops)

        pd.testing.assert_frame_equal(result[structural.columns], structural)
        self.assertEqual(result["homo_adjacency"].tolist(), [0.8, 0.1])
        self.assertEqual(result["homo_int_overl"].tolist(), [0.9, 0.2])
        self.assertEqual(result["homo_pep_overl"].tolist(), [False, True])

    def test_model_and_chain_expansion_share_only_their_original_m3_row(self):
        structural = pd.concat([
            self.structural.loc[[12]],
            self.structural.loc[[7]],
            self.structural.loc[[7]].rename(index={7: "7_1"}),
            self.structural.loc[[12]].rename(index={12: "12_1_1"}),
        ])
        structural["chain_b"] = ["A", "A", "B", "B"]
        result = self.read_pair(structural, self.ops)

        self.assertEqual(result.index.astype(str).tolist(), ["12", "7", "7_1", "12_1_1"])
        self.assertEqual(result["homo_adjacency"].tolist(), [0.8, 0.1, 0.1, 0.8])
        self.assertEqual(result["source_record_id"].tolist(),
                         structural["source_record_id"].tolist())
        self.assertEqual(result["experiment"].tolist(),
                         structural["experiment"].tolist())

    def test_missing_m3_source_is_rejected_instead_of_borrowing_another_row(self):
        with self.assertRaisesRegex(ValueError, "M3"):
            self.read_pair(self.structural, self.ops.loc[[7]])

    def test_mismatched_source_or_position_is_rejected(self):
        for column, value in (
            ("source_record_id", "another_source"),
            ("candidate_id", "another_candidate"),
            ("pos_a", 99),
        ):
            with self.subTest(column=column):
                changed = self.ops.copy()
                changed.loc[7, column] = value
                with self.assertRaisesRegex(ValueError, column):
                    self.read_pair(self.structural, changed)

    def test_duplicate_m3_indices_are_rejected(self):
        duplicated = pd.concat([self.ops, self.ops.loc[[7]]])
        with self.assertRaisesRegex(ValueError, "unique"):
            self.read_pair(self.structural, duplicated)

    def test_existing_ops_values_are_refreshed_from_the_matching_m3_row(self):
        structural = self.structural.assign(
            homo_adjacency=99.0, homo_int_overl=99.0, homo_pep_overl=False,
        )
        result = self.read_pair(structural, self.ops)
        self.assertEqual(result["homo_adjacency"].tolist(), [0.1, 0.8])
        self.assertEqual(result["homo_int_overl"].tolist(), [0.2, 0.9])
        self.assertEqual(result["homo_pep_overl"].tolist(), [True, False])

    def test_empty_af_tables_keep_m3_result_columns(self):
        result = self.read_pair(self.structural.iloc[:0], self.ops.iloc[:0])
        self.assertTrue(result.empty)
        for column in ("homo_adjacency", "homo_int_overl", "homo_pep_overl"):
            self.assertIn(column, result.columns)


if __name__ == "__main__":
    unittest.main()
