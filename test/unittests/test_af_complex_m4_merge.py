"""Tests for linking AF chain candidates to their M3 source rows."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from claudio.module02.src_structure_search.algorithm.af_complex_search import (
    assign_complex_chain_pairs,
)
from claudio.module04.src.io.read_ins import read_inputs


class TestAFComplexM4Merge(unittest.TestCase):
    def test_merges_m3_scores_into_the_correct_chain_candidates(self):
        chains = [
            {"chain_id": "A", "protein_id": "P1", "sequence": "MAAKT"},
            {"chain_id": "C", "protein_id": "P1", "sequence": "MAAKT"},
            {"chain_id": "B", "protein_id": "P2", "sequence": "MGKT"},
            {"chain_id": "D", "protein_id": "P3", "sequence": "MVKT"},
        ]
        data = pd.DataFrame([
            {
                "unip_id_a": "P1",
                "unip_id_b": "P2",
                "seq_a": "MAAKT",
                "seq_b": "MGKT",
                "source_record_id": "record_000008",
            },
            {
                "unip_id_a": "P3",
                "unip_id_b": "P2",
                "seq_a": "MVKT",
                "seq_b": "MGKT",
                "source_record_id": "record_000013",
            },
        ], index=[7, 12])

        structural = assign_complex_chain_pairs(data, chains)

        # Distinct marker values test association, not biological scoring.
        ops = data.copy()
        ops["homo_adjacency"] = [0.123, 0.789]
        ops["homo_int_overl"] = [0.234, 0.890]
        ops["homo_pep_overl"] = [True, True]

        with TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir)
            structural_path = folder / "structural.csv"
            ops_path = folder / "ops.csv"
            structural.to_csv(structural_path)
            ops.to_csv(ops_path)

            result = read_inputs(str(structural_path), str(ops_path))

        self.assertEqual(len(result), 3)
        self.assertTrue(result.index.is_unique)

        expected = {
            "7": ("record_000008", 0.123, 0.234),
            "7_1": ("record_000008", 0.123, 0.234),
            "12": ("record_000013", 0.789, 0.890),
        }
        self.assertEqual(set(map(str, result.index)), set(expected))

        for index, row in result.iterrows():
            source_id, adjacency, overlap = expected[str(index)]
            self.assertEqual(row["source_record_id"], source_id)
            self.assertAlmostEqual(row["homo_adjacency"], adjacency)
            self.assertAlmostEqual(row["homo_int_overl"], overlap)
    def test_preserves_missing_overlap_flag_for_unresolved_positions(self):
        structural = pd.DataFrame([{
            "chain_a": "A",
            "chain_b": "B",
            "position_status_a": "ambiguous",
            "position_status_b": "valid",
            "source_record_id": "record_000001",
        }])

        ops = pd.DataFrame([{
            "source_record_id": "record_000001",
            "homo_adjacency": float("nan"),
            "homo_int_overl": float("nan"),
            "homo_pep_overl": pd.NA,
        }])

        with TemporaryDirectory() as temp_dir:
            folder = Path(temp_dir)
            structural_path = folder / "structural.csv"
            ops_path = folder / "ops.csv"
            structural.to_csv(structural_path)
            ops.to_csv(ops_path)

            result = read_inputs(str(structural_path), str(ops_path))

        self.assertEqual(len(result), 1)
        row = result.iloc[0]
        self.assertTrue(pd.isna(row["homo_adjacency"]))
        self.assertTrue(pd.isna(row["homo_int_overl"]))
        self.assertTrue(
            pd.isna(row["homo_pep_overl"]),
            "M4 must not convert an unevaluated overlap flag to False",
        )

if __name__ == "__main__":
    unittest.main()
