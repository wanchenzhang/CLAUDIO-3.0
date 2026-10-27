"""Keep native distance tokens intact when M4 reads structural results."""
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from claudio.module04.src.io.read_ins import read_inputs


class TestM4RawDistanceTokens(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.structural = pd.DataFrame({
            "source_record_id": ["record_000001", "record_000002", "record_000003"],
            "af_match_status": ["matched"] * 3,
            "chain_a": ["A"] * 3,
            "chain_b": ["A"] * 3,
            "eucl_dist": [18.18, 42.32, float("nan")],
            "topo_dist": [18.22, 42.74, float("nan")],
            "eucl_dist_raw": ["18.180", "42.320", None],
            "topo_dist_raw": ["18.220", "42.740", None],
            "topo_dist_lower_bound": [float("nan")] * 3,
        }, index=[7, 12, 20])
        self.ops = self.structural[["source_record_id"]].assign(
            homo_adjacency=[0.1, 0.2, float("nan")],
            homo_int_overl=[0.0, 0.0, float("nan")],
            homo_pep_overl=[False, True, None],
        )

    def read_and_save(self, structural):
        structural.to_csv(self.root / "structural.csv")
        self.ops.to_csv(self.root / "ops.csv")
        result = read_inputs(
            str(self.root / "structural.csv"), str(self.root / "ops.csv"),
        )
        result.to_csv(self.root / "merged.csv")
        saved = pd.read_csv(
            self.root / "merged.csv", index_col=0,
            dtype={"eucl_dist_raw": str, "topo_dist_raw": str},
        )
        return result, saved

    def test_numeric_tokens_keep_trailing_zeros_through_read_and_write(self):
        result, saved = self.read_and_save(self.structural)
        for frame in (result, saved):
            self.assertEqual(frame["eucl_dist_raw"].iloc[:2].tolist(),
                             ["18.180", "42.320"])
            self.assertEqual(frame["topo_dist_raw"].iloc[:2].tolist(),
                             ["18.220", "42.740"])
            self.assertTrue(pd.api.types.is_numeric_dtype(frame["eucl_dist"]))
            self.assertTrue(pd.api.types.is_numeric_dtype(frame["topo_dist"]))
            self.assertEqual(frame["eucl_dist"].iloc[0], 18.18)
            self.assertEqual(frame["topo_dist"].iloc[1], 42.74)
            self.assertEqual(frame["source_record_id"].tolist(),
                             ["record_000001", "record_000002", "record_000003"])

    def test_bound_marker_and_missing_exact_distance_are_preserved(self):
        structural = self.structural.copy()
        structural["topo_dist_raw"] = [">75.328", "18.220", None]
        structural["topo_dist"] = [float("nan"), 18.22, float("nan")]
        structural["topo_dist_lower_bound"] = [75.328, float("nan"), float("nan")]
        result, saved = self.read_and_save(structural)
        for frame in (result, saved):
            self.assertEqual(frame["topo_dist_raw"].iloc[:2].tolist(),
                             [">75.328", "18.220"])
            self.assertTrue(pd.isna(frame.loc[7, "topo_dist"]))
            self.assertEqual(frame.loc[7, "topo_dist_lower_bound"], 75.328)
            self.assertTrue(pd.isna(frame.loc[20, "topo_dist_raw"]))

    def test_missing_tokens_are_written_as_empty_fields(self):
        structural = self.structural.copy()
        structural["eucl_dist_raw"] = None
        structural["topo_dist_raw"] = None
        result, saved = self.read_and_save(structural)
        for frame in (result, saved):
            self.assertTrue(frame[["eucl_dist_raw", "topo_dist_raw"]].isna().all().all())
        literal = pd.read_csv(self.root / "merged.csv", index_col=0, keep_default_na=False)
        self.assertEqual(literal["eucl_dist_raw"].tolist(), ["", "", ""])
        self.assertEqual(literal["topo_dist_raw"].tolist(), ["", "", ""])

    def test_inputs_without_raw_columns_keep_their_existing_schema(self):
        for af_mode in (False, True):
            with self.subTest(af_mode=af_mode):
                structural = self.structural.drop(
                    columns=["eucl_dist_raw", "topo_dist_raw"],
                )
                if not af_mode:
                    structural = structural.drop(columns=["af_match_status"])
                result, _ = self.read_and_save(structural)
                self.assertNotIn("eucl_dist_raw", result)
                self.assertNotIn("topo_dist_raw", result)
                self.assertEqual(result.index.tolist(), [7, 12, 20])
                self.assertEqual(result["homo_adjacency"].iloc[:2].tolist(), [0.1, 0.2])
                self.assertEqual(result["eucl_dist"].iloc[:2].tolist(), [18.18, 42.32])


if __name__ == "__main__":
    unittest.main()
