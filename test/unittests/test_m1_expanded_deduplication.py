import unittest
from tempfile import TemporaryDirectory

import pandas as pd

from claudio.module01.src.algorithm.check_data import double_check_data
from claudio.utils.utils import build_xl_dataset


class TestM1ExpandedDeduplication(unittest.TestCase):
    def test_removed_source_does_not_reintroduce_candidates(self):
        record = {
            "unip_id_a": "P1", "unip_id_b": "P2",
            "seq_a": "MAKTAKT", "seq_b": "MGKT",
            "pep_a": "AKT", "pep_b": "GKT",
            "pos_a": 3, "pos_b": 3,
            "res_pos_a": 2, "res_pos_b": 2,
            "original_pos_a": 3, "original_pos_b": 3,
        }
        data = pd.DataFrame([
            {**record, "source_record_id": "record_000001"},
            {**record, "source_record_id": "record_000002"},
        ])

        with TemporaryDirectory() as output_dir:
            result = double_check_data(
                data,
                "expanded_duplicates",
                build_xl_dataset("K"),
                output_dir + "/",
                0,
            )

        self.assertEqual(len(result), 2)
        self.assertEqual(
            set(zip(result["pos_a"], result["pos_b"])),
            {(3, 3), (6, 3)},
        )
        self.assertEqual(
            set(result["source_record_id"]),
            {"record_000001"},
        )