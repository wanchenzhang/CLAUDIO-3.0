"""All M4 outputs retain model candidates and clearly scoped relationship labels."""
import importlib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd
from click.testing import CliRunner


m4 = importlib.import_module("claudio.module04.src.main")


class TestModelRelationshipCLI(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        base = {
            "source_record_id": "record_000001", "candidate_id": "candidate_000001",
            "experiment": "condition_A", "source_crosslink_type": "Homodimer",
            "unip_id_a": "P1", "unip_id_b": "P1",
            "seq_a": "KKK", "seq_b": "KKK", "pep_a": "K", "pep_b": "K",
            "pos_a": 1, "pos_b": 2, "res_pos_a": 1, "res_pos_b": 1,
            "position_status_a": "valid", "position_status_b": "valid",
            "pdb_id": "1ABC", "pdb_method": "X-RAY DIFFRACTION",
            "pdb_resolution": 2.0, "path": "/archive/1ABC.cif",
            "chain_a": "A", "chain_b": "A", "pdb_pos_a": 1, "pdb_pos_b": 2,
            "res_crit_a": True, "res_crit_b": True, "pLDDT_a": "-", "pLDDT_b": "-",
            "eucl_dist": 55.0, "topo_dist": 60.0,
            "eucl_dist_raw": "55.000", "topo_dist_raw": "60.000",
        }
        self.data = pd.DataFrame([
            base,
            {**base, "chain_b": "B", "eucl_dist": 15.0, "topo_dist": 20.0},
            {**base, "chain_b": "B", "eucl_dist": 15.0, "topo_dist": 20.0,
             "source_record_id": "record_000002", "experiment": "condition_B"},
            {**base, "position_status_a": "ambiguous",
             "source_record_id": "record_000003", "experiment": "condition_C"},
        ], index=["7", "7_1", "12_1", "20"])

    def run_m4(self, data=None, scoring=True, assessment=True):
        data = self.data if data is None else data
        data.to_csv(self.root / "input.sqcs_structdi.csv")
        ops = data.copy()
        ops.index = [int(str(index).split("_", 1)[0]) for index in ops.index]
        ops = ops.loc[~ops.index.duplicated(), [
            "source_record_id", "candidate_id", "unip_id_a", "unip_id_b",
            "seq_a", "seq_b", "pep_a", "pep_b", "pos_a", "pos_b",
            "res_pos_a", "res_pos_b",
        ]]
        ops["homo_adjacency"] = 0.1
        ops["homo_int_overl"] = 0.0
        ops["homo_pep_overl"] = False
        ops.to_csv(self.root / "input.sqcs_ops.csv")
        policy = ([
            "--euclidean-assessment", "0", "35", "inclusive",
            "--topological-assessment", "0", "35", "inclusive",
        ] if assessment else [])
        with patch(
            "claudio.module04.src.algorithm.retrieve_oligo_state."
            "query_oligo_states_from_swiss", return_value={},
        ), patch.object(m4, "setup_pml_scripts") as pymol, patch.object(
            m4, "create_plots",
        ) as plots:
            result = CliRunner().invoke(m4.main, [
                "-i", str(self.root / "input.sqcs_structdi.csv"),
                "-i2", str(self.root / "input.sqcs_ops.csv"),
                "-o", str(self.root / "out"), "-s", str(scoring), "-v", "0", *policy,
            ])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        pymol.assert_not_called()
        plots.assert_not_called()
        return {
            suffix: pd.read_csv(self.root / "out" / f"input_{suffix}.csv", index_col=0)
            for suffix in ("extended", "final", "candidate_audit")
        }

    def test_pdb_candidates_keep_relationships_and_violations_separate(self):
        outputs = self.run_m4()
        for saved in outputs.values():
            self.assertEqual(len(saved), 4)
            self.assertEqual(saved["XL_type"].tolist(),
                             ["intra", "inter", "inter", "undetermined"])
            self.assertTrue(saved["XL_type_scope"].eq("model_candidate").all())
            self.assertEqual(saved["topo_assessment_status"].tolist(),
                             ["violated", "satisfied", "satisfied", "not_evaluable"])
            self.assertEqual(saved["experiment"].tolist(),
                             ["condition_A", "condition_A", "condition_B", "condition_C"])
            self.assertEqual(saved["source_crosslink_type"].tolist(), ["Homodimer"] * 4)
            self.assertIn("legacy_XL_type", saved)
            self.assertIn("legacy_score_XL_type", saved)
            self.assertNotIn("score_XL_type", saved)
            self.assertNotIn("XL_confirmed", saved)
        self.assertTrue(outputs["candidate_audit"]["m4_selection_status"].eq("retained").all())
        self.assertFalse((self.root / "out" / "homomers").exists())
        self.assertFalse((self.root / "out" / "heteromers").exists())

    def test_relationship_labels_do_not_depend_on_scoring_or_assessment_switches(self):
        for scoring, assessment in ((False, False), (True, False), (False, True)):
            with self.subTest(scoring=scoring, assessment=assessment):
                saved = self.run_m4(scoring=scoring, assessment=assessment)["final"]
                self.assertEqual(saved["XL_type"].tolist(),
                                 ["intra", "inter", "inter", "undetermined"])

    def test_empty_ordinary_and_af_inputs_have_classification_headers(self):
        for af in (False, True):
            with self.subTest(af=af):
                empty = self.data.iloc[:0].copy()
                if af:
                    empty["af_match_status"] = "matched"
                outputs = self.run_m4(empty)
                for saved in outputs.values():
                    self.assertTrue(saved.empty)
                    for column in ("XL_type", "XL_type_scope", "XL_type_reason",
                                   "legacy_XL_type", "topo_assessment_status"):
                        self.assertIn(column, saved)

    def test_rerun_refreshes_primary_relationship_without_changing_reported_type(self):
        self.run_m4()
        changed = self.data.copy()
        changed.loc["7", "chain_b"] = "B"
        saved = self.run_m4(changed)["final"]
        self.assertEqual(saved.loc["7", "XL_type"], "inter")
        self.assertEqual(saved.loc["7", "source_crosslink_type"], "Homodimer")
        self.assertEqual(len(saved), 4)


if __name__ == "__main__":
    unittest.main()
