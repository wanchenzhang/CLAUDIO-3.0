"""Exercise the candidate audit through the M4 command."""
import importlib
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd
from click.testing import CliRunner


m4 = importlib.import_module("claudio.module04.src.main")


class TestM4CandidateAuditCLI(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.row = {
            "source_record_id": "record_000001",
            "candidate_id": "candidate_000001", "experiment": "condition_A",
            "source_crosslink_type": "Inter",
            "unip_id_a": "P1", "unip_id_b": "P2",
            "seq_a": "MAAKTGGK", "seq_b": "MAAKTGGK",
            "pep_a": "AAKT", "pep_b": "GGK",
            "pos_a": 4, "pos_b": 8, "res_pos_a": 3, "res_pos_b": 3,
            "original_pos_a": 4, "original_pos_b": 8,
            "position_status_a": "valid", "position_status_b": "valid",
            "pep_a_repeated": False, "pep_b_repeated": False,
            "pep_a_repeat_count": 1, "pep_b_repeat_count": 1,
            "path": "/archive/1ABC.cif", "pdb_id": "1ABC",
            "pdb_method": "X-RAY DIFFRACTION", "pdb_resolution": 2.0,
            "chain_a": "A", "chain_b": "B", "pdb_pos_a": 4, "pdb_pos_b": 8,
            "pLDDT_a": "-", "pLDDT_b": "-",
            "res_crit_a": True, "res_crit_b": True,
            "eucl_dist": 15.0, "topo_dist": 20.0,
            "eucl_dist_raw": "15.000", "topo_dist_raw": "20.000",
            "is_interfaced": True, "method_a": "realigning",
        }

    def run_m4(self, structural, stem="input"):
        structural.to_csv(self.root / f"{stem}.sqcs_structdi.csv")
        ops = structural.copy()
        ops.index = [int(str(index).split("_", 1)[0]) for index in ops.index]
        ops = ops.loc[~ops.index.duplicated()]
        ops = ops[[
            "source_record_id", "candidate_id", "unip_id_a", "unip_id_b",
            "seq_a", "seq_b", "pep_a", "pep_b", "pos_a", "pos_b",
            "res_pos_a", "res_pos_b",
        ]]
        ops["homo_adjacency"] = 0.1
        ops["homo_int_overl"] = 0.0
        ops["homo_pep_overl"] = False
        ops.to_csv(self.root / f"{stem}.sqcs_ops.csv")
        with patch(
            "claudio.module04.src.algorithm.retrieve_oligo_state."
            "query_oligo_states_from_swiss", return_value={},
        ), patch.object(m4, "setup_pml_scripts"), patch.object(m4, "create_plots"):
            result = CliRunner().invoke(m4.main, [
                "-i", str(self.root / f"{stem}.sqcs_structdi.csv"),
                "-i2", str(self.root / f"{stem}.sqcs_ops.csv"),
                "-o", str(self.root / "out"), "-s", "False", "-v", "0",
            ])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        path = self.root / "out" / f"{stem}_candidate_audit.csv"
        self.assertTrue(path.is_file(), "M4 must write the complete candidate audit")
        return pd.read_csv(path, index_col=0,
                           dtype={"eucl_dist_raw": str, "topo_dist_raw": str})

    def test_preserves_chain_candidates_and_conditions_without_legacy_selection(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "chain_b": "C", "topo_dist": 25.0, "topo_dist_raw": "25.000"},
            {**self.row, "source_record_id": "record_000002",
             "experiment": "condition_B", "source_crosslink_type": "Homodimer"},
        ], index=["7", "7_1", "12_1"])
        saved = self.run_m4(data)
        self.assertEqual(saved.index.astype(str).tolist(), ["7", "7_1", "12_1"])
        self.assertEqual(saved["m4_selection_status"].tolist(),
                         ["retained", "retained", "retained"])
        self.assertEqual(saved["experiment"].tolist(),
                         ["condition_A", "condition_A", "condition_B"])
        self.assertEqual(saved.loc["12_1", "source_crosslink_type"], "Homodimer")
        self.assertEqual(saved.loc["12_1", "source_record_id"], "record_000002")
        self.assertEqual(saved["topo_dist_raw"].tolist(), ["20.000", "25.000", "20.000"])
        self.assertTrue(all(json.loads(value) == [] for value in saved["m4_exclusion_events"]))
        self.assertIn("method_a", saved.columns)
        self.assertEqual(saved["model_chain_relation"].tolist(), ["cross_chain"] * 3)
        final = pd.read_csv(self.root / "out/input_final.csv", index_col=0)
        extended = pd.read_csv(self.root / "out/input_extended.csv", index_col=0)
        self.assertEqual(final.index.astype(str).tolist(), ["7", "7_1", "12_1"])
        self.assertEqual(extended.index.astype(str).tolist(), ["7", "7_1", "12_1"])
        self.assertNotIn("m4_selection_status", final.columns)

    def test_repeated_peptide_sources_keep_their_conditions(self):
        data = pd.DataFrame([
            {**self.row, "seq_a": "MAAKTAAKTGGK"},
            {**self.row, "seq_a": "MAAKTAAKTGGK", "chain_b": "C",
             "topo_dist": 25.0},
            {**self.row, "seq_a": "MAAKTAAKTGGK",
             "source_record_id": "record_000002", "experiment": "condition_B"},
        ], index=["7", "7_1", "12"])
        saved = self.run_m4(data)
        self.assertEqual(saved["m4_selection_status"].tolist(),
                         ["retained", "retained", "retained"])
        self.assertEqual(json.loads(saved.loc["12", "m4_exclusion_events"]), [])
        self.assertEqual(saved.loc["12", "experiment"], "condition_B")

    def test_af_candidates_and_unresolved_positions_are_retained(self):
        data = pd.DataFrame([
            {**self.row, "af_match_status": "matched", "af_model_id": "AF-model",
             "af_model_version": 1},
            {**self.row, "af_match_status": "matched", "af_model_id": "AF-model",
             "af_model_version": 1, "chain_b": "C"},
            {**self.row, "source_record_id": "record_000002",
             "experiment": "condition_B", "af_match_status": "no_candidate_model",
             "af_model_id": None, "af_model_version": None, "path": "-",
             "chain_a": "-", "chain_b": "-", "position_status_a": "ambiguous",
             "distance_status": "not_evaluable", "distance_reason": "no_candidate_model",
             "eucl_dist": None, "topo_dist": None},
        ], index=["7", "7_1", "12"])
        saved = self.run_m4(data)
        self.assertEqual(len(saved), 3)
        self.assertTrue(saved["m4_selection_status"].eq("retained").all())
        self.assertTrue(all(json.loads(value) == []
                            for value in saved["m4_exclusion_events"]))
        self.assertEqual(saved.loc["12", "distance_status"], "not_evaluable")
        self.assertEqual(saved.loc["12", "position_status_a"], "ambiguous")
        final = pd.read_csv(self.root / "out/input_final.csv", index_col=0)
        self.assertEqual(len(final), 3)

    def test_empty_af_input_writes_audit_headers(self):
        data = pd.DataFrame([{**self.row, "af_match_status": "matched",
                              "af_model_id": "AF-model", "af_model_version": 1}]).iloc[:0]
        saved = self.run_m4(data)
        self.assertTrue(saved.empty)
        for column in ("experiment", "source_record_id", "m4_selection_status",
                       "m4_exclusion_events", "m4_audit_rule_version"):
            self.assertIn(column, saved.columns)

    def test_rerun_replaces_an_old_audit_in_the_same_output_directory(self):
        data = pd.DataFrame([
            self.row, {**self.row, "chain_b": "C", "topo_dist": 25.0},
        ], index=["7", "7_1"])
        first = self.run_m4(data)
        self.assertEqual(first["m4_selection_status"].tolist(), ["retained", "retained"])
        second = self.run_m4(data.iloc[[0]])
        self.assertEqual(len(second), 1)
        self.assertEqual(second.iloc[0]["m4_selection_status"], "retained")
        self.assertEqual(json.loads(second.iloc[0]["m4_exclusion_events"]), [])


if __name__ == "__main__":
    unittest.main()
