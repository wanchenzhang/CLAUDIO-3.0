"""AF candidate tables retain their rows and provenance through M4 output."""

import importlib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd
from click.testing import CliRunner

from claudio.module04.src.io.write_outs import write_outputs


m4 = importlib.import_module("claudio.module04.src.main")


class TestAFComplexM4Outputs(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        base = {
            "source_record_id": "record_000001",
            "candidate_id": "candidate_000001", "experiment": "condition_A",
            "unip_id_a": "P1", "unip_id_b": "P1",
            "seq_a": "MAAKTGGK", "seq_b": "MAAKTGGK",
            "pep_a": "AAKT", "pep_b": "GGK",
            "pos_a": 4, "pos_b": 8, "res_pos_a": 3, "res_pos_b": 3,
            "original_pos_a": 2, "original_pos_b": 8,
            "position_status_a": "corrected", "position_status_b": "valid",
            "pep_a_repeated": False, "pep_b_repeated": False,
            "pep_a_repeat_count": 1, "pep_b_repeat_count": 1,
            "af_match_status": "matched", "af_model_id": "AF-model-one",
            "af_model_version": 1, "af_metadata_path": "model-one.json",
            "path": "model-one.cif", "pdb_id": "AF-model-one",
            "pdb_method": "ALPHAFOLD", "pdb_resolution": "ALPHAFOLD",
            "chain_a": "A", "chain_b": "A",
            "pdb_pos_a": 104, "pdb_pos_b": 108,
            "pdb_hetflag_a": " ", "pdb_icode_a": "A",
            "pdb_hetflag_b": " ", "pdb_icode_b": " ",
            "mapping_status_a": "mapped", "mapping_status_b": "mapped",
            "res_crit_a": True, "res_crit_b": True,
            "mapping_reason_a": None, "mapping_reason_b": None,
            "mapped_atom_name_a": "CB", "mapped_atom_name_b": "CB",
            "mapped_residue_id_a": "(' ', 104, 'A')",
            "mapped_residue_id_b": "(' ', 108, ' ')",
            "pLDDT_a": 90.0, "pLDDT_b": 90.0,
            "eucl_dist": 15.0, "topo_dist": 20.0,
            "distance_status": "computed", "distance_reason": None,
        }
        second_source = {
            **base, "source_record_id": "record_000002",
            "experiment": "condition_B", "chain_b": "B",
            "eucl_dist": 20.0, "topo_dist": 25.0,
        }
        unresolved = {
            **base, "source_record_id": "record_000003",
            "candidate_id": None, "position_status_a": "ambiguous",
            "mapping_status_a": "not_evaluable",
            "mapping_reason_a": "unresolved_position",
            "distance_status": "not_evaluable",
            "distance_reason": "unresolved_position",
            "pdb_pos_a": None, "pLDDT_a": None,
            "eucl_dist": None, "topo_dist": None,
        }
        self.structural = pd.DataFrame([
            base,
            {**base, "chain_b": "B", "eucl_dist": 20.0, "topo_dist": 25.0},
            {**base, "chain_b": "C", "eucl_dist": 60.0, "topo_dist": 65.0},
            {
                **base, "chain_b": "B", "af_model_id": "AF-model-two",
                "pdb_id": "AF-model-two", "path": "model-two.cif",
                "af_metadata_path": "model-two.json",
                "eucl_dist": 22.0, "topo_dist": 26.0,
            },
            second_source,
            {**second_source, "chain_b": "C", "eucl_dist": 60.0, "topo_dist": 65.0},
            unresolved,
        ], index=["7", "7_1", "7_2", "7_3", "12", "12_1", "20"])
        self.ops = self.structural.loc[["7", "12", "20"], [
            "source_record_id", "candidate_id", "unip_id_a", "unip_id_b",
            "seq_a", "seq_b", "pep_a", "pep_b", "pos_a", "pos_b",
            "res_pos_a", "res_pos_b",
        ]].copy()
        self.ops.index = [7, 12, 20]
        self.ops["homo_adjacency"] = [0.1, 0.1, float("nan")]
        self.ops["homo_int_overl"] = [0.0, 0.0, float("nan")]
        self.ops["homo_pep_overl"] = [False, False, pd.NA]

    def run_m4(self, structural, ops=None, stem="input", scoring=False):
        input_path = self.root / f"{stem}.sqcs_structdi.csv"
        ops_path = self.root / f"{stem}.sqcs_ops.csv"
        structural.to_csv(input_path)
        (self.ops if ops is None else ops).to_csv(ops_path)
        output = self.root / f"out_{stem}_{scoring}"
        with patch(
            "claudio.module04.src.algorithm.retrieve_oligo_state."
            "query_oligo_states_from_swiss", return_value={},
        ), patch.object(m4, "setup_pml_scripts") as pymol, patch.object(
            m4, "create_plots",
        ) as plots:
            result = CliRunner().invoke(m4.main, [
                "-i", str(input_path), "-i2", str(ops_path),
                "-o", str(output), "-p", "0", "-s", str(scoring), "-v", "0",
            ])
        self.assertEqual(result.exit_code, 0,
                         result.output + repr(result.exception))
        return output, pymol, plots

    def test_extended_and_final_keep_all_candidates_and_source_fields(self):
        output, pymol, plots = self.run_m4(self.structural)
        expected = pd.read_csv(self.root / "input.sqcs_structdi.csv", index_col=0)
        for suffix in ("extended", "final"):
            with self.subTest(output=suffix):
                saved = pd.read_csv(output / f"input_{suffix}.csv", index_col=0)
                self.assertEqual(len(saved), 7)
                pd.testing.assert_frame_equal(saved[expected.columns], expected)
                unresolved = saved.loc["20"]
                self.assertEqual(unresolved["distance_reason"], "unresolved_position")
                self.assertTrue(pd.isna(unresolved["topo_dist"]))
                self.assertTrue(pd.isna(unresolved["homo_pep_overl"]))
        pymol.assert_not_called()
        plots.assert_not_called()
        self.assertFalse((output / "homomers").exists())
        self.assertFalse((output / "heteromers").exists())

    def test_af_writer_keeps_all_columns_without_creating_modeling_inputs(self):
        write_outputs(self.structural, "candidates", False, str(self.root) + "/")
        saved = pd.read_csv(self.root / "candidates_final.csv", index_col=0)
        expected_path = self.root / "expected.csv"
        self.structural.to_csv(expected_path)
        expected = pd.read_csv(expected_path, index_col=0)
        pd.testing.assert_frame_equal(saved, expected)
        self.assertFalse((self.root / "homomers").exists())
        self.assertFalse((self.root / "heteromers").exists())

    def test_no_model_record_keeps_missing_results_and_reason(self):
        data = self.structural.loc[["7"]].copy()
        for column in ("af_model_id", "af_model_version", "af_metadata_path",
                       "pdb_pos_a", "pdb_pos_b", "pLDDT_a", "pLDDT_b",
                       "eucl_dist", "topo_dist"):
            data[column] = None
        data["af_match_status"] = "no_candidate_model"
        data["mapping_status_a"] = "not_evaluable"
        data["mapping_status_b"] = "not_evaluable"
        data["distance_status"] = "not_evaluable"
        data["distance_reason"] = "no_candidate_model"
        for column in ("path", "pdb_id", "chain_a", "chain_b"):
            data[column] = "-"
        output, _, _ = self.run_m4(data)
        saved = pd.read_csv(output / "input_final.csv", index_col=0)
        self.assertEqual(len(saved), 1)
        row = saved.iloc[0]
        self.assertEqual(row["source_record_id"], "record_000001")
        self.assertEqual(row["distance_reason"], "no_candidate_model")
        self.assertTrue(pd.isna(row["eucl_dist"]))
        self.assertTrue(pd.isna(row["topo_dist"]))
        self.assertIn("structure_source", saved.columns)
        self.assertEqual(row["structure_source"], "af_complex")
        self.assertEqual(row["model_chain_relation"], "unknown")
        self.assertEqual(row["model_chain_relation_reason"], "no_candidate_model")

    def test_empty_af_input_writes_result_headers_with_or_without_scoring(self):
        for scoring in (False, True):
            with self.subTest(scoring=scoring):
                output, pymol, plots = self.run_m4(
                    self.structural.iloc[:0], self.ops.iloc[:0], scoring=scoring,
                )
                for suffix in ("extended", "final"):
                    saved = pd.read_csv(output / f"input_{suffix}.csv", index_col=0)
                    self.assertTrue(saved.empty)
                    for column in ("source_record_id", "af_model_id",
                                   "distance_status", "distance_reason",
                                   "homo_pep_overl", "XL_type", "legacy_XL_confirmed",
                                   "XL_type_scope", "XL_type_reason",
                                   "structure_source", "model_chain_relation",
                                   "model_chain_relation_reason",
                                   "chain_relation_rule_version"):
                        self.assertIn(column, saved.columns)
                pymol.assert_not_called()
                plots.assert_not_called()

    def test_non_af_input_keeps_candidates_without_legacy_selection(self):
        data = self.structural.loc[["7"]].drop(columns=[
            column for column in self.structural.columns if column.startswith("af_")
        ])
        with patch.object(
            m4, "clean_dataset", wraps=m4.clean_dataset,
        ) as clean, patch.object(
            m4, "minimize_dataset", wraps=m4.minimize_dataset,
        ) as minimize:
            output, pymol, plots = self.run_m4(data, self.ops.loc[[7]])
        clean.assert_not_called()
        minimize.assert_not_called()
        pymol.assert_not_called()
        plots.assert_not_called()
        self.assertTrue((output / "input_final.csv").is_file())

    def test_input_stem_with_dots_is_preserved_in_output_names(self):
        output, _, _ = self.run_m4(self.structural, stem="sample.run1")
        self.assertTrue((output / "sample.run1_extended.csv").is_file())
        self.assertTrue((output / "sample.run1_final.csv").is_file())


    def test_af_outputs_include_relationships_and_reported_classification(self):
        structural = self.structural.copy()
        labels = ["Intra", "Intra", "Intra", "Intra", "Homodimer", "Homodimer", "Intra"]
        structural["source_crosslink_type"] = labels
        output, _, _ = self.run_m4(structural, stem="af_relations")
        for suffix in ("extended", "final"):
            with self.subTest(output=suffix):
                saved = pd.read_csv(output / f"af_relations_{suffix}.csv", index_col=0)
                for column in ("structure_source", "model_chain_relation",
                               "model_chain_relation_reason", "chain_relation_rule_version"):
                    self.assertIn(column, saved.columns)
                self.assertEqual(saved["structure_source"].tolist(), ["af_complex"] * 7)
                self.assertEqual(saved["model_chain_relation"].tolist(),
                                 ["same_chain", "cross_chain", "cross_chain", "cross_chain",
                                  "cross_chain", "cross_chain", "same_chain"])
                self.assertEqual(saved["source_crosslink_type"].tolist(), labels)
                self.assertTrue(saved["chain_relation_rule_version"].eq("structure-chain-v1").all())
                self.assertEqual(saved["model_chain_relation_reason"].iloc[0], "same_model_chain")
                self.assertEqual(saved["model_chain_relation_reason"].iloc[1], "different_model_chains")
                self.assertEqual(saved["XL_type"].tolist(),
                                 ["intra", "inter", "inter", "inter", "inter", "inter", "undetermined"])
                self.assertEqual(saved["legacy_XL_type"].tolist(),
                                 ["intra", "inter", "inter", "inter", "inter", "inter", "intra"])
                self.assertTrue(saved["XL_type_scope"].eq("model_candidate").all())
                self.assertEqual(saved["distance_reason"].iloc[-1], "unresolved_position")
                self.assertEqual(saved["source_record_id"].tolist(), structural.source_record_id.tolist())

    def test_ordinary_pdb_and_af_monomer_outputs_include_relationships(self):
        cases = [
            ("pdb_same", "1ABC", "X-RAY DIFFRACTION", "A", "P1", "pdb", "same_chain", "Homodimer"),
            ("pdb_cross", "1ABC", "X-RAY DIFFRACTION", "B", "P2", "pdb", "cross_chain", "Inter"),
            ("af_monomer", "afP1", "ALPHAFOLD", "A", "P1", "af_monomer", "same_chain", "Intra"),
        ]
        for stem, model_id, method, chain_b, protein_b, source, relation, label in cases:
            with self.subTest(case=stem):
                structural = self.structural.loc[["7"]].drop(columns=[
                    column for column in self.structural if column.startswith("af_")
                ]).copy()
                structural["pdb_id"] = model_id
                structural["pdb_method"] = method
                structural["path"] = f"/archive/{model_id}.cif"
                structural["chain_b"] = chain_b
                structural["unip_id_b"] = protein_b
                structural["source_crosslink_type"] = label
                if source == "pdb":
                    structural["pLDDT_a"] = "-"
                    structural["pLDDT_b"] = "-"
                ops = self.ops.loc[[7]].copy()
                ops["unip_id_b"] = protein_b
                output, _, _ = self.run_m4(structural, ops, stem=stem)
                for suffix in ("extended", "final"):
                    saved = pd.read_csv(output / f"{stem}_{suffix}.csv", index_col=0)
                    self.assertEqual(len(saved), 1)
                    for column in ("structure_source", "model_chain_relation",
                                   "model_chain_relation_reason", "chain_relation_rule_version",
                                   "source_crosslink_type"):
                        self.assertIn(column, saved.columns)
                    self.assertEqual(saved.iloc[0]["structure_source"], source)
                    self.assertEqual(saved.iloc[0]["model_chain_relation"], relation)
                    self.assertEqual(saved.iloc[0]["source_crosslink_type"], label)
                    self.assertEqual(saved.iloc[0]["model_chain_relation_reason"],
                                     "same_model_chain" if relation == "same_chain" else "different_model_chains")
                    self.assertEqual(saved.iloc[0]["chain_relation_rule_version"], "structure-chain-v1")
                    self.assertEqual(saved.iloc[0]["eucl_dist"], 15.0)
                    self.assertEqual(saved.iloc[0]["topo_dist"], 20.0)

    def test_rerunning_m4_refreshes_stale_relationship_fields(self):
        structural = self.structural.loc[["7"]].copy()
        for column in ("structure_source", "model_chain_relation",
                       "model_chain_relation_reason", "chain_relation_rule_version"):
            structural[column] = "stale"
        output, _, _ = self.run_m4(structural, self.ops.loc[[7]], stem="refresh")
        for suffix in ("extended", "final"):
            saved = pd.read_csv(output / f"refresh_{suffix}.csv", index_col=0)
            self.assertEqual(saved.iloc[0]["structure_source"], "af_complex")
            self.assertEqual(saved.iloc[0]["model_chain_relation"], "same_chain")
            self.assertEqual(saved.iloc[0]["model_chain_relation_reason"], "same_model_chain")
            self.assertEqual(saved.iloc[0]["chain_relation_rule_version"], "structure-chain-v1")
            self.assertNotIn("source_crosslink_type", saved.columns)


if __name__ == "__main__":
    unittest.main()
