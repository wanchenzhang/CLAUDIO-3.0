"""Exercise explicit distance policies through M4 and its three output tables."""
import importlib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd
from click.testing import CliRunner


m4 = importlib.import_module("claudio.module04.src.main")


class TestM4DistanceAssessmentCLI(unittest.TestCase):
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
            "res_crit_a": True, "res_crit_b": True, "res_criteria_fulfilled": True,
            "eucl_dist": 15.0, "topo_dist": 20.0,
            "eucl_dist_raw": "15.000", "topo_dist_raw": "20.000",
        }
        self.policies = [
            "--euclidean-assessment", "5", "30", "inclusive",
            "--topological-assessment", "5", "22", "exclusive",
        ]

    def run_m4(self, data, label, options=None):
        input_path = self.root / "input.sqcs_structdi.csv"
        ops_path = self.root / "input.sqcs_ops.csv"
        data.to_csv(input_path)
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
        ops.to_csv(ops_path)
        output = self.root / label
        with patch(
            "claudio.module04.src.algorithm.retrieve_oligo_state."
            "query_oligo_states_from_swiss", return_value={},
        ), patch.object(m4, "setup_pml_scripts"), patch.object(m4, "create_plots"):
            result = CliRunner().invoke(m4.main, [
                "-i", str(input_path), "-i2", str(ops_path), "-o", str(output),
                "-s", "False", "-v", "0", *(options or []),
            ])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        return {
            suffix: pd.read_csv(
                output / f"input_{suffix}.csv", index_col=0,
                dtype={"eucl_dist_raw": str, "topo_dist_raw": str},
            )
            for suffix in ("extended", "final", "candidate_audit")
        }

    def test_ordinary_outputs_keep_all_assessed_candidates_and_conditions(self):
        data = pd.DataFrame([
            self.row,
            {**self.row, "chain_b": "C", "topo_dist": 25.0, "topo_dist_raw": "25.000"},
            {**self.row, "source_record_id": "record_000002",
             "experiment": "condition_B", "source_crosslink_type": "Homodimer"},
        ], index=["7", "7_1", "12_1"])
        baseline = self.run_m4(data, "baseline")
        assessed = self.run_m4(data, "assessed", self.policies)
        for suffix in baseline:
            with self.subTest(table=suffix):
                pd.testing.assert_frame_equal(
                    assessed[suffix][baseline[suffix].columns], baseline[suffix],
                )
                for prefix in ("eucl", "topo"):
                    for field in ("status", "reason", "stage", "maximum", "boundary",
                                  "plddt_cutoff", "rule_version"):
                        self.assertIn(f"{prefix}_assessment_{field}", assessed[suffix])
        audit = assessed["candidate_audit"]
        self.assertEqual(audit["topo_assessment_status"].tolist(),
                         ["satisfied", "violated", "satisfied"])
        self.assertEqual(audit["m4_selection_status"].tolist(),
                         ["retained", "retained", "retained"])
        self.assertEqual(audit.loc["12_1", "experiment"], "condition_B")
        self.assertEqual(audit.loc["12_1", "source_crosslink_type"], "Homodimer")
        self.assertTrue(audit["topo_assessment_maximum"].eq(22.0).all())
        self.assertTrue(audit["topo_assessment_boundary"].eq("exclusive").all())
        self.assertEqual(assessed["extended"].loc["7_1", "topo_assessment_status"], "violated")
        self.assertEqual(assessed["extended"].loc["7_1", "XL_type"], "inter")

    def test_af_bounds_and_upstream_failures_reach_all_outputs(self):
        base = {
            **self.row, "af_match_status": "matched", "af_model_id": "AF-model",
            "af_model_version": 1, "pdb_id": "AF-model",
            "pdb_method": "ALPHAFOLD", "pLDDT_a": 90.0, "pLDDT_b": 90.0,
            "mapping_status_a": "mapped", "mapping_status_b": "mapped",
            "distance_status": "computed", "distance_reason": None,
        }
        data = pd.DataFrame([
            base,
            {**base, "distance_status": "lower_bound", "topo_dist": None,
             "topo_dist_lower_bound": 75.328, "topo_dist_raw": ">75.328"},
            {**base, "distance_status": "not_evaluable",
             "distance_reason": "identical_atom_endpoints"},
            {**base, "source_record_id": "record_000002", "experiment": "condition_B",
             "af_match_status": "no_candidate_model", "af_model_id": None,
             "af_model_version": None, "path": "-", "chain_a": "-", "chain_b": "-",
             "distance_status": "not_evaluable", "distance_reason": "no_candidate_model"},
        ], index=["7", "7_1", "7_2", "12"])
        outputs = self.run_m4(data, "af", self.policies)
        for suffix, saved in outputs.items():
            with self.subTest(table=suffix):
                self.assertEqual(len(saved), 4)
                self.assertEqual(saved["topo_assessment_status"].tolist(),
                                 ["satisfied", "violated", "not_evaluable", "not_evaluable"])
                self.assertEqual(saved["eucl_assessment_status"].tolist(),
                                 ["satisfied", "satisfied", "not_evaluable", "not_evaluable"])
                self.assertEqual(saved.loc["7_1", "topo_dist_raw"], ">75.328")
                self.assertTrue(pd.isna(saved.loc["7_1", "topo_dist"]))
                self.assertEqual(saved.loc["7_2", "topo_assessment_upstream_reason"],
                                 "identical_atom_endpoints")
                self.assertEqual(saved.loc["12", "topo_assessment_reason"], "no_candidate_model")
        self.assertTrue(outputs["candidate_audit"]["m4_selection_status"].eq("retained").all())

    def test_existing_plddt_option_controls_the_new_endpoint_gate(self):
        data = pd.DataFrame([{
            **self.row, "af_match_status": "matched", "af_model_id": "AF-model",
            "af_model_version": 1, "pdb_method": "ALPHAFOLD",
            "mapping_status_a": "mapped", "mapping_status_b": "mapped",
            "distance_status": "computed", "pLDDT_a": 90.0, "pLDDT_b": 75.0,
        }], index=["7"])
        outputs = self.run_m4(data, "confidence", [*self.policies, "-p", "80"])
        for saved in outputs.values():
            self.assertEqual(saved.iloc[0]["topo_assessment_reason"], "low_plddt")
            self.assertEqual(saved.iloc[0]["topo_assessment_plddt_cutoff"], 80.0)

    def test_legacy_af_monomer_outputs_include_the_new_assessment_columns(self):
        data = pd.DataFrame([{
            **self.row, "unip_id_b": "P1", "chain_b": "A",
            "pdb_id": "afP1", "pdb_method": "ALPHAFOLD",
            "pLDDT_a": 90.0, "pLDDT_b": 90.0,
        }], index=["7"])
        outputs = self.run_m4(data, "monomer", self.policies)
        for saved in outputs.values():
            self.assertEqual(saved.iloc[0]["structure_source"], "af_monomer")
            self.assertEqual(saved.iloc[0]["topo_assessment_status"], "satisfied")
            self.assertTrue(pd.isna(saved.iloc[0]["topo_assessment_upstream_distance_status"]))

    def test_empty_af_tables_still_include_policy_headers(self):
        data = pd.DataFrame([{
            **self.row, "af_match_status": "matched", "af_model_id": "AF-model",
            "af_model_version": 1,
        }]).iloc[:0]
        outputs = self.run_m4(data, "empty", self.policies)
        for saved in outputs.values():
            self.assertTrue(saved.empty)
            self.assertIn("eucl_assessment_status", saved)
            self.assertIn("topo_assessment_rule_version", saved)

    def test_omitting_policies_does_not_enable_distance_assessment(self):
        outputs = self.run_m4(pd.DataFrame([self.row], index=["7"]), "disabled")
        for saved in outputs.values():
            self.assertFalse(any(
                column.startswith(("eucl_assessment_", "topo_assessment_"))
                for column in saved
            ))

    def test_partial_or_invalid_policies_are_rejected_before_input_reads(self):
        invalid = (
            self.policies[:4],
            self.policies[4:],
            ["--euclidean-assessment", "-1", "30", "inclusive", *self.policies[4:]],
            ["--euclidean-assessment", "30", "5", "inclusive", *self.policies[4:]],
            ["--euclidean-assessment", "5", "nan", "inclusive", *self.policies[4:]],
            [*self.policies[:4], "--topological-assessment", "5", "5", "exclusive"],
        )
        for options in invalid:
            with self.subTest(options=options), patch.object(
                m4, "read_inputs", side_effect=AssertionError("Unexpected input read"),
            ), patch.object(
                m4, "create_out_path", side_effect=AssertionError("Unexpected output setup"),
            ):
                result = CliRunner().invoke(m4.main, options)
                self.assertEqual(result.exit_code, 2, result.output + repr(result.exception))
                self.assertIn("assessment", result.output)

    def test_rerunning_with_a_changed_policy_refreshes_output_results(self):
        data = pd.DataFrame([self.row], index=["7"])
        first = self.run_m4(data, "rerun", self.policies)
        self.assertEqual(first["final"].iloc[0]["topo_assessment_status"], "satisfied")
        changed = [*self.policies[:4], "--topological-assessment", "5", "15", "inclusive"]
        second = self.run_m4(data, "rerun", changed)
        for saved in second.values():
            self.assertEqual(saved.iloc[0]["topo_assessment_status"], "violated")
            self.assertEqual(saved.iloc[0]["topo_assessment_maximum"], 15.0)


if __name__ == "__main__":
    unittest.main()
