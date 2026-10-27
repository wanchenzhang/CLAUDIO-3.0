"""Use the same chain-relation output for PDB and AlphaFold candidates."""
import unittest

import pandas as pd

from claudio.module04.src.algorithm.structure_assessment import (
    annotate_structure_relationships,
)
from claudio.module04.src.algorithm.af_complex_assessment import (
    annotate_complex_relationships,
)


class TestStructureRelationships(unittest.TestCase):
    def setUp(self):
        self.pdb = {
            "source_record_id": "record_000001",
            "source_crosslink_type": "Intra",
            "experiment": "condition_A",
            "unip_id_a": "P1", "unip_id_b": "P1",
            "pdb_id": "1ABC", "pdb_method": "X-RAY DIFFRACTION",
            "path": "/archive/1ABC.cif",
            "chain_a": "A", "chain_b": "A",
            "pdb_pos_a": 10, "pdb_pos_b": 20,
            "pLDDT_a": "-", "pLDDT_b": "-",
            "eucl_dist": 15.0, "topo_dist": 20.0,
        }
        self.monomer = {
            **self.pdb, "pdb_id": "afP1", "pdb_method": "ALPHAFOLD",
            "path": "/archive/afP1.pdb", "pLDDT_a": 20.0, "pLDDT_b": 30.0,
        }
        self.complex = {
            **self.pdb, "af_match_status": "matched",
            "af_model_id": "AF-0000000000000001", "af_model_version": 1,
            "path": "/archive/complex.cif", "chain_b": "B",
        }

    def test_pdb_chain_relation_does_not_reclassify_the_source(self):
        data = pd.DataFrame([
            self.pdb,
            {**self.pdb, "source_crosslink_type": "Homodimer", "chain_b": "B"},
            {**self.pdb, "source_crosslink_type": "Inter",
             "unip_id_b": "P2", "chain_b": "B"},
        ], index=[7, "7_1", 12])
        before = data.copy(deep=True)
        result = annotate_structure_relationships(data)
        pd.testing.assert_frame_equal(result[data.columns], before)
        pd.testing.assert_frame_equal(data, before)
        self.assertEqual(result["structure_source"].tolist(), ["pdb"] * 3)
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["same_chain", "cross_chain", "cross_chain"])
        self.assertEqual(result["model_chain_relation_reason"].tolist(),
                         ["same_model_chain", "different_model_chains",
                          "different_model_chains"])

    def test_af_monomer_uses_its_own_structure_identity(self):
        data = pd.DataFrame([
            self.monomer,
            {**self.monomer, "pdb_method": None,
             "method_a": "alphafold", "method_b": "alphafold"},
        ])
        result = annotate_structure_relationships(data)
        self.assertEqual(result["structure_source"].tolist(),
                         ["af_monomer", "af_monomer"])
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["same_chain", "same_chain"])
        pd.testing.assert_frame_equal(result[data.columns], data)

    def test_af_complex_does_not_require_a_pdb_accession(self):
        row = {key: value for key, value in self.complex.items()
               if key not in ("pdb_id", "pdb_method")}
        result = annotate_structure_relationships(pd.DataFrame([row]))
        self.assertEqual(result.iloc[0]["structure_source"], "af_complex")
        self.assertEqual(result.iloc[0]["model_chain_relation"], "cross_chain")

    def test_mixed_tables_are_dispatched_per_row(self):
        data = pd.DataFrame([self.pdb, self.monomer, self.complex])
        result = annotate_structure_relationships(data)
        self.assertEqual(result["structure_source"].tolist(),
                         ["pdb", "af_monomer", "af_complex"])
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["same_chain", "same_chain", "cross_chain"])
        pd.testing.assert_frame_equal(result[data.columns], data)

    def test_failed_af_assignment_is_not_reinterpreted_as_a_pdb_hit(self):
        for status in ("no_candidate_model", "no_matching_chain_pair"):
            with self.subTest(status=status):
                row = {**self.complex, "af_match_status": status}
                result = annotate_structure_relationships(pd.DataFrame([row]))
                self.assertEqual(result.iloc[0]["structure_source"], "af_complex")
                self.assertEqual(result.iloc[0]["model_chain_relation"], "unknown")
                self.assertEqual(result.iloc[0]["model_chain_relation_reason"], status)

    def test_missing_pdb_identity_or_path_is_unknown(self):
        for changes, reason in (
            ({"pdb_id": None}, "missing_model_identity"),
            ({"path": "-"}, "no_usable_structure"),
        ):
            with self.subTest(changes=changes):
                result = annotate_structure_relationships(
                    pd.DataFrame([{**self.pdb, **changes}])
                )
                self.assertEqual(result.iloc[0]["model_chain_relation"], "unknown")
                self.assertEqual(result.iloc[0]["model_chain_relation_reason"], reason)

    def test_complete_chain_ids_and_missing_assignments(self):
        data = pd.DataFrame([
            {**self.pdb, "chain_a": "AA", "chain_b": "AB"},
            {**self.pdb, "chain_a": "AA", "chain_b": "AA"},
            {**self.pdb, "chain_a": "-", "chain_b": "-"},
            {**self.pdb, "chain_b": pd.NA},
        ])
        result = annotate_structure_relationships(data)
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["cross_chain", "same_chain", "unknown", "unknown"])
        self.assertEqual(result["model_chain_relation_reason"].iloc[2:].tolist(),
                         ["missing_chain_assignment", "missing_chain_assignment"])
        pd.testing.assert_frame_equal(result[data.columns], data)

    def test_known_chains_are_independent_of_residue_or_distance_exclusion(self):
        data = pd.DataFrame([
            {**self.pdb, "pdb_pos_a": None, "distance_status": "not_evaluable"},
            {**self.monomer, "distance_status": "not_evaluable",
             "distance_reason": "low_plddt"},
            {**self.complex, "position_status_a": "ambiguous",
             "distance_status": "not_evaluable", "distance_reason": "unresolved_position"},
        ])
        result = annotate_structure_relationships(data)
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["same_chain", "same_chain", "cross_chain"])
        pd.testing.assert_frame_equal(result[data.columns], data)

    def test_af_wrapper_keeps_its_explicit_contract(self):
        data = pd.DataFrame([self.pdb])
        general = annotate_structure_relationships(data)
        af_only = annotate_complex_relationships(data)
        self.assertEqual(general.iloc[0]["model_chain_relation"], "same_chain")
        self.assertEqual(af_only.iloc[0]["model_chain_relation"], "unknown")
        self.assertEqual(af_only.iloc[0]["model_chain_relation_reason"],
                         "missing_af_match_status")
        complex_data = pd.DataFrame([self.complex])
        expected = annotate_structure_relationships(complex_data)
        wrapped = annotate_complex_relationships(complex_data)
        for column in ("structure_source", "model_chain_relation",
                       "model_chain_relation_reason"):
            pd.testing.assert_series_equal(wrapped[column], expected[column])
        self.assertEqual(wrapped.iloc[0]["chain_relation_rule_version"], "af-chain-v1")

    def test_reannotation_preserves_duplicate_indices_and_source_fields(self):
        data = pd.DataFrame([self.pdb, self.complex], index=[7, 7])
        data["structure_source"] = "stale"
        data["model_chain_relation"] = "stale"
        data["model_chain_relation_reason"] = "stale"
        data["chain_relation_rule_version"] = "old"
        before = data.copy(deep=True)
        result = annotate_structure_relationships(data)
        self.assertEqual(result.index.tolist(), [7, 7])
        self.assertEqual(result["structure_source"].tolist(), ["pdb", "af_complex"])
        self.assertEqual(result["model_chain_relation"].tolist(),
                         ["same_chain", "cross_chain"])
        self.assertTrue(result["chain_relation_rule_version"].ne("old").all())
        derived = {"structure_source", "model_chain_relation",
                   "model_chain_relation_reason", "chain_relation_rule_version"}
        columns = [column for column in data if column not in derived]
        pd.testing.assert_frame_equal(result[columns], before[columns])
        pd.testing.assert_frame_equal(data, before)

    def test_missing_structure_metadata_does_not_invent_a_source_type(self):
        data = pd.DataFrame([{"source_record_id": "record_000001"}])
        result = annotate_structure_relationships(data)
        self.assertEqual(result.iloc[0]["structure_source"], "unknown")
        self.assertEqual(result.iloc[0]["model_chain_relation"], "unknown")
        self.assertEqual(result.iloc[0]["model_chain_relation_reason"],
                         "unknown_structure_source")
        self.assertNotIn("source_crosslink_type", result.columns)
        self.assertNotIn("XL_type", result.columns)

    def test_invalid_explicit_source_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "structure_source"):
            annotate_structure_relationships(pd.DataFrame([self.pdb]),
                                             structure_source="unsupported")

    def test_empty_input_keeps_the_common_output_headers(self):
        data = pd.DataFrame([self.pdb]).iloc[:0]
        result = annotate_structure_relationships(data)
        self.assertTrue(result.empty)
        for column in ("structure_source", "model_chain_relation",
                       "model_chain_relation_reason", "chain_relation_rule_version"):
            self.assertIn(column, result.columns)
        pd.testing.assert_frame_equal(result[data.columns], data)


if __name__ == "__main__":
    unittest.main()
