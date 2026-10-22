"""Keep TopoLink distance bounds separate from numeric path results."""

import io
import shlex
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from claudio.module02.src_distance_reevaluation.algorithm import calc_site_distances


class TestM2DistanceBounds(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        model = Model(0)
        chain = Chain("A")
        model.add(chain)
        for number, x in ((1, 0.0), (2, 10.0), (3, 30.0)):
            residue = Residue((" ", number, " "), "LYS", " ")
            residue.add(Atom(
                "CB", np.array([x, 0.0, 0.0]), 90.0, 1.0,
                " ", " CB ", number, element="C",
            ))
            chain.add(residue)
        path = self.root / "structure.pdb"
        writer = PDBIO()
        writer.set_structure(model)
        writer.save(str(path))
        base = {
            "path": str(path), "unip_id_a": "P1", "unip_id_b": "P1",
            "seq_a": "KKK", "seq_b": "KKK", "chain_a": "A", "chain_b": "A",
            "pdb_hetflag_a": " ", "pdb_hetflag_b": " ",
            "pdb_icode_a": " ", "pdb_icode_b": " ",
            "method_a": "realigning", "pLDDT_a": 90.0, "pLDDT_b": 90.0,
            "res_criteria_fulfilled": True, "is_interfaced": True,
        }
        self.data = pd.DataFrame([
            {
                **base, "source_record_id": f"record_{number:06d}",
                "pos_a": a, "pos_b": b, "pdb_pos_a": a, "pdb_pos_b": b,
            }
            for number, (a, b) in enumerate(
                ((1, 2), (1, 2), (2, 1), (2, 3)), start=1,
            )
        ], index=[7, 12, 19, 23])

    def compute(self, euclidean="10.000", topological="12.000", data=None, missing=False):
        def topolink(command):
            lines = Path(shlex.split(command)[-1]).read_text().splitlines()
            observed = [line for line in lines if line.strip().startswith("observed ")]
            self.assertEqual(len(observed), 2)
            if missing:
                return io.StringIO("")
            return io.StringIO(
                f"  LINK: LYS A 1 CB LYS A 2 CB {euclidean} {topological}\n"
                "  LINK: LYS A 2 CB LYS A 3 CB 20.000 22.000\n"
            )

        with patch.object(calc_site_distances.os, "popen", side_effect=topolink):
            return calc_site_distances.compute_dists_with_topolink(
                (self.data if data is None else data).copy(),
                str(self.root / "work") + "/",
                pd.DataFrame([{"res": "K", "atom": "CB"}]), 70.0, None, 0,
            )

    def test_bounded_values_are_backfilled_without_becoming_exact_distances(self):
        cases = [
            ("10.000", ">12.3456", 10.0, None, None, 12.3456),
            (">10.000", "12.000", None, 12.0, 10.0, None),
            (">10.000", ">12.000", None, None, 10.0, 12.0),
            ("10.000", "12.000", 10.0, 12.0, None, None),
        ]
        for euclidean, topological, exact_e, exact_t, bound_e, bound_t in cases:
            with self.subTest(euclidean=euclidean, topological=topological):
                result = self.compute(euclidean, topological)
                self.assertEqual(result.index.tolist(), self.data.index.tolist())
                self.assertEqual(result.source_record_id.tolist(),
                                 self.data.source_record_id.tolist())
                for column, expected in (
                    ("eucl_dist_tplk", exact_e), ("topo_dist_tplk", exact_t),
                    ("eucl_dist_lower_bound", bound_e),
                    ("topo_dist_lower_bound", bound_t),
                ):
                    values = result.loc[[7, 12, 19], column]
                    if expected is None:
                        self.assertTrue(values.isna().all(), column)
                    else:
                        self.assertTrue(values.eq(expected).all(), column)
                self.assertTrue(result.loc[[7, 12, 19], "eucl_dist_raw"].eq(euclidean).all())
                self.assertTrue(result.loc[[7, 12, 19], "topo_dist_raw"].eq(topological).all())
                self.assertEqual(result.loc[23, "topo_dist_tplk"], 22.0)
                self.assertTrue(pd.isna(result.loc[23, "topo_dist_lower_bound"]))

    def test_rerun_replaces_old_bounds_and_raw_tokens(self):
        bounded = self.compute(topological=">12.000")
        result = self.compute(data=bounded)
        self.assertTrue(result["topo_dist_lower_bound"].isna().all())
        self.assertTrue(result.loc[[7, 12, 19], "topo_dist_tplk"].eq(12.0).all())
        self.assertTrue(result.loc[[7, 12, 19], "topo_dist_raw"].eq("12.000").all())

    def test_missing_results_clear_previous_values(self):
        bounded = self.compute(topological=">12.000")
        result = self.compute(data=bounded, missing=True)
        for column in ("eucl_dist_tplk", "topo_dist_tplk",
                       "eucl_dist_lower_bound", "topo_dist_lower_bound",
                       "eucl_dist_raw", "topo_dist_raw"):
            self.assertTrue(result[column].isna().all(), column)

    def test_invalid_tokens_are_preserved_but_not_used_as_distances(self):
        for token in ("NaN", "inf", "-1", "<12", "missing"):
            with self.subTest(token=token):
                result = self.compute(topological=token)
                self.assertTrue(result.loc[[7, 12, 19], "topo_dist_tplk"].isna().all())
                self.assertTrue(result.loc[[7, 12, 19], "topo_dist_lower_bound"].isna().all())
                self.assertTrue(result.loc[[7, 12, 19], "topo_dist_raw"].eq(token).all())


if __name__ == "__main__":
    unittest.main()
