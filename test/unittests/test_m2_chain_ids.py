"""Keep complete chain identities during M2 residue lookup."""
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from claudio.module02.src_distance_reevaluation.algorithm import search_pdb_pos


class TestM2ChainIds(unittest.TestCase):
    def test_matches_complete_chain_ids_for_both_endpoints(self):
        model = Model(0)
        expected = {"AA": 11.0, "AB": 22.0, "A": 33.0}
        for serial, (chain_id, x) in enumerate(expected.items(), start=1):
            chain = Chain(chain_id)
            residue = Residue((" ", 1, " "), "LYS", " ")
            residue.add(
                Atom(
                    "CB", np.array([x, 0.0, 0.0]), 91.0, 1.0,
                    " ", " CB ", serial, element="C",
                )
            )
            chain.add(residue)
            model.add(chain)

        mapping = pd.DataFrame(columns=["PDB", "CHAIN", "SP_PRIMARY"])
        with patch.object(search_pdb_pos, "pdb_cache", {"memory.cif": [model]}):
            for site in ("a", "b"):
                for chain_id, expected_x in expected.items():
                    with self.subTest(site=site, chain=chain_id):
                        row = pd.Series({
                            "path": "memory.cif",
                            "pdb_id": "AFTEST",
                            f"chain_{site}": chain_id,
                            f"unip_id_{site}": "P1",
                            f"seq_{site}": "K",
                            f"pos_{site}": 1,
                        })
                        result = search_pdb_pos.compute_site_pos(
                            7, row, site, "inter", mapping, "", ["K"], 0,
                            return_residue_id=True,
                            xl_atom_map={"K": ("CB",)},
                        )
                        self.assertTrue(
                            result[1], f"Chain {chain_id} must be found exactly"
                        )
                        self.assertEqual(result[0], (" ", 1, " "))
                        np.testing.assert_allclose(
                            result[6], [expected_x, 0.0, 0.0]
                        )


if __name__ == "__main__":
    unittest.main()
