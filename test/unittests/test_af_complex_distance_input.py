import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd
from Bio.PDB import PDBParser

from claudio.module02.src_distance_reevaluation.algorithm import (
    search_pdb_pos,
)
from claudio.module02.src_distance_reevaluation.algorithm.calc_site_distances import (
    isolate_pdb_chain,
)
from claudio.utils.structure_io import read_structure_models

FIXTURE = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "af_complex"
    / "AF-0000000208747029.cif"
)


class TestAFComplexDistanceInput(unittest.TestCase):
    def test_cache_reads_complex_chains_and_residues(self):
        search_pdb_pos.pdb_cache.clear()
        self.addCleanup(search_pdb_pos.pdb_cache.clear)

        path = str(FIXTURE)
        data = pd.DataFrame({"path": [path]})
        search_pdb_pos.populate_pdb_cache(data, 0)

        models = search_pdb_pos.pdb_cache.get(path, [])
        self.assertTrue(models, "AF mmCIF must produce a parsed model")

        model = models[0]
        self.assertEqual({chain.id for chain in model}, {"A", "B"})

        for chain_id, position in (("A", 67), ("B", 11)):
            residue = model[chain_id][position]
            self.assertEqual(residue.get_resname(), "LYS")
            self.assertIn("CB", residue)

    def test_isolation_exports_both_complex_chains(self):
        with TemporaryDirectory() as temp_dir:
            (Path(temp_dir) / "struct").mkdir()

            output_path = isolate_pdb_chain(
                str(FIXTURE),
                "AF-0000000208747029",
                temp_dir + "/",
                ["A", "B"],
            )

            structure = PDBParser(QUIET=True).get_structure(
                "exported", output_path
            )
            model = next(structure.get_models())
            self.assertEqual({chain.id for chain in model}, {"A", "B"})
    def test_reads_single_chain_pdb(self):
        with TemporaryDirectory() as temp_dir:
            (Path(temp_dir) / "struct").mkdir()

            # Export chain A as an ordinary single-chain PDB file.
            pdb_path = isolate_pdb_chain(
                str(FIXTURE),
                "single_chain",
                temp_dir + "/",
                ["A"],
            )

            models = read_structure_models(pdb_path)

            self.assertEqual(len(models), 1)
            self.assertEqual(
                {chain.id for chain in models[0]},
                {"A"},
            )
            residue = models[0]["A"][67]
            self.assertEqual(residue.get_resname(), "LYS")
            self.assertIn("CB", residue)