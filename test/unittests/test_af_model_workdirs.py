"""Prevent parallel AlphaFold models from sharing TopoLink work files."""
import subprocess
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFIO, PDBParser
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from claudio.module02.src_distance_reevaluation.algorithm import (
    calc_site_distances,
)


class TestAFModelWorkdirs(unittest.TestCase):
    @staticmethod
    def write_model(path, distance):
        structure = Structure("complex")
        model = Model(0)
        chain = Chain("A")
        structure.add(model)
        model.add(chain)

        for number, x in ((1, 0.0), (2, distance)):
            residue = Residue((" ", number, " "), "LYS", " ")
            residue.add(
                Atom(
                    "CB",
                    np.array([x, 0.0, 0.0], dtype=float),
                    90.0,
                    1.0,
                    " ",
                    " CB ",
                    number,
                    element="C",
                )
            )
            chain.add(residue)

        writer = MMCIFIO()
        writer.set_structure(structure)
        writer.save(str(path))

    def test_parallel_models_keep_their_own_distances(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            rows = []
            for model_id, distance in (("AF-0001", 10.0), ("AF-0002", 30.0)):
                path = root / f"{model_id}-model_v1.cif"
                self.write_model(path, distance)
                rows.append(
                    {
                        "path": str(path),
                        "unip_id_a": "P1",
                        "unip_id_b": "P1",
                        "seq_a": "KK",
                        "seq_b": "KK",
                        "pos_a": 1,
                        "pos_b": 2,
                        "pdb_pos_a": 1,
                        "pdb_pos_b": 2,
                        "chain_a": "A",
                        "chain_b": "A",
                        "method_a": "alphafold",
                        "pLDDT_a": 90.0,
                        "pLDDT_b": 90.0,
                        "res_criteria_fulfilled": True,
                        "is_interfaced": True,
                    }
                )

            ready = threading.Barrier(2)
            errors = []
            input_paths = []
            structure_paths = []
            lock = threading.Lock()

            def fake_topolink(command, **kwargs):
                try:
                    input_path = Path(kwargs["cwd"]) / command[-1]

                    # Both workers finish writing before either reads.
                    ready.wait(timeout=10)

                    lines = input_path.read_text(encoding="utf-8").splitlines()
                    structure_line = next(
                        line for line in lines if line.startswith("structure ")
                    )
                    structure_path = Path(structure_line.split(maxsplit=1)[1])

                    exported = PDBParser(QUIET=True).get_structure(
                        "exported", (Path(kwargs["cwd"]) / (str(structure_path)))
                    )
                    chain = exported[0]["A"]
                    distance = float(chain[1]["CB"] - chain[2]["CB"])

                    with lock:
                        input_paths.append(input_path)
                        structure_paths.append((Path(kwargs["cwd"]) / structure_path).resolve())

                    # Controlled output at the external-program boundary.
                    return subprocess.CompletedProcess(
                        command, 0,
                        stdout=(
                            f"  LINK: LYS A 1 CB LYS A 2 CB "
                            f"{distance:.3f} {distance + 2.0:.3f}\n"
                        ),
                        stderr="",
                    )
                except Exception as error:
                    with lock:
                        errors.append(repr(error))
                    return subprocess.CompletedProcess(
                        command, 0, stdout="", stderr="",
                    )

            with patch.object(
                calc_site_distances.subprocess, "run", side_effect=fake_topolink
            ) as mocked_topolink:
                result = calc_site_distances.compute_dists_with_topolink(
                    pd.DataFrame(rows),
                    temp_dir=str(root / "work") + "/",
                    df_xl_res=pd.DataFrame([{"res": "K", "atom": "CB"}]),
                    plddt_cutoff=70.0,
                    topolink_bin=None,
                    verbose_level=0,
                )

            self.assertEqual(mocked_topolink.call_count, 2)
            self.assertEqual(errors, [], f"Test setup failed: {errors}")

            np.testing.assert_allclose(
                result["eucl_dist_tplk"].to_numpy(dtype=float),
                [10.0, 30.0],
                err_msg="Parallel models must retain their own distances",
            )
            np.testing.assert_allclose(
                result["topo_dist_tplk"].to_numpy(dtype=float),
                [12.0, 32.0],
            )
            self.assertEqual(len(set(input_paths)), 2)
            self.assertEqual(len(set(structure_paths)), 2)


if __name__ == "__main__":
    unittest.main()
