"""TopoLink failures must not look like successful runs without distances."""
import sys
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


class TestM2TopoLinkExecution(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.bin_dir = self.root / "tools with spaces"
        self.bin_dir.mkdir()
        self.executable = self.bin_dir / "topolink"
        model = Model(0)
        chain = Chain("A")
        model.add(chain)
        for number, x in ((1, 0.0), (2, 10.0)):
            residue = Residue((" ", number, " "), "LYS", " ")
            residue.add(Atom(
                "CA", np.array([x, 0.0, 0.0]), 90.0, 1.0,
                " ", " CA ", number, element="C",
            ))
            chain.add(residue)
        self.structure = self.root / "source model.pdb"
        writer = PDBIO()
        writer.set_structure(model)
        writer.save(str(self.structure))
        self.data = pd.DataFrame([{
            "source_record_id": "record_000001", "path": str(self.structure),
            "unip_id_a": "P1", "unip_id_b": "P1",
            "seq_a": "KK", "seq_b": "KK", "pos_a": 1, "pos_b": 2,
            "chain_a": "A", "chain_b": "A", "pdb_pos_a": 1, "pdb_pos_b": 2,
            "pdb_hetflag_a": " ", "pdb_hetflag_b": " ",
            "pdb_icode_a": " ", "pdb_icode_b": " ",
            "method_a": "alphafold", "pLDDT_a": 90.0, "pLDDT_b": 90.0,
            "res_criteria_fulfilled": True, "is_interfaced": True,
        }], index=[7])

    def write_executable(self, body):
        self.executable.write_text(f"#!{sys.executable}\n" + body)
        self.executable.chmod(0o755)

    def compute(self, work=None, bin_dir=None):
        return calc_site_distances.compute_dists_with_topolink(
            self.data.copy(), str(work or self.root / "work") + "/",
            pd.DataFrame([{"res": "K", "atom": "CA"}]),
            70.0, str(bin_dir or self.bin_dir) + "/", 0,
        )

    def test_long_work_paths_and_spaces_do_not_reach_native_path_buffers(self):
        self.write_executable(
            "import sys\nfrom pathlib import Path\n"
            "argument = sys.argv[1]\n"
            "if len('  TopoLink input file: ' + argument) > 200:\n"
            "    print('input path exceeds native buffer', file=sys.stderr)\n"
            "    sys.exit(23)\n"
            "lines = Path(argument).read_text().splitlines()\n"
            "structure = next(s.split(maxsplit=1)[1] for s in lines if s.startswith('structure '))\n"
            "if len('  Structure file: ' + structure) > 200:\n"
            "    sys.exit(24)\n"
            "assert Path(structure).is_file()\n"
            "print('  LINK: LYS A 1 CA LYS A 2 CA 10.000 12.000')\n"
        )
        work = self.root / ("deep output " * 5) / ("nested validation " * 5)
        cwd_before = Path.cwd()
        result = self.compute(work=work)
        self.assertEqual(Path.cwd(), cwd_before)
        self.assertEqual(result.index.tolist(), [7])
        self.assertEqual(result.loc[7, "source_record_id"], "record_000001")
        self.assertEqual(result.loc[7, "eucl_dist_tplk"], 10.0)
        self.assertEqual(result.loc[7, "topo_dist_tplk"], 12.0)

    def test_nonzero_exit_rejects_partial_stdout_and_keeps_stderr(self):
        self.write_executable(
            "import sys\n"
            "print('  LINK: LYS A 1 CA LYS A 2 CA 10.000 12.000')\n"
            "print('simulated native failure', file=sys.stderr)\n"
            "sys.exit(7)\n"
        )
        with self.assertRaisesRegex(RuntimeError, "TopoLink.*7"):
            self.compute()
        logs = list((self.root / "work").rglob("*.log"))
        self.assertEqual(len(logs), 1)
        text = logs[0].read_text()
        self.assertIn("simulated native failure", text)
        self.assertIn("10.000 12.000", text)
        self.assertIn(str(self.structure), text)

    def test_missing_executable_is_reported_with_a_log(self):
        missing = self.root / "missing tools"
        with self.assertRaisesRegex(RuntimeError, "TopoLink"):
            self.compute(bin_dir=missing)
        logs = list((self.root / "work").rglob("*.log"))
        self.assertEqual(len(logs), 1)
        self.assertIn(str(missing / "topolink"), logs[0].read_text())

    def test_successful_stderr_is_saved_without_losing_distances(self):
        self.write_executable(
            "import sys\n"
            "print('native diagnostic message', file=sys.stderr)\n"
            "print('  LINK: LYS A 1 CA LYS A 2 CA 10.000 12.000')\n"
        )
        result = self.compute()
        self.assertEqual(result.loc[7, "topo_dist_tplk"], 12.0)
        logs = list((self.root / "work").rglob("*.log"))
        self.assertEqual(len(logs), 1)
        self.assertIn("native diagnostic message", logs[0].read_text())

    def test_worker_structure_errors_reach_the_caller(self):
        with patch.object(
            calc_site_distances, "prepare_topolink_structure",
            side_effect=ValueError("bad structure in worker"),
        ):
            with self.assertRaisesRegex(ValueError, "bad structure in worker"):
                self.compute()


if __name__ == "__main__":
    unittest.main()
