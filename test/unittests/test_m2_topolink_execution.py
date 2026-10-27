"""TopoLink failures must not look like successful runs without distances."""
import hashlib
import json
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

    def compute(self, work=None, bin_dir=None, **options):
        return calc_site_distances.compute_dists_with_topolink(
            self.data.copy(), str(work or self.root / "work") + "/",
            pd.DataFrame([{"res": "K", "atom": "CA"}]),
            70.0, str(bin_dir or self.bin_dir) + "/", 0, **options,
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
        manifest = json.loads(logs[0].with_name("topolink_run.json").read_text())
        self.assertEqual(manifest["state"], "failed_to_start")
        self.assertIsNone(manifest["return_code"])
        self.assertIsNone(manifest["executable"]["sha256"])

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

    def test_default_seed_and_effective_native_configuration_are_recorded(self):
        self.write_executable(
            "print('Version 19.156')\n"
            "print('  LINK: LYS A 1 CA LYS A 2 CA 10.000 12.000')\n"
        )
        self.compute()
        files = list((self.root / "work").rglob("topolink_run.json"))
        self.assertEqual(len(files), 1)
        manifest = json.loads(files[0].read_text())
        self.assertIsNone(manifest["requested_seed"])
        self.assertEqual(manifest["effective_seed"], "random")
        self.assertEqual(manifest["state"], "completed")
        self.assertEqual(manifest["return_code"], 0)
        self.assertEqual(manifest["executable"]["version"], "19.156")
        self.assertEqual(manifest["executable"]["sha256"],
                         hashlib.sha256(self.executable.read_bytes()).hexdigest())
        self.assertEqual(manifest["source_structure"]["sha256"],
                         hashlib.sha256(self.structure.read_bytes()).hexdigest())
        self.assertEqual(manifest["parameters"]["plddt_cutoff"], 70.0)
        self.assertEqual(manifest["parameters"]["xl_residues"],
                         [{"res": "K", "atom": "CA"}])
        task = files[0].parent
        native_input = task / manifest["input"]["path"]
        self.assertEqual(manifest["input"]["sha256"],
                         hashlib.sha256(native_input.read_bytes()).hexdigest())
        self.assertEqual(manifest["input"]["text"], native_input.read_text())
        exported = task / manifest["exported_structure"]["path"]
        self.assertEqual(manifest["exported_structure"]["sha256"],
                         hashlib.sha256(exported.read_bytes()).hexdigest())
        self.assertEqual(manifest["command"], [str(self.executable), "in/topolink.tmp"])
        self.assertEqual(Path(manifest["cwd"]), task)
        self.assertTrue(manifest["started_at"].endswith("+00:00"))
        self.assertTrue(manifest["finished_at"].endswith("+00:00"))

    def test_default_preserves_a_seed_already_set_in_the_template(self):
        from io import StringIO

        template = Path(calc_site_distances.__file__).resolve().parents[3] / "data/topolink_inputfile.inp"
        lines = template.read_text().splitlines(keepends=True)
        text = "".join("seed 81\n" if line.split()[:1] == ["seed"] else line
                       for line in lines)
        self.write_executable("print('  LINK: LYS A 1 CA LYS A 2 CA 10.000 12.000')\n")
        real_open = open

        def read_template(path, *args, **kwargs):
            if Path(path) == template:
                return StringIO(text)
            return real_open(path, *args, **kwargs)

        with patch("builtins.open", side_effect=read_template):
            self.compute()
        path = next((self.root / "work").rglob("topolink_run.json"))
        manifest = json.loads(path.read_text())
        self.assertIsNone(manifest["requested_seed"])
        self.assertEqual(manifest["effective_seed"], "81")
        self.assertIsNone(manifest["executable"]["version"])

    def test_explicit_seed_reaches_native_input_and_manifest(self):
        self.write_executable(
            "import sys\nfrom pathlib import Path\n"
            "lines = Path(sys.argv[1]).read_text().splitlines()\n"
            "assert [s.strip() for s in lines if s.split()[:1] == ['seed']] == ['seed 20261026']\n"
            "print('Version 19.156')\n"
            "print('  LINK: LYS A 1 CA LYS A 2 CA 10.000 12.000')\n"
        )
        result = self.compute(topolink_seed=20261026)
        self.assertEqual(result.loc[7, "topo_dist_tplk"], 12.0)
        path = next((self.root / "work").rglob("topolink_run.json"))
        manifest = json.loads(path.read_text())
        self.assertEqual(manifest["requested_seed"], 20261026)
        self.assertEqual(manifest["effective_seed"], "20261026")

    def test_invalid_seed_is_rejected_before_creating_work_files(self):
        for value in (0, -1, 2147483648, True, 1.5, "42"):
            with self.subTest(seed=value):
                with self.assertRaisesRegex(ValueError, "seed"):
                    self.compute(topolink_seed=value)
        self.assertFalse((self.root / "work").exists())

    def test_explicit_seed_requires_one_template_seed_directive(self):
        template = Path(calc_site_distances.__file__).resolve().parents[3] / "data/topolink_inputfile.inp"
        original = template.read_text()
        for replacement in ("", "seed random\nseed random\n"):
            with self.subTest(template_seed=replacement):
                lines = [line for line in original.splitlines(keepends=True)
                         if line.split()[:1] != ["seed"]]
                lines.append("\n" + replacement)
                real_open = open

                def read_template(path, *args, **kwargs):
                    if Path(path) == template:
                        from io import StringIO
                        return StringIO("".join(lines))
                    return real_open(path, *args, **kwargs)

                with patch("builtins.open", side_effect=read_template):
                    with patch.object(calc_site_distances.subprocess, "run") as native:
                        with self.assertRaisesRegex(ValueError, "seed"):
                            self.compute(topolink_seed=42)
                        native.assert_not_called()

    def test_failed_native_run_keeps_its_requested_seed_and_configuration(self):
        self.write_executable(
            "import sys\n"
            "print('Version 19.156')\n"
            "print('seeded run failed', file=sys.stderr)\n"
            "sys.exit(7)\n"
        )
        with self.assertRaisesRegex(RuntimeError, "TopoLink.*7"):
            self.compute(topolink_seed=20261026)
        path = next((self.root / "work").rglob("topolink_run.json"))
        manifest = json.loads(path.read_text())
        self.assertEqual(manifest["state"], "failed")
        self.assertEqual(manifest["return_code"], 7)
        self.assertEqual(manifest["effective_seed"], "20261026")
        self.assertIn("seeded run failed", Path(manifest["log_path"]).read_text())

    def test_worker_structure_errors_reach_the_caller(self):
        with patch.object(
            calc_site_distances, "prepare_topolink_structure",
            side_effect=ValueError("bad structure in worker"),
        ):
            with self.assertRaisesRegex(ValueError, "bad structure in worker"):
                self.compute()


if __name__ == "__main__":
    unittest.main()
