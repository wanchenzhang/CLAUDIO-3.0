"""Forward explicit assessment policies through the full CLI."""
import importlib
import json
import unittest
from contextlib import ExitStack
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from click.testing import CliRunner


pipeline = importlib.import_module("claudio.claudio_m")
STAGES = (
    "db_updater", "run_claudio_lists", "run_claudio_ops",
    "run_claudio_structdi", "run_claudio_xl",
)


class TestDistanceAssessmentPipelineCLI(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.root = Path(folder.name)
        self.input = self.root / "input.csv"
        self.input.write_text("input\n")
        self.output = self.root / "out"
        self.policies = [
            "--euclidean-assessment", "0", "100", "inclusive",
            "--topological-assessment", "0", "70", "exclusive",
        ]
        self.files = {
            "run_claudio_lists": "input.sqcs",
            "run_claudio_ops": "input.sqcs_ops.csv",
            "run_claudio_structdi": "input.sqcs_structdi.csv",
            "run_claudio_xl": "input_final.csv",
        }
        self.calls = []

    def finish_stage(self, name):
        def finish(**params):
            self.calls.append({"stage": name, "params": params})
            if name in self.files:
                self.output.mkdir(parents=True, exist_ok=True)
                (self.output / self.files[name]).write_text("result\n")
            if name == "run_claudio_xl":
                (self.output / "received_parameters.json").write_text(
                    json.dumps(self.calls, indent=2)
                )
        return finish

    def invoke(self, options=()):
        self.calls = []
        with ExitStack() as stack:
            for name in STAGES:
                stack.enter_context(
                    patch.object(pipeline, name, self.finish_stage(name))
                )
            return CliRunner().invoke(pipeline.main, [
                "-i", str(self.input), "-o", str(self.output),
                "-bldb", str(self.root / "missing-db"),
                "-rt", "False", "-s", "False", "-v", "0", *options,
            ])

    def test_policies_reach_only_m4_in_legacy_manual_and_auto_modes(self):
        for mode in (
            [], ["--af-complex-model", "AF-test"], ["--af-complex-auto"],
        ):
            with self.subTest(mode=mode):
                result = self.invoke([
                    *mode, *self.policies, "--topolink-seed", "20261026",
                ])
                self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
                received = json.loads((self.output / "received_parameters.json").read_text())
                by_stage = {item["stage"]: item["params"] for item in received}
                self.assertEqual(
                    by_stage["run_claudio_xl"]["euclidean_assessment"],
                    [0.0, 100.0, "inclusive"],
                )
                self.assertEqual(
                    by_stage["run_claudio_xl"]["topological_assessment"],
                    [0.0, 70.0, "exclusive"],
                )
                for name, params in by_stage.items():
                    if name != "run_claudio_xl":
                        self.assertNotIn("euclidean_assessment", params)
                        self.assertNotIn("topological_assessment", params)
                self.assertEqual(by_stage["run_claudio_structdi"]["topolink_seed"], 20261026)
                self.assertTrue(by_stage["run_claudio_lists"]["preserve_source_records"])
                self.assertEqual("db_updater" in by_stage, not bool(mode))

    def test_omitted_policies_leave_the_existing_m4_call_unchanged(self):
        result = self.invoke(["--af-complex-auto"])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        received = json.loads((self.output / "received_parameters.json").read_text())
        m4 = next(item["params"] for item in received if item["stage"] == "run_claudio_xl")
        self.assertNotIn("euclidean_assessment", m4)
        self.assertNotIn("topological_assessment", m4)

    def test_invalid_policies_stop_before_any_pipeline_stage_or_output_setup(self):
        invalid = (
            self.policies[:4],
            self.policies[4:],
            ["--euclidean-assessment", "-1", "100", "inclusive", *self.policies[4:]],
            ["--euclidean-assessment", "100", "0", "inclusive", *self.policies[4:]],
            ["--euclidean-assessment", "0", "nan", "inclusive", *self.policies[4:]],
            [*self.policies[:4], "--topological-assessment", "5", "5", "exclusive"],
        )
        for options in invalid:
            with self.subTest(options=options):
                result = self.invoke(options)
                self.assertEqual(result.exit_code, 2, result.output + repr(result.exception))
                self.assertIn("assessment", result.output)
                self.assertEqual(self.calls, [])
                self.assertFalse(self.output.exists())

    def test_legacy_config_and_explicit_cli_policies_can_be_used_together(self):
        config = self.root / "legacy.cfg"
        config.write_text(
            "input_filepath=unused.csv\n"
            "input_temppath=None\n"
            "projections=pep_a,pep_b,pos_a,pos_b,res_pos_a,res_pos_b,unip_id_a,unip_id_b\n"
            "read_temps=False\nxl_residues=K:CA:0\nsearch_tool=blastp\n"
            "e_value=1e-5\nquery_id=90\ncoverage=50\nres_cutoff=4\n"
            "plddt_cutoff=75\nlinker_minimum=5\nlinker_maximum=35\n"
            "euclidean_strictness=None\ndistance_maximum=50\ncutoff=0\n"
            "output_directory=unused\nblast_bin=None\nblast_db=unused\n"
            "hhsearch_bin=None\nhhsearch_db=None\ntopolink_bin=None\n"
            "compute_scoring=False\nverbose_level=0\n"
        )
        result = self.invoke([
            "--af-complex-auto", "-c", str(config), *self.policies,
        ])
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        received = json.loads((self.output / "received_parameters.json").read_text())
        by_stage = {item["stage"]: item["params"] for item in received}
        self.assertEqual(float(by_stage["run_claudio_xl"]["plddt_cutoff"]), 75.0)
        self.assertEqual(by_stage["run_claudio_xl"]["euclidean_assessment"],
                         [0.0, 100.0, "inclusive"])
        self.assertEqual(Path(by_stage["run_claudio_lists"]["input_filepath"]), self.input)


if __name__ == "__main__":
    unittest.main()
