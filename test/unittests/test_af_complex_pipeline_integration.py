"""Run the full CLI with local AF geometry and controlled service responses."""

import importlib
import json
import shlex
import unittest
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from urllib.parse import parse_qs, urlsplit

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from click.testing import CliRunner
from requests import Response


pipeline = importlib.import_module("claudio.claudio_m")
FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "af_complex"
MODEL_ID = "AF-0000000208747029"


class TestAFComplexPipelineIntegration(unittest.TestCase):
    def setUp(self):
        folder = TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        self.addCleanup(plt.close, "all")
        self.root = Path(folder.name)
        self.metadata = json.loads((FIXTURES / "metadata.json").read_text())
        chains = {record["chainId"]: record for record in self.metadata}
        self.protein_a = chains["A"]["uniprotAccession"]
        self.protein_b = chains["B"]["uniprotAccession"]
        self.sequences = {
            record["uniprotAccession"]: record["sequence"]
            for record in self.metadata
        }
        sequence_a = chains["A"]["sequence"]
        sequence_b = chains["B"]["sequence"]
        forward = {
            "entry1": self.protein_a, "entry2": self.protein_b,
            "peptide1": sequence_a[63:70], "peptide2": sequence_b[7:14],
            "position1": 76, "position2": 11,
            "k_pos1": 4, "k_pos2": 4, "experiment": "condition_A",
        }
        reverse = {
            "entry1": self.protein_b, "entry2": self.protein_a,
            "peptide1": sequence_b[19:26], "peptide2": sequence_a[63:70],
            "position1": 23, "position2": 67,
            "k_pos1": 4, "k_pos2": 4, "experiment": "condition_B",
        }
        unresolved = {
            "entry1": self.protein_a, "entry2": self.protein_b,
            "peptide1": sequence_a[235:243], "peptide2": sequence_b[66:73],
            "position1": 100, "position2": 70,
            "k_pos1": None, "k_pos2": 4, "experiment": "condition_C",
        }
        self.input = self.root / "crosslinks.run1.csv"
        pd.DataFrame([forward, reverse, unresolved, forward]).to_csv(
            self.input, index=False,
        )
        self.no_candidates = False

    def service_response(self, url, **kwargs):
        response = Response()
        response.url = url
        response.status_code = 200
        response.encoding = "utf-8"
        parsed = urlsplit(url)
        if parsed.netloc == "rest.uniprot.org":
            query = parse_qs(parsed.query)
            protein = query["query"][0]
            if query["format"] == ["fasta"]:
                content = f">sp|{protein}|FIXTURE\n{self.sequences[protein]}\n"
            else:
                content = f"Entry\tProtein names\n{protein}\tFixture protein\n"
            response._content = content.encode()
        elif "/api/complex/" in url:
            self.assertEqual(parsed.path.rsplit("/", 1)[-1], self.protein_a)
            if self.no_candidates:
                response.status_code = 404
                response._content = b"{}"
            else:
                response._content = json.dumps([{
                    "modelEntityId": MODEL_ID, "latestVersion": 1,
                    "isComplex": True,
                    "uniprotAccession": [self.protein_a, self.protein_b],
                }]).encode()
        elif url == f"https://alphafold.ebi.ac.uk/api/prediction/{MODEL_ID}":
            response._content = json.dumps(self.metadata).encode()
        elif url == self.metadata[0]["cifUrl"]:
            response._content = (FIXTURES / f"{MODEL_ID}.cif").read_bytes()
        elif parsed.netloc == "swissmodel.expasy.org":
            response._content = b'{"resultset":[]}'
        else:
            self.fail(f"Unexpected network request: {url}")
        return response

    def run_topolink(self, command):
        lines = Path(shlex.split(command)[1]).read_text().splitlines()
        structure_path = next(
            line.split(maxsplit=1)[1]
            for line in lines if line.startswith("structure ")
        )
        model = next(PDBParser(QUIET=True).get_structure(
            "export", structure_path,
        ).get_models())
        observed = [
            line.split() for line in lines if line.strip().startswith("observed ")
        ]
        self.assertEqual(len(observed), 2)
        results = []
        for link in observed:
            atom_a = model[link[2]][int(link[3])]["CB"]
            atom_b = model[link[5]][int(link[6])]["CB"]
            distance = float(np.linalg.norm(atom_a.coord - atom_b.coord))
            results.append(
                f"  LINK: {' '.join(link[1:4])} CB "
                f"{' '.join(link[4:7])} CB "
                f"{distance:.3f} {distance + 1:.3f} "
                "YES 0.000 100.000 OK: FOUND\n"
            )
        return StringIO("".join(results))

    def run_pipeline(self, output, options):
        with patch(
            "requests.get", side_effect=self.service_response,
        ), patch(
            "claudio.module02.src_distance_reevaluation.algorithm."
            "calc_site_distances.os.popen", side_effect=self.run_topolink,
        ) as topolink, patch.object(
            pipeline, "db_updater",
            side_effect=AssertionError("AF mode must not update structure databases"),
        ) as database:
            result = CliRunner().invoke(pipeline.main, [
                "-i", str(self.input), "-o", str(output),
                "-bldb", str(self.root / "missing-db"),
                "-x", "K", "-pc", "0", "-rt", "False", "-s", "False",
                "-v", "0", *options,
            ])
        self.assertEqual(result.exit_code, 0,
                         result.output + repr(result.exception))
        database.assert_not_called()
        return topolink

    def check_source_outputs(self, output):
        stem = "crosslinks.run1"
        m1_path = output / f"{stem}.sqcs"
        ops_path = output / f"{stem}.sqcs_ops.csv"
        structure_path = output / f"{stem}.sqcs_structdi.csv"
        for path in (m1_path, ops_path, structure_path):
            self.assertTrue(path.is_file(), str(path))
        structural = pd.read_csv(structure_path, index_col=0)
        final = pd.read_csv(output / f"{stem}_final.csv", index_col=0)
        extended = pd.read_csv(output / f"{stem}_extended.csv", index_col=0)
        pd.testing.assert_frame_equal(final[structural.columns], structural)
        pd.testing.assert_frame_equal(extended, final)
        self.assertEqual(final["source_record_id"].tolist(), [
            "record_000001", "record_000002", "record_000003",
        ])
        self.assertEqual(final["experiment"].tolist(), [
            "condition_A", "condition_B", "condition_C",
        ])
        self.assertEqual(final["position_status_a"].tolist(),
                         ["corrected", "valid", "ambiguous"])
        self.assertEqual(final["original_pos_a"].tolist(), [76, 23, 100])
        self.assertEqual(final["pos_a"].tolist(), [67, 23, 100])
        self.assertTrue(pd.isna(final.iloc[2]["candidate_id"]))
        duplicates = pd.read_csv(output / f"{stem}_duplicate_sources.csv")
        self.assertEqual(duplicates.to_dict("records"), [{
            "retained_source_record_id": "record_000001",
            "duplicate_source_record_id": "record_000004",
        }])
        return final

    def test_manual_and_auto_modes_keep_provenance_through_all_modules(self):
        fixture_model = next(MMCIFParser(QUIET=True).get_structure(
            "fixture", str(FIXTURES / f"{MODEL_ID}.cif"),
        ).get_models())
        expected = [
            float(np.linalg.norm(fixture_model["A"][67]["CB"].coord
                                 - fixture_model["B"][11]["CB"].coord)),
            float(np.linalg.norm(fixture_model["B"][23]["CB"].coord
                                 - fixture_model["A"][67]["CB"].coord)),
        ]
        for mode, options in (
            ("manual", ["--af-complex-model", MODEL_ID]),
            ("auto", ["--af-complex-auto"]),
        ):
            with self.subTest(mode=mode):
                output = self.root / mode
                topolink = self.run_pipeline(output, options)
                topolink.assert_called_once()
                final = self.check_source_outputs(output)
                self.assertEqual(final["chain_a"].tolist(), ["A", "B", "A"])
                self.assertEqual(final["chain_b"].tolist(), ["B", "A", "B"])
                self.assertTrue(final["af_model_id"].eq(MODEL_ID).all())
                self.assertEqual(final["distance_status"].tolist(),
                                 ["computed", "computed", "not_evaluable"])
                np.testing.assert_allclose(final.iloc[:2]["eucl_dist"], expected, atol=0.001)
                np.testing.assert_allclose(final.iloc[:2]["topo_dist"],
                                           np.array(expected) + 1, atol=0.001)
                unresolved = final.iloc[2]
                self.assertEqual(unresolved["distance_reason"], "unresolved_position")
                self.assertTrue(pd.isna(unresolved["topo_dist"]))
                self.assertTrue(pd.isna(unresolved["homo_pep_overl"]))
                for column in ("path", "af_metadata_path"):
                    self.assertTrue(Path(final.iloc[0][column]).is_file())
                plt.close("all")

    def test_native_lower_bound_is_preserved_to_final_output(self):
        # These LINK values reproduce the native TopoLink example.
        native_output = (
            "  LINK: LYS A 67 CB LYS B 11 CB 62.907 66.751 "
            "YES 0.000 35.000 BAD: LONG\n"
            "  LINK: LYS A 67 CB LYS B 23 CB 75.957 >75.957 "
            "YES 0.000 35.000 BAD: EUCL\n"
        )
        output = self.root / "lower_bound"
        with patch.object(
            self, "run_topolink", return_value=StringIO(native_output),
        ):
            topolink = self.run_pipeline(output, ["--af-complex-auto"])
        topolink.assert_called_once()
        final = self.check_source_outputs(output)
        first, bounded, unresolved = (final.iloc[i] for i in range(3))
        self.assertEqual(first["distance_status"], "computed")
        self.assertAlmostEqual(first["topo_dist"], 66.751)
        self.assertEqual(bounded["distance_status"], "lower_bound")
        self.assertEqual(bounded["distance_reason"], "topolink_lower_bound")
        self.assertEqual(bounded["topo_dist_raw"], ">75.957")
        self.assertAlmostEqual(bounded["topo_dist_lower_bound"], 75.957)
        self.assertTrue(pd.isna(bounded["topo_dist"]))
        self.assertFalse(bounded["XL_confirmed"])
        self.assertEqual(unresolved["distance_reason"], "unresolved_position")
        self.assertTrue(pd.isna(unresolved["topo_dist_lower_bound"]))

    def test_auto_mode_preserves_sources_when_afdb_has_no_candidate(self):
        self.no_candidates = True
        output = self.root / "no_candidates"
        topolink = self.run_pipeline(output, ["--af-complex-auto"])
        topolink.assert_not_called()
        final = self.check_source_outputs(output)
        self.assertTrue(final["af_match_status"].eq("no_candidate_model").all())
        self.assertTrue(final["distance_status"].eq("not_evaluable").all())
        self.assertTrue(final["distance_reason"].eq("no_candidate_model").all())
        self.assertTrue(final[["eucl_dist", "topo_dist"]].isna().all().all())
        self.assertFalse((output / "homomers").exists())
        self.assertFalse((output / "heteromers").exists())


if __name__ == "__main__":
    unittest.main()
