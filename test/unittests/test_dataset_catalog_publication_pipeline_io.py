"""Exercise publication fields through real pipeline I/O."""

from copy import deepcopy
import json
from pathlib import Path
import unittest

import pandas as pd

from claudio.dataset_catalog import publication_metadata, xifdr
from claudio.module01.src.io.read_in import read_inputfile
from claudio.module01.src.io.write_out import write_outputs as write_m1
from claudio.module02.src_structure_search.io.read_input import read_in as read_structure
from claudio.module02.src_structure_search.io.write_out import write_output as write_structure
from claudio.module02.src_distance_reevaluation.io.read_uniprot_search_out import read_unipsearch_out
from claudio.module02.src_distance_reevaluation.io.write_out import write_output as write_distance
from claudio.module03.src.io.read_in import read_in as read_ops
from claudio.module03.src.io.write_out import write_output as write_ops
from claudio.module04.src.io.read_ins import read_inputs as read_m4
from claudio.module04.src.io.write_outs import write_outputs as write_m4

import test_dataset_catalog_publication_standardize as standardize_tests

PUBLICATION_COLUMNS = (
    "xlms_dataset_accession",
    "xlms_publication_status",
    "xlms_paper_dois_json",
    "xlms_paper_pmids_json",
    "xlms_papers_json",
)
CORE = (
    "peptide1", "peptide2", "position1", "position2",
    "k_pos1", "k_pos2", "entry1", "entry2",
)
INTERNAL = (
    "pep_a", "pep_b", "pos_a", "pos_b",
    "res_pos_a", "res_pos_b", "unip_id_a", "unip_id_b",
)


class TestPublicationPipelineIO(unittest.TestCase):
    def setUp(self):
        self.fixture = standardize_tests.TestPublicationStandardize(
            "test_writer_adds_publication_fields"
        )
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.setUp()
        self.root = self.fixture.root

    def make_input(self, *, mixed=False, legacy=False):
        destination = self.fixture.export(
            bundle=self.fixture.bundle if legacy else None,
            name="first_standardized",
        )
        path = destination / "claudio_input.csv"
        first = pd.read_csv(path)
        if not mixed:
            return path, first

        bundle = deepcopy(self.fixture.bundle)
        prepared = deepcopy(self.fixture.prepared)
        identifier = "PXD000002:source:000001"
        bundle["provenance"]["dataset_accession"] = "PXD000002"
        bundle["records"][0]["source_record_id"] = identifier
        prepared["rows"][0].update(
            source_record_id=identifier,
            unip_id_a="P00002",
            unip_id_b="P00001",
        )
        prepared["audit"]["records"] = deepcopy(bundle["records"])
        bound = publication_metadata.bind_dataset_metadata(
            bundle,
            self.fixture.metadata_path,
            dataset_accession="PXD000002",
        )
        second_directory = Path(xifdr.write_xifdr_outputs(
            bound,
            {},
            prepared,
            self.root / "second_standardized",
            pipeline_policy={
                "input_scope": "source_numeric_selection"
            },
        ))
        second = pd.read_csv(
            second_directory / "claudio_input.csv"
        )
        combined = pd.concat(
            [first, second], ignore_index=True
        )
        path = self.root / "mixed_input.csv"
        combined.to_csv(path, index=False)
        return path, combined

    def expected_metadata(self, original):
        return {
            row["xlms_source_record_id"]: {
                column: row[column]
                for column in PUBLICATION_COLUMNS
            }
            for _, row in original.iterrows()
        }

    def check_metadata(self, data, expected, *, legacy=False):
        self.assertIn(
            "xlms_source_record_id", data.columns
        )
        if legacy:
            for column in PUBLICATION_COLUMNS:
                self.assertNotIn(column, data.columns)
            return

        for column in PUBLICATION_COLUMNS:
            self.assertIn(column, data.columns)

        for _, row in data.iterrows():
            identifier = row["xlms_source_record_id"]
            self.assertIn(identifier, expected)
            for column in PUBLICATION_COLUMNS:
                if column.endswith("_json"):
                    self.assertEqual(
                        json.loads(row[column]),
                        json.loads(
                            expected[identifier][column]
                        ),
                    )
                else:
                    self.assertEqual(
                        row[column],
                        expected[identifier][column],
                    )

    def roundtrip(
        self, input_path, original, *,
        expanded=False, legacy=False,
        scoring=False, homomer=False,
    ):
        expected = (
            {} if legacy
            else self.expected_metadata(original)
        )
        data = read_inputfile(
            str(input_path),
            dict(zip(CORE, INTERNAL)),
        )
        self.check_metadata(
            data, expected, legacy=legacy
        )

        # Supply synthetic fields without sequence or structure lookup.
        for column, value in {
            "seq_a": "AKK",
            "seq_b": "AKK",
            "pep_a_repeated": False,
            "pep_b_repeated": False,
            "pep_a_repeat_count": 0,
            "pep_b_repeat_count": 0,
        }.items():
            data[column] = value

        if homomer:
            data["unip_id_b"] = data["unip_id_a"]

        directory = self.root / "pipeline_io"
        directory.mkdir()
        output = directory.as_posix() + "/"
        name = "publication_case"

        write_m1(
            data,
            pd.DataFrame({
                "unip_id": ["P00001", "P00002"]
            }),
            name,
            output,
            0,
        )
        sqcs = directory / (name + ".sqcs")

        structure, filename = read_structure(str(sqcs))
        self.check_metadata(
            structure, expected, legacy=legacy
        )
        structure["chain_a"] = "A"
        structure["chain_b"] = "B"
        structure["pdb_id"] = "TEST"
        write_structure(structure, filename, output)
        structure_path = (
            directory / (filename + "_structdi.csv")
        )

        distance = read_unipsearch_out(
            str(structure_path)
        )
        self.check_metadata(
            distance, expected, legacy=legacy
        )
        for column, value in {
            "pdb_pos_a": 2,
            "pdb_pos_b": 3,
            "pLDDT_a": 90.0,
            "pLDDT_b": 91.0,
            "is_interfaced": True,
            "topo_dist": 10.0,
            "eucl_dist": 9.0,
        }.items():
            distance[column] = value

        if expanded:
            clone = distance.iloc[[0]].copy()
            clone.index = ["0_1"]
            distance = pd.concat([distance, clone])

        distance_path = (
            directory / (name + "_distance.csv")
        )
        write_distance(
            distance, str(distance_path)
        )

        ops = read_ops(str(sqcs))
        self.check_metadata(
            ops, expected, legacy=legacy
        )
        ops["homo_adjacency"] = 0.5
        ops["homo_int_overl"] = 0.0
        ops["homo_pep_overl"] = False
        write_ops(ops, name, output)

        merged = read_m4(
            str(distance_path),
            str(directory / (name + "_ops.csv")),
        )
        self.check_metadata(
            merged, expected, legacy=legacy
        )
        self.assertEqual(
            len(merged),
            len(original) + int(expanded),
        )

        merged["evidence"] = "Synthetic I/O regression fixture"
        merged["XL_type"] = "inter"
        merged["XL_confirmed"] = False
        merged["swiss_model_homology"] = ""

        final_directory = self.root / "final"
        final_directory.mkdir()

        # Isolate metadata checks from output-directory creation.
        (final_directory / "homomers").mkdir()
        (final_directory / "heteromers").mkdir()

        write_m4(
            merged,
            name,
            scoring,
            final_directory.as_posix() + "/",
        )
        final = pd.read_csv(
            final_directory / (name + "_final.csv"),
            index_col=0,
        )
        self.check_metadata(
            final, expected, legacy=legacy
        )
        self.assertEqual(len(final), len(merged))

        group_root = final_directory / (
            "homomers" if homomer else "heteromers"
        )
        groups = list(group_root.rglob("*.csv"))
        self.assertEqual(len(groups), 1)
        interaction = pd.read_csv(groups[0])

        if legacy:
            for column in PUBLICATION_COLUMNS:
                self.assertNotIn(
                    column, interaction.columns
                )
        else:
            self.check_metadata(
                interaction, expected
            )
            self.assertEqual(
                len(interaction), len(merged)
            )
        return final, interaction

    def test_published_context_survives_all_pipeline_io_stages(self):
        path, original = self.make_input()
        final, _ = self.roundtrip(path, original)
        self.assertEqual(
            final.iloc[0]["xlms_dataset_accession"],
            "PXD000001",
        )

    def test_mixed_datasets_and_expanded_rows_keep_their_own_papers(self):
        path, original = self.make_input(mixed=True)
        final, interaction = self.roundtrip(
            path, original, expanded=True
        )
        for data in (final, interaction):
            counts = (
                data["xlms_dataset_accession"]
                .value_counts().to_dict()
            )
            self.assertEqual(
                counts,
                {"PXD000001": 2, "PXD000002": 1},
            )
            for _, row in data.iterrows():
                expected = (
                    ["10.1234/example"]
                    if row["xlms_dataset_accession"] == "PXD000001"
                    else ["10.1234/other"]
                )
                self.assertEqual(
                    json.loads(
                        row["xlms_paper_dois_json"]
                    ),
                    expected,
                )

    def test_legacy_pipeline_io_adds_no_publication_fields(self):
        path, original = self.make_input(legacy=True)
        self.roundtrip(
            path, original, legacy=True
        )

    def test_scored_output_keeps_publication_fields(self):
        path, original = self.make_input()
        final, _ = self.roundtrip(
            path, original, scoring=True
        )
        self.assertIn(
            "homo_adjacency", final.columns
        )
        self.assertIn(
            "homo_int_overl", final.columns
        )

    def test_homomer_interaction_output_keeps_publication_fields(self):
        path, original = self.make_input()
        final, interaction = self.roundtrip(
            path, original, homomer=True
        )
        self.assertEqual(
            final.iloc[0]["xlms_dataset_accession"],
            "PXD000001",
        )
        self.assertEqual(
            json.loads(
                interaction.iloc[0]["xlms_paper_dois_json"]
            ),
            ["10.1234/example"],
        )


if __name__ == "__main__":
    unittest.main()
