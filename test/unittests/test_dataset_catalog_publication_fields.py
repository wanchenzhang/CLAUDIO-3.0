
"""Synthetic I/O tests for paper metadata in M1 and M4 outputs."""

import importlib
import json
import os
from pathlib import Path
import tempfile
import unittest

import pandas as pd
from pandas.testing import assert_frame_equal

FIELDS = (
    "xlms_dataset_accession", "xlms_publication_status",
    "xlms_paper_dois_json", "xlms_paper_pmids_json",
    "xlms_papers_json",
)
LEGACY = [
    "pdb_id", "chain_a", "chain_b", "pdb_pos_a", "pdb_pos_b",
    "pos_a", "pos_b", "swiss_model_homology", "XL_type",
]


def row(
    dataset="PXD000001",
    doi="10.1234/first",
    left="P11111",
    right="P22222",
):
    return {
        "unip_id_a": left, "unip_id_b": right,
        "pos_a": 2, "pos_b": 3,
        "pep_a": "AKAAK", "pep_b": "AAKAA",
        "res_pos_a": 2.0, "res_pos_b": 3.0,
        "seq_a": "AKAAK", "seq_b": "AAKAA",
        "pep_a_repeated": False, "pep_b_repeated": False,
        "pep_a_repeat_count": 1, "pep_b_repeat_count": 1,
        "pdb_id": "fixture", "pdb_method": "fixture",
        "pdb_resolution": "fixture",
        "chain_a": "A", "chain_b": "B",
        "pdb_pos_a": 2, "pdb_pos_b": 3,
        "pLDDT_a": 90.0, "pLDDT_b": 90.0,
        "topo_dist": 10.0, "eucl_dist": 8.0,
        "homo_adjacency": 0.0, "homo_int_overl": 0.0,
        "homo_pep_overl": 0.0,
        "evidence": "synthetic_fixture",
        "XL_type": "inter", "XL_confirmed": False,
        "swiss_model_homology": "heteromer",
        "dataset_accession": dataset,
        "xlms_source_record_id": dataset + ":fixture:1",
        "source_record_id": dataset + ":internal:1",
        "xlms_dataset_accession": dataset,
        "xlms_publication_status": "published",
        "xlms_paper_dois_json": json.dumps([doi]),
        "xlms_paper_pmids_json": "[]",
        "xlms_papers_json": json.dumps([{
            "doi": doi,
            "pmid": None,
            "title": "Synthetic article",
            "venue": "Fixture journal",
            "year": 2025,
            "authors": [],
        }]),
    }


class TestPublicationOutputFields(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def api(self):
        try:
            return importlib.import_module(
                "claudio.dataset_catalog.publication_fields"
            )
        except ModuleNotFoundError as error:
            if error.name != (
                "claudio.dataset_catalog.publication_fields"
            ):
                raise
            self.fail(
                "Publication field propagation helper "
                "is not implemented"
            )

    def export(self, rows, score=False):
        self.api()
        data = pd.DataFrame(rows)
        original = data.copy(deep=True)
        output = self.root / "output"
        output.mkdir()
        (output / "homomers").mkdir()
        (output / "heteromers").mkdir()

        writer = importlib.import_module(
            "claudio.module04.src.io.write_outs"
        )
        writer.write_outputs(
            data, "fixture", score, str(output) + os.sep
        )
        interactions = [
            path for path in output.rglob("*.csv")
            if path.name != "fixture_final.csv"
        ]
        return data, original, output, interactions

    def assert_metadata_rows(self, actual, expected):
        for field in FIELDS:
            self.assertIn(field, actual.columns)
        for field in (
            "dataset_accession",
            "xlms_source_record_id",
            "source_record_id",
        ):
            self.assertIn(field, actual.columns)

        self.assertEqual(len(actual), len(expected))
        originals = {
            value["xlms_source_record_id"]: value
            for value in expected
        }
        for _, value in actual.iterrows():
            source = originals[
                value["xlms_source_record_id"]
            ]
            for field in FIELDS + (
                "dataset_accession", "source_record_id"
            ):
                self.assertEqual(
                    value[field], source[field]
                )

    def test_interaction_csv_keeps_dataset_specific_paper_rows(self):
        rows = [
            row(),
            row(
                "PXD000002", "10.1234/second",
                "P22222", "P11111",
            ),
        ]
        _, _, _, files = self.export(rows)
        self.assertEqual(len(files), 1)
        actual = pd.read_csv(
            files[0], keep_default_na=False
        )
        self.assert_metadata_rows(actual, rows)

    def test_final_csv_keeps_full_paper_metadata(self):
        rows = [
            row(), row("PXD000002", "10.1234/second")
        ]
        _, _, output, _ = self.export(rows)
        self.assert_metadata_rows(
            pd.read_csv(
                output / "fixture_final.csv",
                keep_default_na=False,
            ),
            rows,
        )

    def test_unrelated_interaction_groups_keep_their_own_papers(self):
        rows = [
            row(),
            row(
                "PXD000002", "10.1234/second",
                "P33333", "P44444",
            ),
        ]
        _, _, _, files = self.export(rows)
        self.assertEqual(len(files), 2)
        observed = set()

        for path in files:
            actual = pd.read_csv(
                path, keep_default_na=False
            )
            self.assertEqual(len(actual), 1)
            source_id = actual.iloc[0][
                "xlms_source_record_id"
            ]
            source = next(
                value for value in rows
                if value["xlms_source_record_id"]
                == source_id
            )
            self.assert_metadata_rows(actual, [source])
            observed.add(source_id)

        self.assertEqual(
            observed,
            {
                value["xlms_source_record_id"]
                for value in rows
            },
        )

    def test_homomer_csv_keeps_paper_fields(self):
        value = row(
            left="P11111", right="P11111"
        )
        value["swiss_model_homology"] = ""
        _, _, _, files = self.export([value])
        self.assertEqual(len(files), 1)
        self.assert_metadata_rows(
            pd.read_csv(
                files[0], keep_default_na=False
            ),
            [value],
        )

    def test_legacy_interaction_columns_are_unchanged(self):
        value = row()
        for field in FIELDS:
            value.pop(field)
        _, _, _, files = self.export([value])
        actual = pd.read_csv(
            files[0], keep_default_na=False
        )
        self.assertEqual(
            list(actual.columns), LEGACY
        )

    def test_partial_metadata_is_preserved_without_inventing_status(self):
        value = row()
        for field in FIELDS:
            if field != "xlms_paper_dois_json":
                value.pop(field)

        _, _, _, files = self.export([value])
        actual = pd.read_csv(
            files[0], keep_default_na=False
        )
        self.assertEqual(
            actual.iloc[0]["xlms_paper_dois_json"],
            '["10.1234/first"]',
        )
        self.assertNotIn(
            "xlms_publication_status", actual.columns
        )
        self.assertNotIn(
            "xlms_papers_json", actual.columns
        )

    def test_m1_reader_preserves_paper_fields_and_source_ids(self):
        self.api()
        metadata = row()
        values = {
            "peptide1": "AKAAK",
            "peptide2": "AAKAA",
            "position1": 2,
            "position2": 3,
            "k_pos1": 2,
            "k_pos2": 3,
            "entry1": "P11111",
            "entry2": "P22222",
        }
        for field in FIELDS + (
            "dataset_accession", "xlms_source_record_id"
        ):
            values[field] = metadata[field]

        path = self.root / "input.csv"
        pd.DataFrame([values]).to_csv(
            path, index=False
        )
        reader = importlib.import_module(
            "claudio.module01.src.io.read_in"
        )
        result = reader.read_inputfile(
            str(path),
            {
                "peptide1": "pep_a",
                "peptide2": "pep_b",
                "position1": "pos_a",
                "position2": "pos_b",
                "k_pos1": "res_pos_a",
                "k_pos2": "res_pos_b",
                "entry1": "unip_id_a",
                "entry2": "unip_id_b",
            },
        )
        for field in FIELDS + (
            "dataset_accession", "xlms_source_record_id"
        ):
            self.assertEqual(
                result.iloc[0][field],
                metadata[field],
            )
        self.assertEqual(
            result.iloc[0]["source_record_id"],
            "record_000001",
        )

    def test_export_does_not_mutate_input_rows(self):
        data, original, _, _ = self.export([row()])
        assert_frame_equal(data, original)

    def test_scored_output_keeps_paper_fields(self):
        rows = [row()]
        _, _, output, files = self.export(
            rows, score=True
        )
        self.assert_metadata_rows(
            pd.read_csv(
                files[0], keep_default_na=False
            ),
            rows,
        )
        final = pd.read_csv(
            output / "fixture_final.csv",
            keep_default_na=False,
        )
        self.assert_metadata_rows(final, rows)
        self.assertIn(
            "homo_adjacency", final.columns
        )


if __name__ == "__main__":
    unittest.main()
