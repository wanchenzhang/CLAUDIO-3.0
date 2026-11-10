"""Offline tests for publication metadata in standardization outputs."""

from copy import deepcopy
import csv
import hashlib
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from click.testing import CliRunner

from claudio.dataset_catalog import cli, publication_metadata, xifdr
from claudio.module01.src.io.read_in import read_inputfile
import test_dataset_catalog_publication_metadata as context_tests

CORE = [
    "peptide1", "peptide2", "position1", "position2",
    "k_pos1", "k_pos2", "entry1", "entry2",
]


class TestPublicationStandardize(unittest.TestCase):
    def setUp(self):
        context = context_tests.TestDatasetPublicationMetadata(
            "test_loads_selected_paper_from_discovery"
        )
        self.addCleanup(context.doCleanups)
        context.setUp()
        self.context = context
        self.root = context.root
        self.metadata_path = context.candidates
        self.bundle = context.bundle()
        self.bundle.update(
            format="xifdr_psm",
            proteins={},
            source_columns=["source_peptide", "source_fdr"],
        )
        self.bundle["records"][0].update(
            record_reason="",
            source_fields={
                "source_peptide": "AKK",
                "source_fdr": "0.005",
            },
            fdr={"psm": "0.005"},
        )
        self.prepared = {
            "rows": [{
                "source_record_id": "PXD000001:source:000001",
                "pep_a": "AKK",
                "pep_b": "AKK",
                "pos_a": 2,
                "pos_b": 3,
                "res_pos_a": 2,
                "res_pos_b": 3,
                "unip_id_a": "P00001",
                "unip_id_b": "P00002",
            }],
            "audit": {
                "records": deepcopy(self.bundle["records"]),
                "assessment_performed": False,
            },
            "fdr_policy": {
                "level": "psm",
                "cutoff": 0.01,
                "kind": "source_reported_numeric_selection",
                "target_decoy_validation_performed": False,
            },
        }
        self.results_path = self.root / "results.csv"
        self.results_path.write_text(
            "offline source\n", encoding="utf-8"
        )
        self.fasta_path = self.root / "search.fasta"
        self.fasta_path.write_text(
            ">P00001\nAKK\n", encoding="utf-8"
        )
        self.mapping_path = self.root / "mapping.json"
        self.mapping_path.write_text("{}", encoding="utf-8")

    def bind(self):
        return publication_metadata.bind_dataset_metadata(
            self.bundle,
            self.metadata_path,
            dataset_accession="PXD000001",
        )

    def export(
        self, bundle=None, prepared=None,
        name="output", pipeline=True,
    ):
        return Path(xifdr.write_xifdr_outputs(
            self.bind() if bundle is None else bundle,
            {},
            self.prepared if prepared is None else prepared,
            self.root / name,
            pipeline_policy=(
                {"input_scope": "source_numeric_selection"}
                if pipeline else None
            ),
        ))

    def read_csv(self, path):
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            return reader.fieldnames, list(reader)

    def read_json(self, path):
        return json.loads(path.read_bytes())

    def invoke(
        self, *, metadata=True,
        accession="PXD000001", audit=False,
    ):
        destination = self.root / "cli_output"
        arguments = [
            "standardize",
            "--results", str(self.results_path),
            "--fasta", str(self.fasta_path),
            "--dataset-accession", accession,
            "--output", str(destination),
        ]
        if metadata:
            arguments += [
                "--dataset-metadata", str(self.metadata_path)
            ]
        if audit:
            arguments += ["--audit-only"]
        else:
            arguments += [
                "--mapping", str(self.mapping_path),
                "--fdr-level", "psm",
                "--fdr-cutoff", "0.01",
                "--crosslinker", "BS3",
                "--residue", "K",
            ]

        def write_input(bundle, mappings, output_dir, **options):
            return xifdr.write_xifdr_outputs(
                bundle,
                mappings,
                deepcopy(self.prepared),
                output_dir,
                pipeline_policy={
                    "input_scope": "source_numeric_selection"
                },
            )

        with patch.object(
            cli, "_detect_psm_format",
            return_value="xifdr_psm",
        ), patch.object(
            cli, "read_xifdr_psms",
            return_value=deepcopy(self.bundle),
            create=True,
        ), patch.object(
            xifdr, "read_xifdr_psms",
            return_value=deepcopy(self.bundle),
        ), patch.object(
            cli, "_read_protein_mapping", return_value={},
        ), patch.object(
            cli, "write_claudio_input",
            side_effect=write_input, create=True,
        ), patch.object(
            xifdr, "write_claudio_input",
            side_effect=write_input,
        ):
            result = CliRunner().invoke(cli.main, arguments)
        return result, destination

    def test_writer_adds_publication_fields(self):
        destination = self.export()
        columns, rows = self.read_csv(
            destination / "claudio_input.csv"
        )
        self.assertIn(
            "xlms_dataset_accession",
            columns,
            "Dataset publication fields are not exported",
        )
        self.assertEqual(columns[:8], CORE)
        self.assertEqual(
            rows[0]["xlms_dataset_accession"], "PXD000001"
        )
        self.assertEqual(
            rows[0]["xlms_publication_status"], "published"
        )
        self.assertEqual(
            json.loads(rows[0]["xlms_paper_dois_json"]),
            ["10.1234/example"],
        )
        self.assertEqual(
            json.loads(rows[0]["xlms_papers_json"])[0]["title"],
            "An example article",
        )

    def test_output_metadata_reloads_without_original_discovery(self):
        destination = self.export()
        self.context.discovery.rename(
            self.root / "discovery_moved"
        )
        packet = publication_metadata.load_dataset_metadata(
            destination / "dataset_metadata.json",
            dataset_accession="PXD000001",
        )
        self.assertEqual(
            packet["record"]["papers"][0]["doi"],
            "10.1234/example",
        )

    def test_writer_preserves_audit_fdr_and_input_objects(self):
        original_bundle = deepcopy(self.bundle)
        original_prepared = deepcopy(self.prepared)
        destination = self.export()
        manifest = self.read_json(
            destination / "manifest.json"
        )
        self.assertEqual(
            manifest["fdr_policy"],
            self.prepared["fdr_policy"],
        )
        self.assertEqual(
            self.read_json(destination / "psm_audit.json"),
            self.prepared["audit"],
        )
        self.assertEqual(self.bundle, original_bundle)
        self.assertEqual(self.prepared, original_prepared)

    def test_prepared_input_also_carries_metadata(self):
        destination = self.export(pipeline=False)
        columns, rows = self.read_csv(
            destination / "af_input.csv"
        )
        self.assertIn("xlms_papers_json", columns)
        self.assertEqual(
            rows[0]["source_record_id"],
            "PXD000001:source:000001",
        )
        self.assertEqual(
            rows[0]["xlms_dataset_accession"], "PXD000001"
        )

    def test_empty_input_keeps_metadata_headers_and_snapshot(self):
        prepared = deepcopy(self.prepared)
        prepared["rows"] = []
        destination = self.export(prepared=prepared)
        columns, rows = self.read_csv(
            destination / "claudio_input.csv"
        )
        self.assertEqual(columns[:8], CORE)
        self.assertIn("xlms_papers_json", columns)
        self.assertEqual(rows, [])
        self.assertTrue(
            (destination / "dataset_metadata.json").is_file()
        )

    def test_rejects_source_snapshot_collision_before_creating_output(self):
        prepared = deepcopy(self.prepared)
        prepared["extra_source_bodies"] = {
            "dataset_metadata.json": b"unrelated"
        }
        with self.assertRaises(ValueError):
            self.export(prepared=prepared)
        self.assertFalse((self.root / "output").exists())

    def test_rejects_reserved_paper_column_before_creating_output(self):
        prepared = deepcopy(self.prepared)
        prepared["rows"][0]["xlms_papers_json"] = "unrelated"
        with self.assertRaises(ValueError):
            self.export(prepared=prepared)
        self.assertFalse((self.root / "output").exists())

    def test_changed_bound_metadata_creates_no_output(self):
        bound = self.bind()
        self.metadata_path.write_bytes(
            self.metadata_path.read_bytes() + b" "
        )
        with self.assertRaises(ValueError):
            self.export(bundle=bound)
        self.assertFalse((self.root / "output").exists())

    def test_legacy_export_has_no_added_metadata(self):
        destination = self.export(bundle=self.bundle)
        columns, rows = self.read_csv(
            destination / "claudio_input.csv"
        )
        self.assertEqual(
            columns,
            [
                "peptide1", "peptide2",
                "position1", "position2",
                "k_pos1", "k_pos2",
                "entry1", "entry2",
                "xlms_source_record_id",
            ],
        )
        self.assertEqual(
            rows[0]["xlms_source_record_id"],
            "PXD000001:source:000001",
        )
        self.assertNotIn(
            "dataset_metadata",
            self.read_json(destination / "manifest.json"),
        )
        self.assertFalse(
            (destination / "dataset_metadata.json").exists()
        )

    def test_manifest_pins_every_packaged_file(self):
        destination = self.export()
        manifest = self.read_json(
            destination / "manifest.json"
        )
        self.assertEqual(
            manifest["dataset_metadata"]["dataset_accession"],
            "PXD000001",
        )
        self.assertIn(
            "dataset_metadata.json",
            manifest["output_sha256"],
        )
        for name, expected in manifest["output_sha256"].items():
            self.assertEqual(
                hashlib.sha256(
                    (destination / name).read_bytes()
                ).hexdigest(),
                expected,
            )

    def test_m1_reader_preserves_metadata_from_the_real_writer(self):
        destination = self.export()
        renamed = [
            "pep_a", "pep_b", "pos_a", "pos_b",
            "res_pos_a", "res_pos_b",
            "unip_id_a", "unip_id_b",
        ]
        data = read_inputfile(
            str(destination / "claudio_input.csv"),
            dict(zip(CORE, renamed)),
        )
        self.assertEqual(
            data.iloc[0]["xlms_dataset_accession"],
            "PXD000001",
        )
        self.assertEqual(
            json.loads(
                data.iloc[0]["xlms_paper_dois_json"]
            ),
            ["10.1234/example"],
        )
        self.assertEqual(
            data.iloc[0]["xlms_source_record_id"],
            "PXD000001:source:000001",
        )

    def test_cli_help_exposes_dataset_metadata_option(self):
        result = CliRunner().invoke(
            cli.main, ["standardize", "--help"]
        )
        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn(
            "--dataset-metadata",
            result.output,
            "Dataset metadata CLI option is not implemented",
        )

    def test_cli_passes_metadata_to_the_common_writer(self):
        result, destination = self.invoke()
        self.assertEqual(
            result.exit_code, 0,
            f"{result.output}\n{result.exception}",
        )
        _, rows = self.read_csv(
            destination / "claudio_input.csv"
        )
        self.assertEqual(
            rows[0]["xlms_dataset_accession"], "PXD000001"
        )
        self.assertTrue(
            (destination / "dataset_metadata.json").is_file()
        )
        api_destination = self.export(name="api_output")
        for filename in (
            "claudio_input.csv", "psm_audit.json"
        ):
            self.assertEqual(
                (destination / filename).read_bytes(),
                (api_destination / filename).read_bytes(),
            )

    def test_cli_without_metadata_retains_legacy_output(self):
        result, destination = self.invoke(metadata=False)
        self.assertEqual(
            result.exit_code, 0,
            f"{result.output}\n{result.exception}",
        )
        columns, _ = self.read_csv(
            destination / "claudio_input.csv"
        )
        self.assertNotIn("xlms_papers_json", columns)
        self.assertFalse(
            (destination / "dataset_metadata.json").exists()
        )

    def test_cli_rejects_wrong_dataset_without_creating_output(self):
        result, destination = self.invoke(
            accession="PXD999999"
        )
        self.assertNotEqual(result.exit_code, 0)
        self.assertFalse(destination.exists())

    def test_cli_rejects_changed_evidence_without_creating_output(self):
        self.metadata_path.write_bytes(
            self.metadata_path.read_bytes() + b" "
        )
        result, destination = self.invoke()
        self.assertNotEqual(result.exit_code, 0)
        self.assertFalse(destination.exists())

    def test_cli_rejects_metadata_in_audit_only_mode(self):
        result, destination = self.invoke(audit=True)
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("--dataset-metadata", result.output)
        self.assertFalse(destination.exists())


if __name__ == "__main__":
    unittest.main()
