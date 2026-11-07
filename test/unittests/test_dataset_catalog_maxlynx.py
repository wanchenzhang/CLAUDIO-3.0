import copy
import csv
import importlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from collections import Counter
from pathlib import Path

from click.testing import CliRunner

from claudio.dataset_catalog import cli
from claudio.dataset_catalog.cli import _read_reference_projection


class TestMaxLynxStandardize(unittest.TestCase):
    def setUp(self):
        self.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures" / "maxlynx_pxd044630"
        )
        self.expected = json.loads(
            (self.fixture / "expected.json").read_text(encoding="utf-8")
        )
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def api(self):
        self.assertIsNotNone(
            importlib.util.find_spec("claudio.dataset_catalog.maxlynx"),
            "Missing MaxLynx adapter",
        )
        return importlib.import_module("claudio.dataset_catalog.maxlynx")

    def read_bundle(self, results=None, parameters=None, backend="pyxlms"):
        return self.api().read_maxlynx_csms(
            self.fixture / "csms.txt" if results is None else results,
            self.fixture / "search.fasta",
            dataset_accession="PXD044630",
            parameters_path=(
                self.fixture / "parameters.xml"
                if parameters is None else parameters
            ),
            backend=backend,
        )

    def mapping(self):
        return _read_reference_projection(
            self.fixture / "evidence" / "reference_mapping.json"
        )

    def export(self, bundle=None, mapping=None, name="input", **options):
        arguments = {
            "allowed_crosslinkers": {"DSSO"},
            "allowed_residues": {"K"},
            "export_unverified_candidates": True,
        }
        arguments.update(options)
        return Path(self.api().write_maxlynx_claudio_input(
            self.read_bundle() if bundle is None else bundle,
            self.mapping() if mapping is None else mapping,
            self.root / name,
            **arguments,
        ))

    def load_output(self, directory):
        with (directory / "claudio_input.csv").open(
            encoding="utf-8", newline=""
        ) as handle:
            reader = csv.DictReader(handle)
            columns = reader.fieldnames
            rows = list(reader)
        audit = json.loads(
            (directory / "psm_audit.json").read_text(encoding="utf-8")
        )
        manifest = json.loads(
            (directory / "manifest.json").read_text(encoding="utf-8")
        )
        return columns, rows, audit, manifest

    def mutated_source(self, transform):
        with (self.fixture / "csms.txt").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames
            rows = list(reader)
        first_pair = next(
            row for row in rows
            if row["Crosslink Product Type"] in ("IntraProXL", "InterProXL")
        )
        transform(first_pair)
        destination = self.root / "changed.txt"
        with destination.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)
        return destination

    def test_import_preserves_source_records_and_endpoints(self):
        bundle = self.read_bundle()
        self.assertEqual(bundle["format"], "maxlynx_csm")
        self.assertEqual(len(bundle["records"]), 1139)
        self.assertEqual(bundle["parser"]["parsed_records"], 59)
        self.assertEqual(
            Counter(record["record_kind"] for record in bundle["records"]),
            {
                "paired_crosslink_candidate": 59,
                "linear_peptide": 930,
                "dead_end": 143,
                "loop_link": 7,
            },
        )
        with (self.fixture / "csms.txt").open(
            encoding="utf-8-sig", newline=""
        ) as handle:
            source = list(csv.reader(handle, delimiter="\t"))
        self.assertEqual(bundle["source_columns"], source[0])
        for record, values in zip(bundle["records"], source[1:]):
            self.assertEqual(record["source_values"], values)
            self.assertEqual(
                record["source_fields"], dict(zip(source[0], values))
            )
        pairs = [
            record for record in bundle["records"]
            if record["record_kind"] == "paired_crosslink_candidate"
        ]
        self.assertTrue(all(record["record_status"] == "valid" for record in pairs))
        self.assertEqual(
            Counter(record["parser_orientation"] for record in pairs),
            {"same": 40, "swapped": 19},
        )
        first = pairs[0]
        self.assertEqual(first["source_fields"]["Pro_InterLink1"], "202;")
        self.assertEqual(first["endpoints"]["a"]["peptide"], "SPKPSANSVTSPHSK")
        self.assertEqual(first["endpoints"]["a"]["protein_position"], 202)
        self.assertEqual(first["endpoints"]["b"]["protein_position"], 146)

    def test_multiple_localizations_are_retained_without_first_match_selection(self):
        def change(row):
            for name in (
                "Pro_InterLink1", "Pep_InterLink1", "AA_InterLink1",
                "Pro_InterLink2", "Pep_InterLink2", "AA_InterLink2",
                "InterLinks",
            ):
                value = row[name].rstrip(";")
                row[name] = value + ";" + value + ";"
        bundle = self.read_bundle(
            self.mutated_source(change), backend="native"
        )
        record = next(
            item for item in bundle["records"]
            if item["record_kind"] == "paired_crosslink_candidate"
        )
        self.assertEqual(len(record["localization"]["candidates"]), 2)
        self.assertEqual(record["record_status"], "not_evaluable")

    def test_mismatched_localization_lists_remain_unresolved(self):
        def change(row):
            row["Pep_InterLink1"] += "3;"
        bundle = self.read_bundle(
            self.mutated_source(change), backend="native"
        )
        record = next(
            item for item in bundle["records"]
            if item["record_kind"] == "paired_crosslink_candidate"
        )
        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertTrue(record.get("source_error"))
        self.assertEqual(len(bundle["records"]), 1139)

    def test_export_matches_verified_projection_and_preserves_uncertainty(self):
        bundle = self.read_bundle()
        original = copy.deepcopy(bundle)
        directory = self.export(bundle=bundle)
        columns, rows, audit, manifest = self.load_output(directory)
        self.assertEqual(columns[:8], [
            "peptide1", "peptide2", "position1", "position2",
            "k_pos1", "k_pos2", "entry1", "entry2",
        ])
        self.assertEqual(len(rows), 59)
        self.assertEqual(len(audit), 1139)
        for actual, expected in zip(rows, self.expected["rows"]):
            for name, value in expected.items():
                self.assertEqual(actual[name], value)
            self.assertEqual(actual["input_scope"], "unverified_candidate")
            self.assertEqual(actual["fdr_verification_status"], "configured_only")
            for name in (
                "selected_fdr_level", "selected_fdr_value", "selected_fdr_cutoff"
            ):
                self.assertEqual(actual[name], "")
        for source, exported in zip(bundle["records"], audit):
            self.assertEqual(source["source_fields"], exported["source_fields"])
            self.assertEqual(source["source_values"], exported["source_values"])
            self.assertEqual(source["source_record_id"], exported["source_record_id"])
        self.assertEqual(manifest["kind"], "maxlynx_claudio_input")
        self.assertEqual(
            manifest["fdr_policy"]["verification_status"], "configured_only"
        )
        self.assertFalse(manifest["fdr_policy"]["per_record_fdr_applied"])
        self.assertFalse(manifest["assessment_performed"])
        self.assertEqual(bundle, original)

    def test_default_export_requires_explicit_candidate_mode(self):
        with self.assertRaisesRegex(ValueError, "(?i)unverified"):
            self.export(export_unverified_candidates=False)
        self.assertFalse((self.root / "input").exists())

    def test_native_decoy_is_not_promoted_by_target_sequence_matching(self):
        path = self.mutated_source(
            lambda row: row.update(Decoy="decoy")
        )
        bundle = self.read_bundle(path, backend="native")
        first = next(
            record for record in bundle["records"]
            if record["record_kind"] == "paired_crosslink_candidate"
        )
        self.assertNotEqual(first["decoy_status"], "target_target")
        directory = self.export(bundle=bundle)
        _, rows, audit, _ = self.load_output(directory)
        self.assertEqual(len(rows), 58)
        self.assertNotIn(
            first["source_record_id"],
            {row["xlms_source_record_id"] for row in rows},
        )
        excluded = next(
            record for record in audit
            if record["source_record_id"] == first["source_record_id"]
        )
        self.assertIn("decoy_record", excluded["input_reasons"])

    def test_changed_parameters_are_rejected_before_output(self):
        parameters = self.root / "parameters.xml"
        shutil.copyfile(self.fixture / "parameters.xml", parameters)
        bundle = self.read_bundle(parameters=parameters)
        parameters.write_bytes(parameters.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum|changed"):
            self.export(bundle=bundle)
        self.assertFalse((self.root / "input").exists())

    def test_changed_reference_evidence_is_rejected_before_output(self):
        mapping = self.mapping()
        receipt = mapping["references"]["tr|D1D2|D1D2_RABIT"]["sources"][0]
        changed = self.root / "changed_reference.json"
        changed.write_bytes(Path(receipt["absolute_path"]).read_bytes() + b"\n")
        receipt["absolute_path"] = str(changed)
        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            self.export(mapping=mapping)
        self.assertFalse((self.root / "input").exists())

    def cli_arguments(self, output, projection=None):
        return [
            "standardize",
            "--results", str(self.fixture / "csms.txt"),
            "--fasta", str(self.fixture / "search.fasta"),
            "--maxlynx-parameters", str(self.fixture / "parameters.xml"),
            "--dataset-accession", "PXD044630",
            "--reference-projection", str(
                self.fixture / "evidence" / "reference_mapping.json"
                if projection is None else projection
            ),
            "--crosslinker", "DSSO",
            "--residue", "K",
            "--export-unverified-candidates",
            "--output", str(output),
        ]

    def test_cli_reproduces_verified_api_input(self):
        api_output = self.export(name="api_input")
        output = self.root / "cli_input"
        result = CliRunner().invoke(
            cli.main, self.cli_arguments(output)
        )
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertEqual(
            (output / "claudio_input.csv").read_bytes(),
            (api_output / "claudio_input.csv").read_bytes(),
        )

    def test_cli_reuses_moved_reference_evidence(self):
        api_output = self.export(name="api_input")
        moved = self.root / "moved"
        api_output.rename(moved)
        output = self.root / "cli_input"
        result = CliRunner().invoke(
            cli.main,
            self.cli_arguments(
                output, projection=moved / "reference_projection.json"
            ),
        )
        self.assertEqual(result.exit_code, 0, result.output + repr(result.exception))
        self.assertEqual(
            (output / "claudio_input.csv").read_bytes(),
            (moved / "claudio_input.csv").read_bytes(),
        )
