import copy
import hashlib
import importlib
import importlib.util
import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


def object_sha256(value):
    return hashlib.sha256(json.dumps(
        value, ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")).hexdigest()


def local_tag(element):
    return element.tag.rsplit("}", 1)[-1]


class TestOpenMSLabeledMzIdentML(unittest.TestCase):
    def setUp(self):
        self.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures/openms_mzidentml_pxd014523"
        )
        self.expected = json.loads(
            (self.fixture / "expected.json").read_text(encoding="utf-8")
        )
        self.results = self.fixture / "OpenPepXL_BSA_PDH_CID.mzid"
        self.fasta = self.fixture / "search.fasta"

    def api(self):
        name = "claudio.dataset_catalog.openms_mzidentml"
        self.assertIsNotNone(
            importlib.util.find_spec(name),
            "Missing OpenMS labeled mzIdentML adapter",
        )
        module = importlib.import_module(name)
        self.assertTrue(callable(module.read_openms_labeled_mzidentml))
        return module

    def read_bundle(self, results=None, fasta=None):
        return self.api().read_openms_labeled_mzidentml(
            results or self.results,
            fasta or self.fasta,
            dataset_accession="PXD014523",
        )

    def test_import_retains_all_native_groups_and_items(self):
        bundle = self.read_bundle()
        self.assertEqual(bundle["format"], "mzidentml_csm")
        self.assertEqual(bundle["source_profile"], "openms_labeled")
        self.assertEqual(len(bundle["records"]), 70)
        self.assertEqual(len(bundle["source_identifications"]), 70)
        self.assertEqual(sum(
            len(record["source_items"]) for record in bundle["records"]
        ), 204)
        self.assertEqual(
            Counter(record["record_kind"] for record in bundle["records"]),
            {"crosslink_group": 32, "loop_link_group": 38},
        )
        self.assertTrue(all(
            record["record_status"] == "valid"
            for record in bundle["records"]
        ))
        self.assertEqual(bundle["parser"]["backend"], "native")
        self.assertIsNone(bundle["parser"].get("library"))

    def test_native_items_and_rank_zero_are_not_rewritten(self):
        bundle = self.read_bundle()
        for actual, expected in zip(
            bundle["records"], self.expected["groups"]
        ):
            for key in (
                "source_record_id", "source_record_number",
                "source_sir_attributes", "source_group_value",
                "source_rank_values",
            ):
                self.assertEqual(actual[key], expected[key])
            self.assertEqual(
                object_sha256(actual["source_items"]),
                expected["source_items_sha256"],
            )
        self.assertEqual(
            Counter(
                rank for record in bundle["records"]
                for rank in record["source_rank_values"]
            ),
            {"0": 204},
        )

    def test_light_and_heavy_observations_keep_both_native_ids(self):
        bundle = self.read_bundle()
        for actual, expected in zip(
            bundle["records"], self.expected["groups"]
        ):
            self.assertEqual(
                actual["source_native_ids"], expected["source_native_ids"]
            )
            self.assertEqual(len(actual["source_native_ids"]), 2)
            self.assertEqual(
                ",".join(actual["source_native_ids"]),
                actual["source_sir_attributes"]["spectrumID"],
            )
            self.assertEqual(
                actual["source_precursor_variants"],
                expected["source_precursor_variants"],
            )
            self.assertIsNone(actual.get("scan_number"))

    def test_zero_based_evidence_coordinates_are_verified_against_fasta(self):
        bundle = self.read_bundle()
        self.assertEqual(
            bundle["coordinate_policy"]["source_peptide_evidence_base"], 0
        )
        for actual, expected in zip(
            bundle["records"], self.expected["groups"]
        ):
            self.assertEqual(
                actual["normalized_sites"], expected["normalized_sites"]
            )
            for site in actual["normalized_sites"].values():
                sequence = bundle["proteins"][
                    site["source_protein_id"]
                ]["sequence"]
                self.assertEqual(
                    sequence[site["peptide_start"] - 1:site["peptide_end"]],
                    site["peptide"],
                )
                self.assertEqual(
                    site["protein_position"],
                    site["peptide_start"] + site["link_position"] - 1,
                )

        root = ET.fromstring(self.results.read_bytes())
        native = {
            element.attrib["id"]: element
            for element in root.iter()
            if local_tag(element) == "DBSequence"
        }
        copied = {
            item["attributes"]["id"]: item
            for item in bundle["source_sequence_collection"]["db_sequences"]
        }
        self.assertEqual(set(native), set(copied))
        for identifier, element in native.items():
            self.assertEqual(copied[identifier]["attributes"], element.attrib)
            self.assertFalse(any(
                local_tag(child) == "Seq" for child in element
            ))
            self.assertIsNone(copied[identifier]["sequence"])

    def test_negative_native_fdr_is_preserved_without_clamping(self):
        bundle = self.read_bundle()
        self.assertEqual(
            Counter(
                record["native_fdr"]["status"]
                for record in bundle["records"]
            ),
            {"valid": 63, "out_of_range": 7},
        )
        for actual, expected in zip(
            bundle["records"], self.expected["groups"]
        ):
            self.assertEqual(
                actual["native_fdr"]["raw_values"],
                expected["native_fdr_raw_values"],
            )
            if any(
                float(value) < 0 for value in expected["native_fdr_raw_values"]
            ):
                self.assertEqual(
                    actual["native_fdr"]["status"], "out_of_range"
                )
                self.assertIsNone(actual["native_fdr"]["value"])
            self.assertIs(actual["native_fdr_scope_verified"], False)
        self.assertIs(bundle["fdr_filter_applied"], False)

    def test_loop_links_keep_one_peptide_and_two_isotope_channels(self):
        bundle = self.read_bundle()
        loops = [
            record for record in bundle["records"]
            if record["record_kind"] == "loop_link_group"
        ]
        self.assertEqual(len(loops), 38)
        for record in loops:
            donor = record["normalized_sites"]["donor"]
            acceptor = record["normalized_sites"]["acceptor"]
            self.assertEqual(donor["peptide_ref"], acceptor["peptide_ref"])
            self.assertEqual(len(record["source_peptide_refs"]), 1)
            self.assertEqual(len(record["source_items"]), 2)
            self.assertEqual(len(record["source_precursor_variants"]), 2)
            self.assertTrue(all(
                len(variant["item_ids"]) == 1
                for variant in record["source_precursor_variants"]
            ))

    def test_missing_source_protein_does_not_use_another_fasta_entry(self):
        with tempfile.TemporaryDirectory() as temporary:
            fasta = Path(temporary) / "decoy_only.fasta"
            fasta.write_bytes((self.fixture / "BSA_rev.fasta").read_bytes())
            bundle = self.read_bundle(fasta=fasta)
        self.assertEqual(len(bundle["records"]), 70)
        self.assertTrue(all(
            record["record_status"] != "valid"
            for record in bundle["records"]
        ))
        self.assertEqual(sum(
            len(record["source_items"]) for record in bundle["records"]
        ), 204)

    def test_missing_heavy_chain_member_is_retained_as_unresolved(self):
        expected = next(
            record for record in self.expected["groups"]
            if record["record_kind"] == "crosslink_group"
        )
        item_to_remove = expected["source_precursor_variants"][1]["item_ids"][0]
        root = ET.fromstring(self.results.read_bytes())
        removed = 0
        for element in root.iter():
            for child in list(element):
                if child.attrib.get("id") == item_to_remove:
                    element.remove(child)
                    removed += 1
        self.assertEqual(removed, 1)

        with tempfile.TemporaryDirectory() as temporary:
            changed = Path(temporary) / "changed.mzid"
            changed.write_bytes(ET.tostring(root, encoding="utf-8"))
            bundle = self.read_bundle(results=changed)

        self.assertEqual(len(bundle["records"]), 70)
        self.assertEqual(sum(
            len(record["source_items"]) for record in bundle["records"]
        ), 203)
        record = next(
            item for item in bundle["records"]
            if item["source_sir_attributes"]["id"]
            == expected["source_sir_attributes"]["id"]
            and item["source_group_value"] == expected["source_group_value"]
        )
        self.assertEqual(record["record_status"], "not_evaluable")
        self.assertTrue(record["source_errors"])

    def test_placeholder_crosslinker_term_is_not_renamed(self):
        bundle = self.read_bundle()
        for actual, expected in zip(
            bundle["records"], self.expected["groups"]
        ):
            self.assertEqual(
                actual["source_crosslinker_terms"],
                expected["source_crosslinker_terms"],
            )
            self.assertIs(actual["crosslinker_identity_verified"], False)
        self.assertTrue(any(
            term.get("accession") == "XLMOD:XXXXX"
            and term.get("name") == "DSS"
            for record in bundle["records"]
            for term in record["source_crosslinker_terms"]
        ))

class TestOpenMSLabeledAuditCLI(unittest.TestCase):
    def setUp(self):
        self.fixture = (
            Path(__file__).resolve().parents[1]
            / "fixtures/openms_mzidentml_pxd014523"
        )
        self.results = self.fixture / "OpenPepXL_BSA_PDH_CID.mzid"
        self.fasta = self.fixture / "search.fasta"
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def arguments(self, destination):
        return [
            "standardize",
            "--results", str(self.results),
            "--fasta", str(self.fasta),
            "--dataset-accession", "PXD014523",
            "--mzidentml-profile", "openms-labeled",
            "--output", str(destination),
        ]

    def test_cli_audit_matches_the_native_profile(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main
        from claudio.dataset_catalog.openms_mzidentml import (
            read_openms_labeled_mzidentml,
        )

        destination = self.root / "audit"
        result = CliRunner().invoke(
            main, self.arguments(destination) + ["--audit-only"]
        )
        self.assertEqual(
            result.exit_code, 0, result.output + repr(result.exception)
        )
        report = json.loads(
            (destination / "import_report.json").read_text(encoding="utf-8")
        )
        expected = read_openms_labeled_mzidentml(
            self.results, self.fasta, dataset_accession="PXD014523"
        )
        self.assertEqual(report, expected)
        self.assertEqual(len(report["records"]), 70)
        self.assertEqual(sum(
            len(record["source_items"]) for record in report["records"]
        ), 204)
        self.assertEqual(
            Counter(record["native_fdr"]["status"]
                    for record in report["records"]),
            {"valid": 63, "out_of_range": 7},
        )
        manifest = json.loads(
            (destination / "manifest.json").read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["counts"]["source_records"], 70)
        self.assertEqual(
            manifest["counts"]["validated_source_records"], 70
        )
        self.assertIs(manifest["fdr_filter_applied"], False)
        self.assertIs(manifest["input_exported"], False)
        self.assertFalse((destination / "claudio_input.csv").exists())

    def test_profile_requires_audit_mode_and_rejects_selection_options(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main

        cases = [
            ([], "--audit-only"),
            (["--audit-only", "--residue", "D"], "--residue"),
            (["--audit-only", "--crosslinker", "PDH"], "--crosslinker"),
            (
                ["--audit-only", "--fdr-level", "psm",
                 "--fdr-cutoff", "0.05"],
                "--fdr-level",
            ),
            (
                ["--audit-only", "--export-unverified-candidates"],
                "--export-unverified-candidates",
            ),
        ]
        for number, (options, expected) in enumerate(cases):
            with self.subTest(options=options):
                destination = self.root / f"rejected_{number}"
                result = CliRunner().invoke(
                    main, self.arguments(destination) + options
                )
                self.assertNotEqual(result.exit_code, 0)
                self.assertIn(expected, result.output)
                self.assertFalse(destination.exists())

    def test_profile_option_rejects_other_result_formats(self):
        from click.testing import CliRunner
        from claudio.dataset_catalog.cli import main

        other = (
            self.fixture.parent / "kojak_pxd015044"
            / "20160620_TauHis_TDIG_monomer_CBDPS-H8-D81mM_OT_MT_2.kojak.txt"
        )
        self.assertTrue(other.is_file())
        destination = self.root / "wrong_format"
        arguments = self.arguments(destination) + ["--audit-only"]
        arguments[arguments.index("--results") + 1] = str(other)
        result = CliRunner().invoke(main, arguments)
        self.assertNotEqual(result.exit_code, 0)
        self.assertIn("--mzidentml-profile", result.output)
        self.assertFalse(destination.exists())

    def test_changed_sources_are_rejected_before_audit_output(self):
        from claudio.dataset_catalog.openms_mzidentml import (
            read_openms_labeled_mzidentml,
        )
        from claudio.dataset_catalog.xifdr import write_import_report

        results = self.root / "source.mzid"
        fasta = self.root / "search.fasta"
        results.write_bytes(self.results.read_bytes())
        fasta.write_bytes(self.fasta.read_bytes())
        bundle = read_openms_labeled_mzidentml(
            results, fasta, dataset_accession="PXD014523"
        )
        for name, path in (("psm", results), ("fasta", fasta)):
            with self.subTest(source=name):
                original = path.read_bytes()
                destination = self.root / ("changed_" + name)
                try:
                    path.write_bytes(original + b"\n")
                    with self.assertRaisesRegex(
                        ValueError, "(?i)checksum|changed"
                    ):
                        write_import_report(bundle, destination)
                    self.assertFalse(destination.exists())
                finally:
                    path.write_bytes(original)
