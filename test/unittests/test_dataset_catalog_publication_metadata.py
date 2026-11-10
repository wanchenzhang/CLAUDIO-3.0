"""Offline tests for dataset-scoped publication evidence."""

from copy import deepcopy
import hashlib
import importlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from claudio.dataset_catalog import published_discovery
from claudio.dataset_catalog.discovery_history import initialize_discovery_history
from claudio.dataset_catalog.publication_lookup import PublicationLookup
from test_dataset_catalog_publication_lookup import Clock
from test_dataset_catalog_published_discovery import ProviderSession, page, record


def encoded(value):
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode()


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


class TestDatasetPublicationMetadata(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.calls = []
        guard = patch(
            "requests.sessions.Session.request",
            side_effect=AssertionError("Live HTTP is forbidden in this test"),
        )
        guard.start()
        self.addCleanup(guard.stop)
        records = [
            record(),
            record("PXD000002", '<a href="https://doi.org/10.1234/other">Paper</a>'),
        ]
        self.discovery = self.make_discovery(records, "discovery")
        self.candidates = self.discovery / "candidates.json"

    def make_discovery(self, records, name):
        history = self.root / (name + "_history")
        initialize_discovery_history(history)

        def factory(directory, **options):
            clock = Clock()
            return PublicationLookup(
                directory, session=ProviderSession(self.calls),
                clock=clock, sleeper=clock.sleep, **options,
            )

        with patch.object(
            published_discovery, "fetch_dataset_pages",
            side_effect=lambda term, *, page_size: [page(term, records, page_size)],
        ), patch.object(published_discovery, "PublicationLookup", side_effect=factory):
            return published_discovery.run_published_discovery(
                ["XL-MS"], self.root / name, history_dir=history,
                as_of="2026-11-09",
            )

    def api(self):
        name = "claudio.dataset_catalog.publication_metadata"
        try:
            return importlib.import_module(name)
        except ModuleNotFoundError as error:
            if error.name != name:
                raise
            self.fail("Dataset publication metadata loader is not implemented")

    def load(self, path=None, accession="PXD000001"):
        return self.api().load_dataset_metadata(
            path or self.candidates, dataset_accession=accession,
        )

    def read(self, path):
        return json.loads(Path(path).read_bytes())

    def rewrite_output(self, path, document):
        path.write_bytes(encoded(document))
        manifest_path = path.parent / "manifest.json"
        manifest = self.read(manifest_path)
        manifest["output_sha256"][path.name] = digest(path.read_bytes())
        if path.name == "dataset_metadata.json":
            manifest["dataset_metadata"]["sha256"] = digest(path.read_bytes())
        manifest_path.write_bytes(encoded(manifest))

    def bundle(self):
        return {
            "provenance": {"dataset_accession": "PXD000001"},
            "records": [{
                "source_record_id": "PXD000001:source:000001",
                "fdr": {"psm": "0.02"},
            }],
            "proteins": {},
            "fdr_policy": {"verification_status": "unverified"},
        }

    def snapshot(self):
        api = self.api()
        bound = api.bind_dataset_metadata(
            self.bundle(), self.candidates,
            dataset_accession="PXD000001",
        )
        prepared = api.prepare_publication_export(bound)
        destination = self.root / "export"
        destination.mkdir()
        for name, raw in prepared["source_bodies"].items():
            target = destination / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(raw)
        (destination / "manifest.json").write_bytes(encoded({
            "schema_version": 1,
            "kind": "synthetic_claudio_input",
            "status": "complete",
            "dataset_metadata": prepared["manifest_entry"],
            "output_sha256": {
                name: digest(raw)
                for name, raw in prepared["source_bodies"].items()
            },
        }))
        return destination / "dataset_metadata.json", prepared

    def test_loads_selected_paper_from_discovery(self):
        packet = self.load()
        self.assertEqual(packet["record"]["dataset_accession"], "PXD000001")
        self.assertEqual(packet["record"]["publication_status"], "published")
        paper = packet["record"]["papers"][0]
        self.assertEqual(paper["doi"], "10.1234/example")
        self.assertEqual(paper["title"], "An example article")
        self.assertTrue(paper["verification_evidence"])

    def test_merges_independently_verified_doi_and_pmid(self):
        source = self.make_discovery([record(publications=(
            '<a href="https://doi.org/10.1234/example">DOI</a>'
            '<a href="https://pubmed.ncbi.nlm.nih.gov/1234/">PMID</a>'
        ))], "paired_ids")
        packet = self.load(source / "candidates.json")
        self.assertEqual(len(packet["record"]["papers"]), 1)
        paper = packet["record"]["papers"][0]
        self.assertEqual(
            (paper["doi"], paper["pmid"]),
            ("10.1234/example", "1234"),
        )
        self.assertEqual(len(paper["verification_evidence"]), 2)

    def test_rejects_dataset_not_emitted_in_discovery(self):
        with self.assertRaises(ValueError):
            self.load(accession="PXD999999")

    def test_rejects_changed_candidate_bytes(self):
        self.candidates.write_bytes(self.candidates.read_bytes() + b" ")
        with self.assertRaises(ValueError):
            self.load()

    def test_rejects_changed_http_body(self):
        proof = self.read(self.candidates)[0]["papers"][0]["verification_evidence"][0]
        receipt = self.read(
            self.discovery / "publication_evidence" / proof["receipt_path"]
        )
        body = (
            self.discovery / "publication_evidence"
            / receipt["attempts"][-1]["body_path"]
        )
        body.write_bytes(body.read_bytes() + b" ")
        with self.assertRaises(ValueError):
            self.load()

    def test_rejects_paper_without_verification_evidence(self):
        candidates = self.read(self.candidates)
        candidates[0]["papers"][0]["verification_evidence"] = []
        self.rewrite_output(self.candidates, candidates)
        with self.assertRaises(ValueError):
            self.load()

    def test_rejects_unpublished_candidate_status(self):
        candidates = self.read(self.candidates)
        candidates[0]["publication_status"] = "pending"
        self.rewrite_output(self.candidates, candidates)
        with self.assertRaises(ValueError):
            self.load()

    def test_packages_only_selected_paper_http_evidence(self):
        packet = self.load()
        receipts = [
            self.read(
                self.discovery / "publication_evidence" / proof["receipt_path"]
            )
            for proof in self.read(self.candidates)[1]["papers"][0]["verification_evidence"]
        ]
        for receipt in receipts:
            for attempt in receipt["attempts"]:
                self.assertNotIn(
                    "publication_sources/http/" + attempt["body_path"],
                    packet["source_bodies"],
                )
        self.assertEqual(len(packet["record"]["papers"]), 1)

    def test_loading_does_not_change_source_files(self):
        before = {
            str(path.relative_to(self.discovery)): digest(path.read_bytes())
            for path in self.discovery.rglob("*") if path.is_file()
        }
        self.load()
        after = {
            str(path.relative_to(self.discovery)): digest(path.read_bytes())
            for path in self.discovery.rglob("*") if path.is_file()
        }
        self.assertEqual(after, before)

    def test_packaged_metadata_reloads_without_original_discovery(self):
        path, prepared = self.snapshot()
        self.discovery.rename(self.root / "original_moved")
        packet = self.load(path)
        self.assertEqual(packet["record"], prepared["record"])

    def test_packaged_metadata_rejects_wrong_dataset(self):
        path, _ = self.snapshot()
        with self.assertRaises(ValueError):
            self.load(path, "PXD000002")

    def test_rejects_forged_paper_title_even_with_updated_outer_checksum(self):
        path, _ = self.snapshot()
        document = self.read(path)
        document["papers"][0]["title"] = "An unrelated article"
        self.rewrite_output(path, document)
        with self.assertRaises(ValueError):
            self.load(path)

    def test_rejects_publication_evidence_path_escape(self):
        path, _ = self.snapshot()
        document = self.read(path)
        document["papers"][0]["verification_evidence"][0]["root"] = "../outside"
        self.rewrite_output(path, document)
        with self.assertRaises(ValueError):
            self.load(path)

    def test_binding_preserves_source_records_and_confidence(self):
        original = self.bundle()
        before = deepcopy(original)
        bound = self.api().bind_dataset_metadata(
            original, self.candidates, dataset_accession="PXD000001",
        )
        prepared = self.api().prepare_publication_export(bound)
        self.assertEqual(original, before)
        self.assertEqual(bound["records"], before["records"])
        self.assertEqual(bound["fdr_policy"], before["fdr_policy"])
        fields = prepared["row_fields"]
        self.assertEqual(fields["xlms_dataset_accession"], "PXD000001")
        self.assertEqual(
            json.loads(fields["xlms_paper_dois_json"]),
            ["10.1234/example"],
        )
        self.assertEqual(
            json.loads(fields["xlms_papers_json"]),
            prepared["record"]["papers"],
        )

    def test_rejects_bundle_dataset_conflict(self):
        bundle = self.bundle()
        bundle["provenance"]["dataset_accession"] = "PXD000002"
        with self.assertRaises(ValueError):
            self.api().bind_dataset_metadata(
                bundle, self.candidates, dataset_accession="PXD000001",
            )

    def test_rejects_foreign_source_record_dataset(self):
        bundle = self.bundle()
        bundle["records"][0]["source_record_id"] = "PXD000002:source:000001"
        with self.assertRaises(ValueError):
            self.api().bind_dataset_metadata(
                bundle, self.candidates, dataset_accession="PXD000001",
            )

    def test_rejects_metadata_changed_after_binding(self):
        bound = self.api().bind_dataset_metadata(
            self.bundle(), self.candidates,
            dataset_accession="PXD000001",
        )
        self.candidates.write_bytes(self.candidates.read_bytes() + b" ")
        with self.assertRaises(ValueError):
            self.api().prepare_publication_export(bound)

    def test_rejects_changed_packaged_http_body(self):
        path, _ = self.snapshot()
        body = next(
            (path.parent / "publication_sources/http/bodies").iterdir()
        )
        body.write_bytes(body.read_bytes() + b" ")
        with self.assertRaises(ValueError):
            self.load(path)

    def test_legacy_bundle_without_metadata_is_unchanged(self):
        bundle = self.bundle()
        self.assertIs(
            self.api().bind_dataset_metadata(
                bundle, None, dataset_accession=None
            ),
            bundle,
        )
        self.assertIsNone(self.api().prepare_publication_export(bundle))

    def test_rejects_legacy_seed_as_a_new_discovery_context(self):
        path = self.discovery / "manifest.json"
        manifest = self.read(path)
        manifest["operation"] = "legacy_history_seed"
        path.write_bytes(encoded(manifest))
        with self.assertRaises(ValueError):
            self.load()

    def test_rejects_changed_packaged_registry_association(self):
        path, _ = self.snapshot()
        document = self.read(path)
        document["source_candidate"]["publications"] = "10.1234/unrelated"
        self.rewrite_output(path, document)
        with self.assertRaises(ValueError):
            self.load(path)


if __name__ == "__main__":
    unittest.main()
