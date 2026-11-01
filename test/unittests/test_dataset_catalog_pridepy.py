import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import requests

from claudio.dataset_catalog import repositories
from claudio.dataset_catalog.inventory import inspect_candidate, inspect_catalog, summarize_inventories


def registry_response(accession, repository, native):
    code = {"MassIVE": "MS:1002487", "jPOST": "MS:1002632", "iProX": "MS:1002836"}[repository]
    data = {"datasetSummary": {"hostingRepository": repository}, "identifiers": [
        {"accession": "MS:1001919", "value": accession},
        {"accession": code, "value": native},
    ]}
    response = requests.Response()
    response.status_code = 200
    response.url = f"https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets/{accession}?resultType=extended"
    response.headers["Content-Type"] = "application/json"
    response._content = json.dumps(data).encode()
    return response


def provider_file(native, name="cross-linked_sites.csv", size=16, folder="results"):
    host = {"MSV": "massive.ucsd.edu", "JPS": "ftp.jpostdb.org", "IPX": "download.iprox.org"}[native[:3]]
    return {"accession": native, "fileName": name, "relativePath": f"{folder}/{name}",
            "fileCategory": {"value": "RESULT"}, "fileSizeBytes": size,
            "publicFileLocations": [{"name": "FTP Protocol", "value": f"ftp://{host}/{native}/{folder}/{name}"}]}


class TestPridepyInventories(unittest.TestCase):
    def setUp(self):
        self.version = patch.object(repositories, "_pridepy_version", return_value="0.0.16").start()
        self.provider = Mock()
        self.provider_patch = patch.object(repositories, "_pridepy_provider", return_value=self.provider).start()
        self.addCleanup(patch.stopall)

    def candidate(self, repository="PRIDE"):
        return {"accession": "PXD000001", "repository": repository, "review_status": "pending", "search_terms": ["XL-MS"]}

    @patch.object(repositories, "_pridepy_project")
    def test_pride_uses_library_and_records_export_provenance(self, project):
        project.return_value = {"accession": "PXD000001", "title": "Example"}
        records = [{"accession": "file-1", "projectAccessions": ["PXD000001"],
                    "fileName": "cross-linked_sites.csv", "fileSizeBytes": 16,
                    "publicFileLocations": [{"name": "FTP", "value": "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/PXD000001/sites.csv"}]}]
        self.provider.list_files.return_value = records
        original = deepcopy(records)
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "inventory"
            candidate = self.candidate()
            result = inspect_candidate(candidate, target)
            self.assertEqual(result["inspection_status"], "collected", result)
            manifest = json.loads((target / "manifest.json").read_bytes())
            self.assertEqual(manifest["backend"], "pridepy")
            self.assertEqual(manifest["backend_version"], "0.0.16")
            self.assertEqual(manifest["file_responses"][0]["source_kind"], "library_export")
            self.assertIsNone(manifest["reported_file_count"])
            self.assertEqual(json.loads((target / "files.json").read_bytes())[0]["accession"], "file-1")
            raw = manifest["file_responses"][0]
            self.assertEqual(hashlib.sha256((target / raw["path"]).read_bytes()).hexdigest(), raw["sha256"])
            self.assertEqual(result["search_terms"], candidate["search_terms"])
            self.assertEqual(result["review_status"], "pending")
        self.provider.list_files.assert_called_once_with("PXD000001")
        self.assertEqual(records, original)

    @patch.object(repositories.requests, "get")
    def test_native_repositories_resolve_pxd_before_library_listing(self, get):
        for repository, native in [("MassIVE", "MSV000079237"), ("jPOST", "JPST000391"), ("iProX", "IPX0001609000")]:
            with self.subTest(repository=repository), TemporaryDirectory() as tmp:
                get.return_value = registry_response("PXD000001", repository, native)
                self.provider.list_files.return_value = [provider_file(native)]
                self.provider.list_files.reset_mock()
                result = inspect_candidate(self.candidate(repository), Path(tmp) / "inventory")
                self.assertEqual(result["inspection_status"], "collected", result)
                self.provider.list_files.assert_called_once_with(native)
                manifest = json.loads((Path(result["inventory_path"]) / "manifest.json").read_bytes())
                self.assertEqual(manifest["repository_accession"], native)
                self.assertEqual(manifest["project_response"]["source_kind"], "http_response")

    @patch.object(repositories.requests, "get")
    def test_dataset_accession_is_not_used_as_every_file_identity(self, get):
        native = "JPST000391"
        get.return_value = registry_response("PXD000001", "jPOST", native)
        self.provider.list_files.return_value = [provider_file(native, folder="rep1"), provider_file(native, folder="rep2")]
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "inventory"
            result = inspect_candidate(self.candidate("jPOST"), target)
            self.assertEqual(result["inspection_status"], "collected", result)
            records = json.loads((target / "files.json").read_bytes())
            self.assertEqual(len({r["accession"] for r in records}), 2)
            self.assertEqual({r["repositoryAccessions"][0] for r in records}, {native})
            self.assertEqual({r["projectAccessions"][0] for r in records}, {"PXD000001"})

    @patch.object(repositories.requests, "get")
    def test_identical_duplicate_links_keep_their_provider_indices(self, get):
        native = "JPST000391"
        get.return_value = registry_response("PXD000001", "jPOST", native)
        record = provider_file(native)
        self.provider.list_files.return_value = [record, deepcopy(record)]
        with TemporaryDirectory() as tmp:
            result = inspect_candidate(self.candidate("jPOST"), Path(tmp) / "inventory")
            self.assertEqual(result["inspection_status"], "collected", result)
            saved = json.loads((Path(result["inventory_path"]) / "files.json").read_bytes())
            self.assertEqual(len(saved), 1)
            self.assertEqual(saved[0]["providerRecordIndices"], [0, 1])

    @patch.object(repositories.requests, "get")
    def test_conflicting_duplicate_links_fail_without_creating_snapshot(self, get):
        native = "JPST000391"
        get.return_value = registry_response("PXD000001", "jPOST", native)
        first = provider_file(native)
        second = deepcopy(first)
        second["fileSizeBytes"] = 17
        self.provider.list_files.return_value = [first, second]
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "inventory"
            result = inspect_candidate(self.candidate("jPOST"), target)
            self.assertEqual(result["inspection_status"], "failed")
            self.assertIn("conflicting", result["error_message"].lower())
            self.assertFalse(target.exists())

    @patch.object(repositories.requests, "get")
    def test_registry_mismatch_does_not_list_another_dataset(self, get):
        get.return_value = registry_response("PXD999999", "jPOST", "JPST000391")
        with TemporaryDirectory() as tmp:
            result = inspect_candidate(self.candidate("jPOST"), Path(tmp) / "inventory")
            self.assertEqual(result["inspection_status"], "failed")
        self.provider.list_files.assert_not_called()

    @patch.object(repositories.requests, "get")
    def test_wrong_hosting_repository_is_not_silently_routed(self, get):
        get.return_value = registry_response("PXD000001", "MassIVE", "MSV000079237")
        with TemporaryDirectory() as tmp:
            result = inspect_candidate(self.candidate("jPOST"), Path(tmp) / "inventory")
            self.assertEqual(result["inspection_status"], "failed")
        self.provider.list_files.assert_not_called()

    @patch.object(repositories, "_pridepy_project", return_value={"accession": "PXD000001"})
    def test_none_from_library_is_a_failure_not_an_empty_success(self, project):
        self.provider.list_files.return_value = None
        with TemporaryDirectory() as tmp:
            result = inspect_candidate(self.candidate(), Path(tmp) / "inventory")
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["review_status"], "pending")

    @patch.object(repositories, "_pridepy_project", return_value={"accession": "PXD000001"})
    def test_library_exception_is_recorded_without_excluding_dataset(self, project):
        self.provider.list_files.side_effect = RuntimeError("listing unavailable")
        with TemporaryDirectory() as tmp:
            result = inspect_candidate(self.candidate(), Path(tmp) / "inventory")
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_type"], "RuntimeError")
            self.assertEqual(result["review_status"], "pending")

    @patch.object(repositories, "_pridepy_project", return_value={"accession": "PXD000001"})
    def test_library_interrupt_is_not_swallowed(self, project):
        self.provider.list_files.side_effect = KeyboardInterrupt()
        with TemporaryDirectory() as tmp, self.assertRaises(KeyboardInterrupt):
            inspect_candidate(self.candidate(), Path(tmp) / "inventory")

    def test_unknown_repository_remains_pending(self):
        with TemporaryDirectory() as tmp:
            result = inspect_candidate(self.candidate("Other"), Path(tmp) / "inventory")
            self.assertEqual(result["inspection_status"], "adapter_pending")
        self.provider.list_files.assert_not_called()

    @patch.object(repositories, "_pridepy_project", return_value={"accession": "PXD000001"})
    def test_resume_and_summary_read_saved_exports_without_network(self, project):
        self.provider.list_files.return_value = [{"accession": "file-1", "projectAccessions": ["PXD000001"],
            "fileName": "cross-linked_sites.csv", "fileSizeBytes": 16, "fileCategory": {"value": "RESULT"},
            "publicFileLocations": [{"name": "FTP", "value": "ftp://ftp.pride.ebi.ac.uk/PXD000001/sites.csv"}]}]
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            catalog = root / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]))
            original = inspect_catalog(catalog, root / "run")
            self.provider.list_files.reset_mock()
            self.assertEqual(inspect_catalog(catalog, root / "run", resume=True), original)
            summarize_inventories(catalog, [root / "run"], root / "summary")
            self.provider.list_files.assert_not_called()
            report = json.loads((root / "summary/manifest.json").read_bytes())
            self.assertEqual(report["status_counts"]["collected"], 1)

    @patch.object(repositories, "_pridepy_project", return_value={"accession": "PXD000001"})
    def test_resume_rejects_changed_provider_export(self, project):
        self.provider.list_files.return_value = []
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            catalog = root / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]))
            results = inspect_catalog(catalog, root / "run")
            manifest_path = root / "run" / results[0]["inventory_path"] / "manifest.json"
            manifest = json.loads(manifest_path.read_bytes())
            export = manifest_path.parent / manifest["file_responses"][0]["path"]
            export.write_text("{}")
            with self.assertRaises(ValueError):
                inspect_catalog(catalog, root / "run", resume=True)
