import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import requests

from claudio.dataset_catalog.inventory import inspect_candidate, inspect_catalog
from claudio.dataset_catalog.repositories import (
    RepositoryResponse, fetch_panorama_project, fetch_panorama_file_pages,
    write_repository_outputs,
)


ORIGIN = "https://panoramaweb.org"
QUERY_URL = ORIGIN + "/Panorama%20Public/query-selectRows.api"
ROOT = ORIGIN + "/_webdav/Panorama%20Public/2020/Study/"


def metadata():
    return {
        "schemaName": "panoramapublic", "queryName": "ExperimentAnnotations",
        "rowCount": 1,
        "rows": [{
            "pxid": "PXD000001", "Id": 123, "Title": "Cross-linking dataset",
            "_labkeyurl_Title": "/Panorama%20Public/2020/Study/project-begin.view",
        }],
    }


def directory(url):
    return {"href": url, "collection": True, "canRead": True}


def file_entry(url, size):
    return {
        "href": url, "collection": False, "canRead": True,
        "contentlength": size, "etag": 'W/"size-time"',
        "dataFileUrl": "file:///server/private/path/never-use-this",
    }


def listings():
    return {
        ROOT: {"files": [
            directory(ROOT + "@files/"), directory(ROOT + "replicate/"),
            directory(ROOT + "@pipeline/"),
        ]},
        ROOT + "@files/": {"files": [
            directory(ROOT + "@files/raw/"), file_entry(ROOT + "@files/results.tsv", 0),
        ]},
        ROOT + "@files/raw/": {"files": [file_entry(ROOT + "@files/raw/run.raw", 123)]},
        ROOT + "replicate/": {"files": [directory(ROOT + "replicate/@files/")]},
        ROOT + "replicate/@files/": {"files": [file_entry(ROOT + "replicate/@files/run.raw", 456)]},
    }


def response(url, data, params=None):
    result = requests.Response()
    result.status_code = 200
    result.url = requests.Request("GET", url, params=params).prepare().url
    result.headers["Content-Type"] = "application/json"
    result._content = (json.dumps(data) + "\n\n").encode("utf-8")
    return result


def snapshot(url, data, params=None):
    raw = response(url, data, params)
    return RepositoryResponse(
        raw.url, "2026-10-31T17:00:00+00:00", raw.content,
        raw.headers.copy(), data,
    )


class TestPanoramaInventory(unittest.TestCase):
    def candidate(self):
        return {
            "accession": "PXD000001", "repository": "PanoramaPublic",
            "title": "Test dataset", "search_terms": ["cross-linking"],
            "review_status": "pending",
        }

    def inputs(self):
        return snapshot(QUERY_URL, metadata()), [
            snapshot(url, data, {"method": "json"})
            for url, data in listings().items()
        ]

    def respond(self, url, **kwargs):
        if url == QUERY_URL:
            return response(url, metadata(), kwargs.get("params"))
        data = listings().get(url)
        if data is None:
            raise AssertionError(f"Unexpected request: {url}")
        return response(url, data, kwargs.get("params"))

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_finds_the_current_project_by_pxd_and_preserves_the_query_response(self, mock_get):
        mock_get.side_effect = self.respond
        project = fetch_panorama_project("PXD000001")
        self.assertEqual(project.data, metadata())
        self.assertEqual(json.loads(project.raw_body), metadata())
        params = mock_get.call_args.kwargs["params"]
        self.assertEqual(params["query.pxid~eq"], "PXD000001")
        self.assertEqual(params["query.Version~contains"], "Current")
        self.assertEqual(params["containerFilter"], "CurrentAndSubfolders")
        self.assertEqual(params["maxRows"], 2)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_missing_wrong_ambiguous_or_unsafe_project_identity(self, mock_get):
        for problem in ("empty", "wrong_pxd", "multiple", "outside", "invalid_id"):
            with self.subTest(problem=problem):
                data = metadata()
                if problem == "empty":
                    data["rows"] = []
                    data["rowCount"] = 0
                elif problem == "wrong_pxd":
                    data["rows"][0]["pxid"] = "PXD000002"
                elif problem == "multiple":
                    data["rows"] *= 2
                    data["rowCount"] = 2
                elif problem == "outside":
                    data["rows"][0]["_labkeyurl_Title"] = "https://example.org/project-begin.view"
                else:
                    data["rows"][0]["Id"] = None
                mock_get.return_value = response(QUERY_URL, data)
                with self.assertRaises(ValueError):
                    fetch_panorama_project("PXD000001")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_invalid_accessions_before_requesting_data(self, mock_get):
        for accession in (None, "", "PXD1", "PXD000001/other"):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError):
                    fetch_panorama_project(accession)
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_reads_file_stores_and_subprojects_without_following_server_local_paths(self, mock_get):
        mock_get.side_effect = self.respond
        project, _ = self.inputs()

        pages = fetch_panorama_file_pages("PXD000001", project)

        self.assertEqual(len(pages), 5)
        self.assertEqual({url.split("?", 1)[0] for url in [p.request_url for p in pages]}, set(listings()))
        self.assertEqual(mock_get.call_count, 5)
        for call in mock_get.call_args_list:
            self.assertEqual(call.kwargs["params"], {"method": "json"})
            self.assertTrue(call.args[0].startswith(ROOT))

    def test_writes_unique_file_paths_sizes_and_raw_sources_without_using_etag_as_checksum(self):
        project, pages = self.inputs()
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("PanoramaPublic", "PXD000001", project, pages, output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            files = json.loads((output_dir / "files.json").read_bytes())
            self.assertEqual(manifest["repository_accession"], "123")
            self.assertEqual(manifest["file_inventory_scope"], "recursive_webdav_listing")
            self.assertEqual(manifest["file_count"], 3)
            self.assertIsNone(manifest["reported_file_count"])
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(len({f["accession"] for f in files}), 3)
            self.assertEqual(sorted(f["fileSizeBytes"] for f in files), [0, 123, 456])
            self.assertEqual(sum(f["fileName"] == "run.raw" for f in files), 2)
            for file in files:
                self.assertIsNone(file["checksum"])
                self.assertEqual(file["repositoryMetadata"]["etag"], 'W/"size-time"')
                self.assertTrue(file["publicFileLocations"][0]["value"].startswith(ROOT))
            for source in [manifest["project_response"], *manifest["file_responses"]]:
                body = (output_dir / source["path"]).read_bytes()
                self.assertEqual(hashlib.sha256(body).hexdigest(), source["sha256"])

    def test_requires_every_discovered_directory_response_before_writing(self):
        for problem in ("missing", "duplicate", "edited"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                project, pages = self.inputs()
                if problem == "missing":
                    pages.pop()
                elif problem == "duplicate":
                    pages.append(pages[0])
                else:
                    pages[0].data["files"].pop()
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("PanoramaPublic", "PXD000001", project, pages, output_dir)
                self.assertFalse(output_dir.exists())

    def test_accepts_an_empty_file_store_but_not_a_missing_file_store(self):
        for has_store in (True, False):
            with self.subTest(has_store=has_store), TemporaryDirectory() as temp_dir:
                project, _ = self.inputs()
                data = {"files": [directory(ROOT + "@files/")] if has_store else []}
                pages = [snapshot(ROOT, data, {"method": "json"})]
                if has_store:
                    pages.append(snapshot(ROOT + "@files/", {"files": []}, {"method": "json"}))
                output_dir = Path(temp_dir) / "project"
                if has_store:
                    write_repository_outputs("PanoramaPublic", "PXD000001", project, pages, output_dir)
                    manifest = json.loads((output_dir / "manifest.json").read_bytes())
                    self.assertEqual(manifest["file_count"], 0)
                else:
                    with self.assertRaises(ValueError):
                        write_repository_outputs("PanoramaPublic", "PXD000001", project, pages, output_dir)
                    self.assertFalse(output_dir.exists())

    def test_rejects_outside_paths_malformed_listings_and_bad_file_sizes(self):
        for problem in ("outside", "parent", "missing_files", "unknown_kind", "bad_size"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                project, pages = self.inputs()
                data = {"files": [file_entry(ROOT + "@files/run.raw", 1)]}
                if problem == "outside":
                    data["files"][0]["href"] = "https://example.org/run.raw"
                elif problem == "parent":
                    data["files"][0]["href"] = ROOT + "@files/%2e%2e/run.raw"
                elif problem == "missing_files":
                    data = {}
                elif problem == "unknown_kind":
                    data["files"][0].pop("collection")
                else:
                    data["files"][0]["contentlength"] = -1
                pages[1] = snapshot(ROOT + "@files/", data, {"method": "json"})
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("PanoramaPublic", "PXD000001", project, pages, output_dir)
                self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_collects_a_candidate_and_preserves_its_review_status(self, mock_get):
        mock_get.side_effect = self.respond
        candidate = self.candidate()
        before = deepcopy(candidate)
        with TemporaryDirectory() as temp_dir:
            result = inspect_candidate(candidate, Path(temp_dir) / "project")
            self.assertEqual(result["inspection_status"], "collected")
            self.assertEqual(result["file_count"], 3)
            self.assertEqual(result["review_status"], "pending")
            self.assertEqual(candidate, before)
            self.assertEqual(mock_get.call_count, 6)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_records_rate_limiting_as_a_failure_instead_of_an_empty_inventory(self, mock_get):
        limited = response(ROOT, {"error": "Too many requests"})
        limited.status_code = 429
        mock_get.side_effect = [response(QUERY_URL, metadata()), limited]
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            result = inspect_candidate(self.candidate(), output_dir)
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "file_list")
            self.assertIsNone(result["file_count"])
            self.assertEqual(result["review_status"], "pending")
            self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.time.sleep")
    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_honors_a_short_retry_after_and_stops_after_three_attempts(self, mock_get, mock_sleep):
        limited = response(QUERY_URL, {"error": "Too many requests"})
        limited.status_code = 429
        limited.headers["Retry-After"] = "1"
        mock_get.side_effect = [limited, response(QUERY_URL, metadata())]

        project = fetch_panorama_project("PXD000001")

        self.assertEqual(project.data, metadata())
        self.assertEqual(mock_get.call_count, 2)
        mock_sleep.assert_called_once_with(1)
        mock_get.reset_mock()
        mock_sleep.reset_mock()
        mock_get.side_effect = [limited, limited, limited]
        with self.assertRaises(requests.exceptions.HTTPError):
            fetch_panorama_project("PXD000001")
        self.assertEqual(mock_get.call_count, 3)
        self.assertEqual(mock_sleep.call_count, 2)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_retains_a_candidate_when_a_required_directory_is_not_readable(self, mock_get):
        data = {"files": [dict(directory(ROOT + "@files/"), canRead=False)]}
        mock_get.side_effect = [response(QUERY_URL, metadata()), response(ROOT, data)]
        with TemporaryDirectory() as temp_dir:
            result = inspect_candidate(self.candidate(), Path(temp_dir) / "project")
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "file_list")
            self.assertEqual(result["review_status"], "pending")
            self.assertIsNone(result["file_count"])

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_an_error_response_even_when_it_contains_an_empty_files_array(self, mock_get):
        root_data = {"files": [directory(ROOT + "@files/")]}
        failed = {"files": [], "success": False, "exception": "Listing unavailable"}
        mock_get.side_effect = [
            response(QUERY_URL, metadata()), response(ROOT, root_data),
            response(ROOT + "@files/", failed),
        ]
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            result = inspect_candidate(self.candidate(), output_dir)
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "file_list")
            self.assertIsNone(result["file_count"])
            self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_resumes_complete_directory_sources_and_rejects_a_changed_native_identity(self, mock_get):
        mock_get.side_effect = self.respond
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            results = inspect_catalog(catalog, output_dir)
            self.assertEqual(results[0]["inspection_status"], "collected")
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError("Saved Panorama inventory must be reused")

            self.assertEqual(inspect_catalog(catalog, output_dir, resume=True), results)
            mock_get.assert_not_called()
            path = output_dir / results[0]["inventory_path"] / "manifest.json"
            manifest = json.loads(path.read_bytes())
            manifest["repository_accession"] = "999"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            before = (output_dir / "checkpoint.json").read_bytes()
            with self.assertRaisesRegex(ValueError, "Cannot resume inventory"):
                inspect_catalog(catalog, output_dir, resume=True)
            mock_get.assert_not_called()
            self.assertEqual((output_dir / "checkpoint.json").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
