import csv
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
    RepositoryResponse, fetch_iprox_project, write_repository_outputs,
)


def metadata():
    return {
        "accession": {"accession": "MS:1001919", "value": "PXD000001"},
        "title": "iProX cross-linking dataset",
        "datasetLink": {
            "accession": "MS:1002837", "name": "iProX dataset URI",
            "value": "https://www.iprox.org/page/project.html?id=IPX0000001000",
        },
        "dataFiles": [
            {
                "accession": "MS:1002846", "name": "Associated raw file URI",
                "value": "https://download.iprox.org/IPX0000001000/IPX0000001001/sample.raw",
            },
            {
                "accession": "MS:1002849", "name": "Search engine output file URI",
                "value": "http://download.iprox.org/IPX0000001000/IPX0000001002/sample%20results.tsv",
            },
        ],
    }


def response(data):
    result = requests.Response()
    result.status_code = 200
    result.url = "https://www.iprox.cn/proxi/datasets/PXD000001"
    result.headers["Content-Type"] = "application/json"
    result._content = (json.dumps(data) + "\n\n").encode("utf-8")
    return result


def snapshot(data):
    raw = response(data)
    return RepositoryResponse(
        raw.url, "2026-10-31T16:00:00+00:00", raw.content,
        raw.headers.copy(), data,
    )


class TestIproxInventory(unittest.TestCase):
    def candidate(self):
        return {
            "accession": "PXD000001", "repository": "iProX", "title": "Test dataset",
            "search_terms": ["crosslinking"], "review_status": "pending",
        }

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_fetches_object_shaped_identifiers_and_preserves_the_response(self, mock_get):
        data = metadata()
        raw = response(data)
        mock_get.return_value = raw

        project = fetch_iprox_project("PXD000001")

        self.assertEqual(project.data, data)
        self.assertEqual(project.raw_body, raw.content)
        self.assertEqual(project.headers, raw.headers)
        mock_get.assert_called_once_with(
            "https://www.iprox.cn/proxi/datasets/PXD000001",
            headers={"Accept": "application/json"}, timeout=(15, 120),
        )

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_accepts_list_shaped_fields_and_a_list_wrapped_project(self, mock_get):
        data = metadata()
        data["accession"] = [data["accession"]]
        data["datasetLink"] = [data["datasetLink"]]
        mock_get.return_value = response([data])
        project = fetch_iprox_project("PXD000001")

        self.assertEqual(project.data, [data])
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("iProX", "PXD000001", project, [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["repository_accession"], "IPX0000001000")
            self.assertEqual(manifest["file_count"], 2)
            self.assertEqual(json.loads((output_dir / "project.json").read_bytes()), [data])

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_wrong_missing_or_ambiguous_project_identity(self, mock_get):
        for problem in ("pxd", "missing_ipx", "invalid_ipx", "conflicting_ipx", "duplicate_project"):
            with self.subTest(problem=problem):
                data = metadata()
                if problem == "pxd":
                    data["accession"]["value"] = "PXD000002"
                elif problem == "missing_ipx":
                    data.pop("datasetLink")
                elif problem == "invalid_ipx":
                    data["datasetLink"]["value"] = "https://www.iprox.org/page/project.html?id=IPX123"
                elif problem == "conflicting_ipx":
                    other = deepcopy(data["datasetLink"])
                    other["value"] = "https://www.iprox.org/page/project.html?id=IPX0000002000"
                    data["datasetLink"] = [data["datasetLink"], other]
                else:
                    data = [data, deepcopy(data)]
                mock_get.return_value = response(data)
                with self.assertRaises(ValueError):
                    fetch_iprox_project("PXD000001")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_invalid_accessions_before_requesting_data(self, mock_get):
        for accession in (None, "", "PXD1", "IPX0000001000", "PXD000001/other"):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError):
                    fetch_iprox_project(accession)
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_propagates_http_and_timeout_errors_and_rejects_non_json_content(self, mock_get):
        raw = response({})
        raw.status_code = 503
        mock_get.return_value = raw
        with self.assertRaises(requests.exceptions.HTTPError):
            fetch_iprox_project("PXD000001")
        mock_get.side_effect = requests.exceptions.Timeout("timed out")
        with self.assertRaises(requests.exceptions.Timeout):
            fetch_iprox_project("PXD000001")
        mock_get.side_effect = None
        raw.status_code = 200
        raw._content = b"<html>Service unavailable</html>"
        with self.assertRaises(ValueError):
            fetch_iprox_project("PXD000001")

    def test_writes_reported_links_with_subproject_paths_and_unknown_sizes(self):
        data = metadata()
        project = snapshot(data)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("iProX", "PXD000001", project, [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            files = json.loads((output_dir / "files.json").read_bytes())

            self.assertEqual(manifest["source"], "iProX")
            self.assertEqual(manifest["repository_accession"], "IPX0000001000")
            self.assertEqual(manifest["file_count"], 2)
            self.assertIsNone(manifest["reported_file_count"])
            self.assertEqual(manifest["file_inventory_scope"], "repository_reported_links")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_responses"], [])
            self.assertEqual(files[0]["relativePath"], "IPX0000001001/sample.raw")
            self.assertEqual(files[1]["fileName"], "sample results.tsv")
            self.assertEqual(files[1]["repositoryMetadata"], data["dataFiles"][1])
            self.assertEqual(files[1]["publicFileLocations"][0]["value"], data["dataFiles"][1]["value"])
            self.assertIsNone(files[0]["fileSizeBytes"])
            source = manifest["project_response"]
            self.assertEqual((output_dir / source["path"]).read_bytes(), project.raw_body)
            self.assertEqual(source["sha256"], hashlib.sha256(project.raw_body).hexdigest())
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["file_size_bytes"], "")
            self.assertEqual(rows[0]["checksum"], "")

    def test_keeps_explicitly_empty_file_lists_distinct_from_missing_metadata(self):
        with TemporaryDirectory() as temp_dir:
            data = metadata()
            data["dataFiles"] = []
            output_dir = Path(temp_dir) / "empty"
            write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_count"], 0)
            self.assertIsNone(manifest["reported_file_count"])
            data.pop("dataFiles")
            output_dir = Path(temp_dir) / "missing"
            with self.assertRaisesRegex(ValueError, "dataFiles"):
                write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
            self.assertFalse(output_dir.exists())

    def test_accepts_a_single_file_object_without_changing_the_original_metadata(self):
        data = metadata()
        data["dataFiles"] = data["dataFiles"][0]
        before = deepcopy(data)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_count"], 1)
            self.assertEqual(json.loads((output_dir / "project.json").read_bytes()), before)
            self.assertEqual(data, before)

    def test_rejects_wrong_roots_invalid_links_and_duplicates_before_writing(self):
        for problem in ("root", "host", "traversal", "duplicate"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                data = metadata()
                if problem == "root":
                    data["dataFiles"][0]["value"] = "https://download.iprox.org/IPX0000002000/sample.raw"
                elif problem == "host":
                    data["dataFiles"][0]["value"] = "https://example.org/IPX0000001000/sample.raw"
                elif problem == "traversal":
                    data["dataFiles"][0]["value"] = "https://download.iprox.org/IPX0000001000/%2e%2e/sample.raw"
                else:
                    data["dataFiles"].append(deepcopy(data["dataFiles"][0]))
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
                self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_collects_a_candidate_once_and_preserves_review_status(self, mock_get):
        mock_get.return_value = response(metadata())
        candidate = self.candidate()
        before = deepcopy(candidate)
        with TemporaryDirectory() as temp_dir:
            result = inspect_candidate(candidate, Path(temp_dir) / "project")
            self.assertEqual(result["inspection_status"], "collected")
            self.assertEqual(result["file_count"], 2)
            self.assertEqual(result["review_status"], "pending")
            self.assertEqual(candidate, before)
            self.assertEqual(mock_get.call_count, 1)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_records_missing_file_metadata_as_a_file_list_failure(self, mock_get):
        data = metadata()
        data.pop("dataFiles")
        mock_get.return_value = response(data)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            result = inspect_candidate(self.candidate(), output_dir)
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "file_list")
            self.assertIsNone(result["file_count"])
            self.assertEqual(result["review_status"], "pending")
            self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_resumes_saved_outputs_and_rejects_a_changed_native_identity(self, mock_get):
        mock_get.return_value = response(metadata())
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            results = inspect_catalog(catalog, output_dir)
            self.assertEqual(results[0]["inspection_status"], "collected")
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError("Saved iProX inventory must be reused")

            self.assertEqual(inspect_catalog(catalog, output_dir, resume=True), results)
            mock_get.assert_not_called()
            manifest_path = output_dir / results[0]["inventory_path"] / "manifest.json"
            manifest = json.loads(manifest_path.read_bytes())
            manifest["repository_accession"] = "IPX0000002000"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            before = (output_dir / "checkpoint.json").read_bytes()
            with self.assertRaisesRegex(ValueError, "Cannot resume inventory"):
                inspect_catalog(catalog, output_dir, resume=True)
            mock_get.assert_not_called()
            self.assertEqual((output_dir / "checkpoint.json").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
