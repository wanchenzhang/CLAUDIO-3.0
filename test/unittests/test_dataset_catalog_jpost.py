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
    RepositoryResponse, fetch_jpost_project, write_repository_outputs,
)


def metadata():
    return {
        "title": "Cross-linking dataset",
        "accession": [
            {"accession": "MS:1001919", "value": "PXD000001"},
            {"accession": "MS:1002632", "value": "JPST000001"},
        ],
        "datasetFiles": [
            {
                "accession": "PRIDE:0000404", "name": "Associated raw file URI",
                "value": "ftp://ftp.jpostdb.org/JPST000001/raw/sample.raw",
            },
            {
                "accession": "PRIDE:0000408", "name": "Search engine output file URI",
                "value": "ftp://ftp.jpostdb.org/JPST000001/results/sample%20%26%20%231.tsv",
            },
        ],
    }


def response(data):
    result = requests.Response()
    result.status_code = 200
    result.url = (
        "https://repository.jpostdb.org/proxi/datasets"
        "?accession=PXD000001&resultType=full"
    )
    result.headers["Content-Type"] = "application/json"
    result._content = (json.dumps(data) + "\n\n").encode("utf-8")
    return result


def snapshot(data):
    result = response(data)
    return RepositoryResponse(
        result.url, "2026-10-31T15:00:00+00:00", result.content,
        result.headers.copy(), data,
    )


class TestJpostInventory(unittest.TestCase):
    def candidate(self):
        return {
            "accession": "PXD000001", "repository": "jPOST", "title": "Test dataset",
            "review_status": "pending", "search_terms": ["cross-linking"],
        }

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_fetches_the_project_and_preserves_its_original_response(self, mock_get):
        data = metadata()
        raw_response = response(data)
        mock_get.return_value = raw_response

        project = fetch_jpost_project("PXD000001")

        self.assertEqual(project.data, data)
        self.assertEqual(project.raw_body, raw_response.content)
        self.assertEqual(project.request_url, raw_response.url)
        self.assertEqual(project.headers, raw_response.headers)
        mock_get.assert_called_once_with(
            "https://repository.jpostdb.org/proxi/datasets",
            params={"accession": "PXD000001", "resultType": "full"},
            headers={"Accept": "application/json"}, timeout=(15, 120),
        )

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_wrong_missing_or_conflicting_project_identifiers(self, mock_get):
        for problem in ("different_pxd", "missing_jpst", "conflicting_jpst"):
            with self.subTest(problem=problem):
                data = metadata()
                if problem == "different_pxd":
                    data["accession"][0]["value"] = "PXD000002"
                elif problem == "missing_jpst":
                    data["accession"].pop()
                else:
                    data["accession"].append({"accession": "MS:1002632", "value": "JPST000002"})
                mock_get.return_value = response(data)
                with self.assertRaises(ValueError):
                    fetch_jpost_project("PXD000001")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_invalid_accessions_without_requesting_data(self, mock_get):
        for accession in (None, "", "PXD1", "JPST000001", "PXD000001/other"):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError):
                    fetch_jpost_project(accession)
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_propagates_request_errors_and_rejects_invalid_json(self, mock_get):
        failed = response({})
        failed.status_code = 503
        mock_get.return_value = failed
        with self.assertRaises(requests.exceptions.HTTPError):
            fetch_jpost_project("PXD000001")
        mock_get.side_effect = requests.exceptions.Timeout("timed out")
        with self.assertRaises(requests.exceptions.Timeout):
            fetch_jpost_project("PXD000001")
        mock_get.side_effect = None
        failed.status_code = 200
        failed._content = b"<html>Service unavailable</html>"
        with self.assertRaises(ValueError):
            fetch_jpost_project("PXD000001")

    def test_writes_reported_file_links_and_keeps_unknown_counts_and_sizes_empty(self):
        data = metadata()
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            files = json.loads((output_dir / "files.json").read_bytes())

            self.assertEqual(manifest["source"], "jPOST")
            self.assertEqual(manifest["repository_accession"], "JPST000001")
            self.assertEqual(manifest["file_count"], 2)
            self.assertIsNone(manifest["reported_file_count"])
            self.assertEqual(manifest["file_inventory_scope"], "repository_reported_links")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_responses"], [])
            source = manifest["project_response"]
            self.assertEqual((output_dir / source["path"]).read_bytes(), snapshot(data).raw_body)
            self.assertEqual(source["sha256"], hashlib.sha256(snapshot(data).raw_body).hexdigest())
            self.assertEqual(files[1]["fileName"], "sample & #1.tsv")
            self.assertEqual(files[1]["repositoryMetadata"], data["datasetFiles"][1])
            self.assertEqual(files[1]["publicFileLocations"][0]["value"], data["datasetFiles"][1]["value"])
            self.assertIsNone(files[1]["fileSizeBytes"])
            self.assertNotEqual(files[0]["accession"], files[1]["accession"])
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["file_size_bytes"], "")
            self.assertEqual(rows[0]["checksum"], "")

    def test_accepts_an_explicitly_empty_link_list_but_rejects_a_missing_one(self):
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "empty"
            data = metadata()
            data["datasetFiles"] = []
            write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_count"], 0)
            self.assertIsNone(manifest["reported_file_count"])
            data.pop("datasetFiles")
            missing = Path(temp_dir) / "missing"
            with self.assertRaisesRegex(ValueError, "datasetFiles"):
                write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], missing)
            self.assertFalse(missing.exists())

    def test_rejects_unusable_file_links_before_creating_outputs(self):
        invalid = [
            "ftp://ftp.jpostdb.org/JPST000002/sample.raw",
            "ftp://example.org/JPST000001/sample.raw",
            "ftp://ftp.jpostdb.org/JPST000001/%2e%2e/sample.raw",
            "ftp://ftp.jpostdb.org/JPST000001/",
            None,
        ]
        for url in invalid:
            with self.subTest(url=url), TemporaryDirectory() as temp_dir:
                data = metadata()
                data["datasetFiles"][0]["value"] = url
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], output_dir)
                self.assertFalse(output_dir.exists())

    def test_rejects_duplicate_links_and_data_that_differs_from_its_source(self):
        for problem in ("duplicate", "edited"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                data = metadata()
                if problem == "duplicate":
                    data["datasetFiles"].append(deepcopy(data["datasetFiles"][0]))
                project = snapshot(data)
                if problem == "edited":
                    project.data["datasetFiles"][0]["value"] += ".changed"
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("jPOST", "PXD000001", project, [], output_dir)
                self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_collects_a_candidate_with_one_request_and_preserves_its_review_status(self, mock_get):
        mock_get.return_value = response(metadata())
        candidate = self.candidate()
        before = deepcopy(candidate)
        with TemporaryDirectory() as temp_dir:
            result = inspect_candidate(candidate, Path(temp_dir) / "project")
            self.assertEqual(result["inspection_status"], "collected")
            self.assertEqual(result["review_status"], "pending")
            self.assertEqual(result["file_count"], 2)
            self.assertEqual(candidate, before)
            self.assertEqual(mock_get.call_count, 1)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_records_missing_file_links_as_a_file_list_failure(self, mock_get):
        data = metadata()
        data.pop("datasetFiles")
        mock_get.return_value = response(data)
        with TemporaryDirectory() as temp_dir:
            result = inspect_candidate(self.candidate(), Path(temp_dir) / "project")
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "file_list")
            self.assertEqual(result["review_status"], "pending")
            self.assertIsNone(result["file_count"])

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_resumes_a_saved_project_and_rejects_a_changed_native_identity(self, mock_get):
        mock_get.return_value = response(metadata())
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            original = inspect_catalog(catalog, output_dir)
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError("Saved jPOST project must be reused")

            restored = inspect_catalog(catalog, output_dir, resume=True)

            self.assertEqual(restored, original)
            self.assertEqual(original[0]["attempt_number"], 1)
            mock_get.assert_not_called()
            path = output_dir / original[0]["inventory_path"] / "manifest.json"
            manifest = json.loads(path.read_bytes())
            manifest["repository_accession"] = "JPST000002"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            before = (output_dir / "checkpoint.json").read_bytes()
            with self.assertRaisesRegex(ValueError, "Cannot resume inventory"):
                inspect_catalog(catalog, output_dir, resume=True)
            mock_get.assert_not_called()
            self.assertEqual((output_dir / "checkpoint.json").read_bytes(), before)


    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_collects_and_resumes_a_list_wrapped_project_without_losing_the_raw_array(self, mock_get):
        data = [metadata()]
        raw_response = response(data)
        mock_get.return_value = raw_response
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"

            results = inspect_catalog(catalog, output_dir)

            self.assertEqual(results[0]["inspection_status"], "collected")
            self.assertEqual(results[0]["file_count"], 2)
            project_dir = output_dir / results[0]["inventory_path"]
            self.assertEqual(json.loads((project_dir / "project.json").read_bytes()), data)
            self.assertEqual((project_dir / "raw/project.json").read_bytes(), raw_response.content)
            manifest = json.loads((project_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["repository_accession"], "JPST000001")
            mock_get.assert_called_once()
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError("Saved array response must be reused")

            self.assertEqual(inspect_catalog(catalog, output_dir, resume=True), results)
            mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_selects_the_unique_matching_pxd_instead_of_the_first_project(self, mock_get):
        other = metadata()
        other["accession"][0]["value"] = "PXD000002"
        other["accession"][1]["value"] = "JPST000002"
        other["datasetFiles"] = []
        data = [other, metadata()]
        mock_get.return_value = response(data)
        project = fetch_jpost_project("PXD000001")

        self.assertEqual(project.data, data)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("jPOST", "PXD000001", project, [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["repository_accession"], "JPST000001")
            self.assertEqual(manifest["file_count"], 2)
            self.assertEqual(json.loads((output_dir / "project.json").read_bytes()), data)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_empty_unmatched_ambiguous_or_malformed_project_arrays(self, mock_get):
        other = metadata()
        other["accession"][0]["value"] = "PXD000002"
        for data in ([], [other], [metadata(), metadata()], [None], [metadata(), {}]):
            with self.subTest(data=data):
                mock_get.return_value = response(data)
                with self.assertRaises(ValueError):
                    fetch_jpost_project("PXD000001")


if __name__ == "__main__":
    unittest.main()
