import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import requests

from claudio.dataset_catalog.inventory import inspect_candidate


def make_response(url, data, params=None, headers=None):
    response = requests.Response()
    response.status_code = 200
    response.url = requests.Request("GET", url, params=params).prepare().url
    response.headers["Content-Type"] = "application/json"
    response.headers.update(headers or {})
    response._content = json.dumps(data).encode("utf-8")
    return response


class TestDatasetCatalogInventory(unittest.TestCase):
    def _candidate(self, repository="PRIDE"):
        return {
            "accession": "PXD000001",
            "repository": repository,
            "title": "Test dataset",
            "search_terms": ["cross-linking"],
            "review_status": "pending",
        }

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_collects_pride_inventory_and_preserves_candidate_provenance(self, mock_get):
        base_url = (
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000001"
        )
        project = {
            "accession": "PXD000001",
            "title": "Test dataset",
            "references": [],
        }
        files = [{
            "accession": "file-1",
            "projectAccessions": ["PXD000001"],
            "fileName": "sample.raw",
        }]

        def respond(url, **kwargs):
            if url == base_url:
                return make_response(url, project)
            if url == f"{base_url}/files":
                return make_response(
                    url, files, params=kwargs.get("params"),
                    headers={"total_records": "1"},
                )
            raise AssertionError(f"Unexpected request: {url}")

        mock_get.side_effect = respond
        candidate = self._candidate()
        original = deepcopy(candidate)

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "PXD000001"
            result = inspect_candidate(candidate, output_dir)

            self.assertIsNotNone(result)
            self.assertEqual(result["inspection_status"], "collected")
            self.assertEqual(result["accession"], "PXD000001")
            self.assertEqual(result["review_status"], "pending")
            self.assertEqual(result["search_terms"], ["cross-linking"])
            self.assertEqual(result["file_count"], 1)
            self.assertEqual(result["error_stage"], "")
            self.assertTrue((output_dir / "manifest.json").is_file())
            self.assertEqual(
                json.loads(
                    (output_dir / "project.json").read_text(encoding="utf-8")
                ),
                project,
            )
            self.assertEqual(candidate, original)

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_retains_candidates_from_repositories_without_an_adapter(self, mock_get):
        candidate = self._candidate(repository="MassIVE")

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "PXD000001"
            result = inspect_candidate(candidate, output_dir)

            self.assertIsNotNone(result)
            self.assertEqual(result["inspection_status"], "adapter_pending")
            self.assertEqual(result["repository"], "MassIVE")
            self.assertEqual(result["review_status"], "pending")
            self.assertIsNone(result["file_count"])
            self.assertFalse(output_dir.exists())

        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_records_request_failure_without_excluding_the_candidate(self, mock_get):
        mock_get.side_effect = requests.exceptions.Timeout("read timed out")

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "PXD000001"
            result = inspect_candidate(self._candidate(), output_dir)

            self.assertIsNotNone(result)
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["review_status"], "pending")
            self.assertIsNone(result["file_count"])
            self.assertEqual(result["error_stage"], "project_metadata")
            self.assertEqual(result["error_type"], "Timeout")
            self.assertIn("read timed out", result["error_message"])
            self.assertFalse((output_dir / "manifest.json").exists())
    def _pride_responses(self, files):
        base_url = (
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000001"
        )
        return [
            make_response(base_url, {
                "accession": "PXD000001",
                "title": "Test dataset",
                "submissionType": "PARTIAL",
            }),
            make_response(
                f"{base_url}/files",
                files,
                params={"pageSize": 100, "page": 0},
                headers={"total_records": str(len(files))},
            ),
        ]

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_records_file_list_failure_at_the_correct_stage(self, mock_get):
        mock_get.side_effect = [
            self._pride_responses([])[0],
            requests.exceptions.Timeout("file list timed out"),
        ]

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "PXD000001"
            result = inspect_candidate(self._candidate(), output_dir)

            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "file_list")
            self.assertEqual(result["error_type"], "Timeout")
            self.assertIn("file list timed out", result["error_message"])
            self.assertEqual(result["review_status"], "pending")
            self.assertIsNone(result["file_count"])
            self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_records_write_failure_without_marking_partial_output_collected(self, mock_get):
        mock_get.side_effect = self._pride_responses([{
            "accession": "file-1",
            "projectAccessions": ["PXD000001"],
            "fileName": "sample.raw",
        }])

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "PXD000001"
            with patch(
                "claudio.dataset_catalog.pride_outputs.csv.DictWriter.writeheader",
                side_effect=OSError("No space left on device"),
            ):
                result = inspect_candidate(self._candidate(), output_dir)

            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "output")
            self.assertEqual(result["error_type"], "OSError")
            self.assertIn("No space left on device", result["error_message"])
            self.assertEqual(result["review_status"], "pending")
            self.assertIsNone(result["file_count"])
            self.assertTrue((output_dir / "files.json").is_file())
            self.assertFalse((output_dir / "manifest.json").exists())

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_preserves_confirmed_zero_files_as_a_collected_inventory(self, mock_get):
        mock_get.side_effect = self._pride_responses([])

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "PXD000001"
            result = inspect_candidate(self._candidate(), output_dir)

            self.assertEqual(result["inspection_status"], "collected")
            self.assertEqual(result["file_count"], 0)
            self.assertEqual(result["review_status"], "pending")
            self.assertEqual(result["error_stage"], "")
            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["file_count"], 0)

if __name__ == "__main__":
    unittest.main()