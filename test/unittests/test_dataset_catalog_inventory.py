import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import csv
import hashlib
import requests

from claudio.dataset_catalog.inventory import inspect_candidate, inspect_catalog


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
                "claudio.dataset_catalog.pride.csv.DictWriter.writeheader",
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
    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_batch_retains_failed_unsupported_and_collected_candidates(self, mock_get):
        base_url = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects"
        candidates = [
            self._candidate(),
            dict(self._candidate("MassIVE"), accession="PXD000002"),
            dict(self._candidate(), accession="PXD000003"),
        ]

        def respond(url, **kwargs):
            if url == f"{base_url}/PXD000001":
                raise requests.exceptions.Timeout("first project timed out")
            if url == f"{base_url}/PXD000003":
                return make_response(url, {
                    "accession": "PXD000003", "title": "Third dataset",
                })
            if url == f"{base_url}/PXD000003/files":
                return make_response(
                    url,
                    [{
                        "accession": "file-3",
                        "projectAccessions": ["PXD000003"],
                        "fileName": "sample.raw",
                    }],
                    params=kwargs.get("params"),
                    headers={"total_records": "1"},
                )
            raise AssertionError(f"Unexpected request: {url}")

        mock_get.side_effect = respond

        with TemporaryDirectory() as temp_dir:
            catalog_path = Path(temp_dir) / "candidates.json"
            catalog_body = (json.dumps(candidates) + "\n").encode("utf-8")
            catalog_path.write_bytes(catalog_body)
            output_dir = Path(temp_dir) / "inventory"

            results = inspect_catalog(catalog_path, output_dir)

            self.assertEqual(
                [row["inspection_status"] for row in results],
                ["failed", "adapter_pending", "collected"],
            )
            self.assertEqual(
                [row["accession"] for row in results],
                ["PXD000001", "PXD000002", "PXD000003"],
            )
            self.assertTrue(
                all(row["review_status"] == "pending" for row in results)
            )
            self.assertEqual(
                (output_dir / "input_candidates.json").read_bytes(),
                catalog_body,
            )
            checkpoint = json.loads(
                (output_dir / "checkpoint.json").read_text(encoding="utf-8")
            )
            self.assertEqual(checkpoint["status"], "finished")
            self.assertEqual(checkpoint["results"], results)
            self.assertEqual(
                json.loads(
                    (output_dir / "projects.json").read_text(encoding="utf-8")
                ),
                results,
            )
            with (output_dir / "projects.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)
            self.assertEqual(rows[0]["file_count"], "")
            self.assertEqual(rows[2]["file_count"], "1")

            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "finished")
            self.assertEqual(manifest["candidate_count"], 3)
            self.assertEqual(
                manifest["status_counts"],
                {"collected": 1, "failed": 1, "adapter_pending": 1},
            )
            self.assertEqual(
                manifest["source_catalog"]["sha256"],
                hashlib.sha256(catalog_body).hexdigest(),
            )
            relative_path = Path(results[2]["inventory_path"])
            self.assertFalse(relative_path.is_absolute())
            self.assertTrue(
                (output_dir / relative_path / "manifest.json").is_file()
            )

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_batch_saves_progress_before_starting_the_next_candidate(self, mock_get):
        files = [{
            "accession": "file-1",
            "projectAccessions": ["PXD000001"],
            "fileName": "sample.raw",
        }]
        mock_get.side_effect = [
            *self._pride_responses(files),
            KeyboardInterrupt(),
        ]
        candidates = [
            self._candidate(),
            dict(self._candidate(), accession="PXD000002"),
        ]

        with TemporaryDirectory() as temp_dir:
            catalog_path = Path(temp_dir) / "candidates.json"
            catalog_path.write_text(json.dumps(candidates), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"

            with self.assertRaises(KeyboardInterrupt):
                inspect_catalog(catalog_path, output_dir)

            checkpoint = json.loads(
                (output_dir / "checkpoint.json").read_text(encoding="utf-8")
            )
            self.assertEqual(checkpoint["status"], "running")
            self.assertEqual(checkpoint["candidate_count"], 2)
            self.assertEqual(len(checkpoint["results"]), 1)
            self.assertEqual(
                checkpoint["results"][0]["inspection_status"], "collected"
            )
            self.assertFalse((output_dir / "manifest.json").exists())
    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_resume_reuses_completed_inventory_and_continues_after_interrupt(self, mock_get):
        first_files = [{
            "accession": "file-1",
            "projectAccessions": ["PXD000001"],
            "fileName": "first.raw",
        }]
        mock_get.side_effect = [
            *self._pride_responses(first_files),
            KeyboardInterrupt(),
        ]
        candidates = [
            self._candidate(),
            dict(self._candidate(), accession="PXD000002"),
        ]

        with TemporaryDirectory() as temp_dir:
            catalog_path = Path(temp_dir) / "candidates.json"
            catalog_path.write_text(json.dumps(candidates), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"

            with self.assertRaises(KeyboardInterrupt):
                inspect_catalog(catalog_path, output_dir)

            checkpoint = json.loads(
                (output_dir / "checkpoint.json").read_text(encoding="utf-8")
            )
            first_result = checkpoint["results"][0]
            saved_manifest = (
                output_dir / first_result["inventory_path"] / "manifest.json"
            )
            original_manifest = saved_manifest.read_bytes()

            second_url = (
                "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000002"
            )

            def respond(url, **kwargs):
                if url == second_url:
                    return make_response(url, {
                        "accession": "PXD000002",
                        "title": "Second dataset",
                    })
                if url == f"{second_url}/files":
                    return make_response(
                        url,
                        [{
                            "accession": "file-2",
                            "projectAccessions": ["PXD000002"],
                            "fileName": "second.raw",
                        }],
                        params=kwargs.get("params"),
                        headers={"total_records": "1"},
                    )
                raise AssertionError(f"Unexpected request during resume: {url}")

            mock_get.reset_mock()
            mock_get.side_effect = respond

            results = inspect_catalog(
                catalog_path, output_dir, resume=True
            )

            self.assertEqual(
                [row["inspection_status"] for row in results],
                ["collected", "collected"],
            )
            self.assertEqual(results[0], first_result)
            self.assertEqual(saved_manifest.read_bytes(), original_manifest)
            self.assertEqual(
                [call.args[0] for call in mock_get.call_args_list],
                [second_url, f"{second_url}/files"],
            )

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_resume_rejects_a_changed_catalog_without_modifying_the_run(self, mock_get):
        files = [{
            "accession": "file-1",
            "projectAccessions": ["PXD000001"],
            "fileName": "sample.raw",
        }]
        mock_get.side_effect = self._pride_responses(files)

        with TemporaryDirectory() as temp_dir:
            catalog_path = Path(temp_dir) / "candidates.json"
            catalog_path.write_text(
                json.dumps([self._candidate()]), encoding="utf-8"
            )
            output_dir = Path(temp_dir) / "inventory"
            inspect_catalog(catalog_path, output_dir)

            checkpoint_path = output_dir / "checkpoint.json"
            manifest_path = output_dir / "manifest.json"
            input_path = output_dir / "input_candidates.json"
            original_checkpoint = checkpoint_path.read_bytes()
            original_manifest = manifest_path.read_bytes()
            original_input = input_path.read_bytes()

            changed = self._candidate()
            changed["title"] = "Changed catalog metadata"
            catalog_path.write_text(json.dumps([changed]), encoding="utf-8")
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError(
                "Changed catalog must not trigger network requests"
            )

            with self.assertRaisesRegex(ValueError, "catalog"):
                inspect_catalog(catalog_path, output_dir, resume=True)

            mock_get.assert_not_called()
            self.assertEqual(checkpoint_path.read_bytes(), original_checkpoint)
            self.assertEqual(manifest_path.read_bytes(), original_manifest)
            self.assertEqual(input_path.read_bytes(), original_input)
    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_resume_rejects_damaged_saved_files_before_making_requests(self, mock_get):
        files = [{
            "accession": "file-1",
            "projectAccessions": ["PXD000001"],
            "fileName": "sample.raw",
        }]
        for damaged in ("input_candidates.json", "raw", "files.json", "files.csv"):
            with self.subTest(damaged=damaged), TemporaryDirectory() as temp_dir:
                catalog = Path(temp_dir) / "candidates.json"
                catalog.write_text(json.dumps([self._candidate()]), encoding="utf-8")
                output_dir = Path(temp_dir) / "inventory"
                mock_get.side_effect = self._pride_responses(files)
                result = inspect_catalog(catalog, output_dir)[0]
                project_dir = output_dir / result["inventory_path"]
                if damaged == "input_candidates.json":
                    target = output_dir / damaged
                elif damaged == "raw":
                    target = project_dir / "raw/files_page_0000.json"
                else:
                    target = project_dir / damaged
                target.write_text("changed data\n", encoding="utf-8")
                before = {
                    p.relative_to(output_dir): p.read_bytes()
                    for p in output_dir.rglob("*") if p.is_file()
                }
                mock_get.reset_mock()
                mock_get.side_effect = AssertionError("Resume must validate saved files first")

                with self.assertRaisesRegex(ValueError, "Cannot resume inventory"):
                    inspect_catalog(catalog, output_dir, resume=True)

                mock_get.assert_not_called()
                self.assertEqual({
                    p.relative_to(output_dir): p.read_bytes()
                    for p in output_dir.rglob("*") if p.is_file()
                }, before)

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_resume_keeps_partial_attempt_and_writes_to_a_new_directory(self, mock_get):
        mock_get.side_effect = KeyboardInterrupt()
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self._candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            with self.assertRaises(KeyboardInterrupt):
                inspect_catalog(catalog, output_dir)

            partial = output_dir / "projects/candidate_000001/attempt_0001"
            partial.mkdir(parents=True)
            old_file = partial / "files.json"
            old_file.write_bytes(b"unfinished previous attempt\n")
            mock_get.reset_mock()
            mock_get.side_effect = self._pride_responses([])

            results = inspect_catalog(catalog, output_dir, resume=True)

            self.assertEqual(results[0]["inspection_status"], "collected")
            self.assertEqual(results[0]["file_count"], 0)
            self.assertEqual(results[0]["attempt_number"], 2)
            self.assertEqual(
                results[0]["inventory_path"],
                "projects/candidate_000001/attempt_0002",
            )
            self.assertEqual(old_file.read_bytes(), b"unfinished previous attempt\n")
            self.assertTrue((
                output_dir / results[0]["inventory_path"] / "manifest.json"
            ).is_file())

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_resume_keeps_recorded_failures_and_unsupported_candidates(self, mock_get):
        candidates = [
            self._candidate(),
            dict(self._candidate("MassIVE"), accession="PXD000002"),
            dict(self._candidate(), accession="PXD000003"),
        ]
        mock_get.side_effect = [requests.exceptions.Timeout("timed out"), KeyboardInterrupt()]
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps(candidates), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            with self.assertRaises(KeyboardInterrupt):
                inspect_catalog(catalog, output_dir)
            previous = json.loads((output_dir / "checkpoint.json").read_bytes())["results"]

            base_url = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000003"
            mock_get.reset_mock()
            mock_get.side_effect = [
                make_response(base_url, {"accession": "PXD000003"}),
                make_response(base_url + "/files", [], headers={"total_records": "0"}),
            ]
            results = inspect_catalog(catalog, output_dir, resume=True)

            self.assertEqual(results[:2], previous)
            self.assertEqual(
                [r["inspection_status"] for r in results],
                ["failed", "adapter_pending", "collected"],
            )
            self.assertEqual(
                [call.args[0] for call in mock_get.call_args_list],
                [base_url, base_url + "/files"],
            )

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_resuming_a_finished_inventory_preserves_existing_reports(self, mock_get):
        mock_get.side_effect = self._pride_responses([])
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self._candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            original = inspect_catalog(catalog, output_dir)
            before = {
                p.relative_to(output_dir): p.read_bytes()
                for p in output_dir.rglob("*") if p.is_file()
            }
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError("Completed inventory must be reused")

            results = inspect_catalog(catalog, output_dir, resume=True)

            self.assertEqual(results, original)
            mock_get.assert_not_called()
            self.assertEqual({
                p.relative_to(output_dir): p.read_bytes()
                for p in output_dir.rglob("*") if p.is_file()
            }, before)


if __name__ == "__main__":
    unittest.main()
