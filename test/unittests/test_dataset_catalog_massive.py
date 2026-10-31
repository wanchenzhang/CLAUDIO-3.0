import json
import csv
import hashlib
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from urllib.parse import parse_qs, urlsplit

import requests

from claudio.dataset_catalog.repositories import (
    RepositoryResponse, fetch_massive_project, fetch_massive_file_pages,
    write_repository_outputs,
)
from claudio.dataset_catalog.inventory import inspect_candidate, inspect_catalog


def massive_metadata():
    return {
        "title": "Cross-linking dataset",
        "accession": [
            {"accession": "MS:1001919", "value": "PXD000001"},
            {"accession": "MS:1002487", "value": "MSV000012345"},
        ],
    }


def massive_file(path="raw/sample.raw", **fields):
    return {
        "file_descriptor": "f.MSV000012345/" + path,
        "name": path.rsplit("/", 1)[-1],
        "collection": path.split("/", 1)[0],
        "size": 123,
        "type": "File",
        **fields,
    }


def http_response(data=None, files=None, total=None):
    response = requests.Response()
    response.status_code = 200
    if files is None:
        response.url = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/PXD000001"
        response.headers["Content-Type"] = "application/json"
        response._content = (json.dumps(data) + "\n").encode("utf-8")
    else:
        response.url = "https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?accession=MSV000012345"
        response.headers["Content-Type"] = "text/html; charset=UTF-8"
        payload = {
            "source": "ProteoSAFe",
            "total_rows": len(files) if total is None else total,
            "row_data": files,
        }
        response._content = (
            "<html><script>\nvar dataset_files = "
            + json.dumps(payload) + ";\n</script></html>\n"
        ).encode("utf-8")
    return response


def snapshot(response, data):
    return RepositoryResponse(
        response.url, "2026-10-31T15:00:00+00:00", response.content,
        response.headers.copy(), data,
    )


class TestMassiveProject(unittest.TestCase):
    def setUp(self):
        self.url = (
            "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/PXD000001"
        )
        self.metadata = {
            "title": "Cross-linking dataset",
            "accession": [
                {"accession": "MS:1001919", "value": "PXD000001"},
                {"accession": "MS:1002487", "value": "MSV000012345"},
            ],
            "datasetLink": [{
                "name": "Dataset FTP location",
                "value": "ftp://massive-ftp.ucsd.edu/v04/MSV000012345/",
            }],
            "species": [[{"name": "taxonomy: scientific name", "value": "Homo sapiens"}]],
        }

    def response(self, data):
        response = requests.Response()
        response.status_code = 200
        response.url = self.url
        response.headers["Content-Type"] = "application/json"
        response._content = (json.dumps(data) + "\n\n").encode("utf-8")
        return response

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_preserves_metadata_and_the_original_response(self, mock_get):
        response = self.response(self.metadata)
        mock_get.return_value = response

        project = fetch_massive_project("PXD000001")

        self.assertIsInstance(project, RepositoryResponse)
        self.assertEqual(project.data, self.metadata)
        self.assertEqual(project.raw_body, response.content)
        self.assertEqual(project.request_url, self.url)
        self.assertEqual(project.headers, response.headers)
        self.assertTrue(project.retrieved_at_utc)
        mock_get.assert_called_once_with(self.url, timeout=(15, 120))

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_a_response_for_another_pxd(self, mock_get):
        metadata = deepcopy(self.metadata)
        metadata["accession"][0]["value"] = "PXD000002"
        mock_get.return_value = self.response(metadata)

        with self.assertRaisesRegex(ValueError, "accession mismatch"):
            fetch_massive_project("PXD000001")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_missing_conflicting_or_invalid_identifiers(self, mock_get):
        for entries in (
            [],
            [self.metadata["accession"][0]],
            [self.metadata["accession"][1]],
            self.metadata["accession"] + [
                {"accession": "MS:1002487", "value": "MSV000054321"}
            ],
            [self.metadata["accession"][0],
             {"accession": "MS:1002487", "value": "MSV123"}],
        ):
            with self.subTest(entries=entries):
                mock_get.return_value = self.response({"accession": entries})
                with self.assertRaisesRegex(ValueError, "identifier"):
                    fetch_massive_project("PXD000001")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_invalid_json_or_response_shapes(self, mock_get):
        for data in ([], None, "not metadata", {"accession": None}):
            with self.subTest(data=data):
                mock_get.return_value = self.response(data)
                with self.assertRaises(ValueError):
                    fetch_massive_project("PXD000001")
        response = self.response({})
        response._content = b"<html>Service unavailable</html>"
        mock_get.return_value = response
        with self.assertRaises(ValueError):
            fetch_massive_project("PXD000001")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_invalid_accessions_before_requesting_data(self, mock_get):
        for accession in (None, "", "MSV000012345", "PXD1", "PXD000001/other"):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError):
                    fetch_massive_project(accession)
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_propagates_request_failures(self, mock_get):
        response = self.response({})
        response.status_code = 503
        mock_get.return_value = response
        with self.assertRaises(requests.exceptions.HTTPError):
            fetch_massive_project("PXD000001")

        mock_get.side_effect = requests.exceptions.Timeout("read timed out")
        with self.assertRaises(requests.exceptions.Timeout):
            fetch_massive_project("PXD000001")


class TestMassiveFiles(unittest.TestCase):
    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_fetches_the_complete_main_site_listing_and_preserves_its_html(self, mock_get):
        records = [massive_file("raw/a/sample.raw"), massive_file("raw/b/sample.raw")]
        response = http_response(files=records)
        mock_get.return_value = response

        pages = fetch_massive_file_pages("MSV000012345")

        self.assertEqual(len(pages), 1)
        self.assertEqual(pages[0].raw_body, response.content)
        self.assertEqual(pages[0].data["row_data"], records)
        self.assertEqual(pages[0].data["total_rows"], 2)
        mock_get.assert_called_once_with(
            "https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp",
            params={"accession": "MSV000012345"}, timeout=(15, 120),
        )

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_accepts_a_confirmed_empty_listing(self, mock_get):
        mock_get.return_value = http_response(files=[])
        pages = fetch_massive_file_pages("MSV000012345")
        self.assertEqual(pages[0].data["row_data"], [])
        self.assertEqual(pages[0].data["total_rows"], 0)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_wrong_repeated_or_incomplete_file_records(self, mock_get):
        cases = [
            ([massive_file()], 2),
            ([massive_file()] * 2, 2),
            ([massive_file(file_descriptor="f.MSV000054321/raw/sample.raw")], 1),
            ([massive_file(file_descriptor="f.MSV000012345/../sample.raw")], 1),
            ([massive_file(type="Directory")], 1),
            ([{"name": "sample.raw"}], 1),
            ([massive_file(size=-1)], 1),
        ]
        for records, total in cases:
            with self.subTest(records=records, total=total):
                mock_get.return_value = http_response(files=records, total=total)
                with self.assertRaises(ValueError):
                    fetch_massive_file_pages("MSV000012345")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_missing_or_truncated_embedded_json(self, mock_get):
        for body in (
            b"<html>Service unavailable</html>",
            b'<script>var dataset_files = {"total_rows": 3, "row_data": [',
            b"<script>var dataset_files = null;</script>",
        ):
            with self.subTest(body=body):
                response = http_response(files=[])
                response._content = body
                mock_get.return_value = response
                with self.assertRaises(ValueError):
                    fetch_massive_file_pages("MSV000012345")

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_validates_the_native_accession_before_requesting_files(self, mock_get):
        for accession in (None, "", "PXD000001", "MSV123", "MSV000012345/other"):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError):
                    fetch_massive_file_pages(accession)
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_propagates_file_request_errors(self, mock_get):
        response = http_response(files=[])
        response.status_code = 503
        mock_get.return_value = response
        with self.assertRaises(requests.exceptions.HTTPError):
            fetch_massive_file_pages("MSV000012345")
        mock_get.side_effect = requests.exceptions.Timeout("file list timed out")
        with self.assertRaises(requests.exceptions.Timeout):
            fetch_massive_file_pages("MSV000012345")


class TestMassiveInventory(unittest.TestCase):
    def candidate(self):
        return {
            "accession": "PXD000001", "repository": "MassIVE",
            "title": "Cross-linking dataset", "search_terms": ["XL-MS"],
            "review_status": "pending",
        }

    def inputs(self, records=None):
        metadata = massive_metadata()
        records = [massive_file()] if records is None else records
        project = snapshot(http_response(metadata), metadata)
        data = {"source": "ProteoSAFe", "total_rows": len(records), "row_data": records}
        page = snapshot(http_response(files=records), data)
        return project, [page]

    def test_writes_common_outputs_with_native_identity_and_original_html(self):
        records = [massive_file("raw/a/sample.raw"), massive_file("result/b/sample & #1.tsv")]
        project, pages = self.inputs(records)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("MassIVE", "PXD000001", project, pages, output_dir)

            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            files = json.loads((output_dir / "files.json").read_bytes())
            self.assertEqual(manifest["source"], "MassIVE")
            self.assertEqual(manifest["repository_accession"], "MSV000012345")
            self.assertEqual(manifest["project_accession"], "PXD000001")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_count"], 2)
            self.assertEqual([f["accession"] for f in files], [r["file_descriptor"] for r in records])
            self.assertEqual(files[1]["projectAccessions"], ["PXD000001"])
            self.assertEqual(files[1]["relativePath"], "result/b/sample & #1.tsv")
            location = files[1]["publicFileLocations"][0]["value"]
            self.assertEqual(parse_qs(urlsplit(location).query)["file"], [records[1]["file_descriptor"]])
            saved = manifest["file_responses"][0]
            self.assertTrue(saved["path"].endswith(".html"))
            self.assertEqual((output_dir / saved["path"]).read_bytes(), pages[0].raw_body)
            self.assertEqual(saved["sha256"], hashlib.sha256(pages[0].raw_body).hexdigest())
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[1]["file_name"], "sample & #1.tsv")

    def test_writes_a_header_for_a_confirmed_empty_inventory(self):
        project, pages = self.inputs([])
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("MassIVE", "PXD000001", project, pages, output_dir)
            self.assertEqual(json.loads((output_dir / "files.json").read_bytes()), [])
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                self.assertIn("file_accession", reader.fieldnames)
                self.assertEqual(list(reader), [])

    def test_rejects_bad_sources_before_creating_outputs(self):
        for problem in ("project", "missing_files", "raw_html", "different_msv"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                project, pages = self.inputs()
                if problem == "project":
                    project.data["accession"][0]["value"] = "PXD000002"
                elif problem == "missing_files":
                    pages = []
                elif problem == "raw_html":
                    pages[0].data["row_data"][0]["size"] = 999
                else:
                    metadata = deepcopy(project.data)
                    metadata["accession"][1]["value"] = "MSV000054321"
                    project = snapshot(http_response(metadata), metadata)
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("MassIVE", "PXD000001", project, pages, output_dir)
                self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_collects_massive_candidates_without_changing_review_status(self, mock_get):
        mock_get.side_effect = [http_response(massive_metadata()), http_response(files=[massive_file()])]
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            result = inspect_candidate(self.candidate(), output_dir)
            self.assertEqual(result["inspection_status"], "collected")
            self.assertEqual(result["repository"], "MassIVE")
            self.assertEqual(result["review_status"], "pending")
            self.assertEqual(result["file_count"], 1)
            self.assertEqual(mock_get.call_count, 2)
            self.assertTrue((output_dir / "manifest.json").is_file())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_records_massive_file_request_failure_without_excluding_the_candidate(self, mock_get):
        mock_get.side_effect = [http_response(massive_metadata()), requests.exceptions.Timeout("timed out")]
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            result = inspect_candidate(self.candidate(), output_dir)
            self.assertEqual(result["inspection_status"], "failed")
            self.assertEqual(result["error_stage"], "file_list")
            self.assertEqual(result["review_status"], "pending")
            self.assertIsNone(result["file_count"])
            self.assertFalse(output_dir.exists())

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_resume_checks_saved_massive_outputs_without_refetching_them(self, mock_get):
        mock_get.side_effect = [http_response(massive_metadata()), http_response(files=[massive_file()])]
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            results = inspect_catalog(catalog, output_dir)
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError("Completed MassIVE project must be reused")

            restored = inspect_catalog(catalog, output_dir, resume=True)

            self.assertEqual(restored, results)
            self.assertEqual(results[0]["attempt_number"], 1)
            mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_resume_continues_after_a_massive_project_is_collected(self, mock_get):
        candidates = [
            self.candidate(),
            dict(self.candidate(), accession="PXD000002", repository="PRIDE"),
        ]
        mock_get.side_effect = [
            http_response(massive_metadata()), http_response(files=[massive_file()]),
            KeyboardInterrupt(),
        ]
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps(candidates), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            with self.assertRaises(KeyboardInterrupt):
                inspect_catalog(catalog, output_dir)
            completed = json.loads((output_dir / "checkpoint.json").read_bytes())["results"][0]

            base_url = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000002"
            project = http_response({"accession": "PXD000002", "title": "PRIDE dataset"})
            project.url = base_url
            files = http_response([])
            files.url = base_url + "/files?pageSize=100&page=0"
            files.headers["total_records"] = "0"
            mock_get.reset_mock()
            mock_get.side_effect = [project, files]

            restored = inspect_catalog(catalog, output_dir, resume=True)

            self.assertEqual(restored[0], completed)
            self.assertEqual([row["inspection_status"] for row in restored], ["collected", "collected"])
            self.assertEqual(restored[1]["file_count"], 0)
            self.assertEqual(
                [call.args[0] for call in mock_get.call_args_list],
                [base_url, base_url + "/files"],
            )

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_resume_rejects_a_damaged_massive_response_without_changing_reports(self, mock_get):
        mock_get.side_effect = [http_response(massive_metadata()), http_response(files=[massive_file()])]
        with TemporaryDirectory() as temp_dir:
            catalog = Path(temp_dir) / "candidates.json"
            catalog.write_text(json.dumps([self.candidate()]), encoding="utf-8")
            output_dir = Path(temp_dir) / "inventory"
            result = inspect_catalog(catalog, output_dir)[0]
            raw = output_dir / result["inventory_path"] / "raw/files_page_0000.html"
            raw.write_bytes(b"changed saved page")
            before = (output_dir / "checkpoint.json").read_bytes()
            mock_get.reset_mock()
            mock_get.side_effect = AssertionError("Damaged response must be rejected first")

            with self.assertRaisesRegex(ValueError, "Cannot resume inventory"):
                inspect_catalog(catalog, output_dir, resume=True)

            mock_get.assert_not_called()
            self.assertEqual((output_dir / "checkpoint.json").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
