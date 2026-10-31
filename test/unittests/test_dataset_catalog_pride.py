import json
import unittest
from datetime import datetime
from unittest.mock import patch

import requests

from claudio.dataset_catalog.pride import (
    fetch_pride_file_pages,
    fetch_pride_project,
)


class TestDatasetCatalogPride(unittest.TestCase):
    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_fetches_project_metadata_and_preserves_raw_response(self, mock_get):
        payload = {
            "accession": "PXD035844",
            "title": "Test cross-linking dataset",
            "projectDescription": "Test project description",
            "sampleProcessingProtocol": "Test sample preparation",
            "dataProcessingProtocol": "Test identification workflow",
            "organisms": [
                {
                    "@type": "CvParam",
                    "cvLabel": "NEWT",
                    "accession": "NEWT:9606",
                    "name": "Homo sapiens (human)",
                }
            ],
            "references": [
                {"doi": "10.1073/pnas.2219418120", "pubmedID": 37071682}
            ],
            "submissionType": "PARTIAL",
        }
        raw_body = (
            json.dumps(payload, separators=(",", ":")) + "\n\n"
        ).encode("utf-8")
        response = requests.Response()
        response.status_code = 200
        response.url = (
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD035844"
        )
        response.headers["Content-Type"] = "application/json"
        response._content = raw_body
        mock_get.return_value = response

        project = fetch_pride_project("PXD035844")

        self.assertIsNotNone(project)
        self.assertEqual(project.data, payload)
        self.assertEqual(project.raw_body, raw_body)
        self.assertEqual(project.request_url, response.url)
        self.assertEqual(project.headers["Content-Type"], "application/json")
        retrieved_at = datetime.fromisoformat(project.retrieved_at_utc)
        self.assertIsNotNone(retrieved_at.utcoffset())
        self.assertEqual(retrieved_at.utcoffset().total_seconds(), 0)
        self.assertEqual(mock_get.call_args.args[0], response.url)
    def _response(self, body, status_code=200):
        response = requests.Response()
        response.status_code = status_code
        response.url = (
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD035844"
        )
        response.headers["Content-Type"] = "application/json"
        response._content = body
        return response

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_propagates_missing_project_error(self, mock_get):
        mock_get.return_value = self._response(
            b'{"status":404,"error":"Not Found"}', status_code=404
        )

        with self.assertRaises(requests.exceptions.HTTPError):
            fetch_pride_project("PXD035844")

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_rejects_wrong_or_missing_project_accession(self, mock_get):
        cases = [
            {"accession": "PXD013928", "title": "Another project"},
            {"title": "Response without an accession"},
        ]
        for payload in cases:
            with self.subTest(payload=payload):
                mock_get.return_value = self._response(
                    json.dumps(payload).encode("utf-8")
                )

                with self.assertRaisesRegex(ValueError, "accession mismatch"):
                    fetch_pride_project("PXD035844")

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_rejects_invalid_project_responses(self, mock_get):
        bodies = [
            b"",
            b"<html>Service unavailable</html>",
            b'{"accession":',
            b"[]",
            b"null",
        ]
        for body in bodies:
            with self.subTest(body=body):
                mock_get.return_value = self._response(body)

                with self.assertRaises(ValueError):
                    fetch_pride_project("PXD035844")
    def _file_response(self, page_number, records, total_records=3):
        response = requests.Response()
        response.status_code = 200
        response.headers["Content-Type"] = "application/json"
        response.headers["total_records"] = str(total_records)
        response.url = requests.Request(
            "GET",
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD035844/files",
            params={"pageSize": 2, "page": page_number},
        ).prepare().url
        response._content = (
            json.dumps(records, separators=(",", ":")) + "\n"
        ).encode("utf-8")
        return response

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_fetches_all_file_pages_without_losing_file_metadata(self, mock_get):
        first_records = [
            {
                "accession": "file-1",
                "projectAccessions": ["PXD035844"],
                "fileName": "sample.raw",
                "fileCategory": {"value": "RAW"},
                "fileSizeBytes": 100,
                "checksum": "a" * 40,
                "publicFileLocations": [
                    {"value": "ftp://example.org/sample.raw"}
                ],
            },
            {
                "accession": "file-2",
                "projectAccessions": ["PXD035844"],
                "fileName": "identifications.mzid",
                "fileCategory": {"value": "RESULT"},
                "fileSizeBytes": 200,
            },
        ]
        last_records = [
            {
                "accession": "file-3",
                "projectAccessions": ["PXD035844"],
                "fileName": "search.fasta",
                "fileCategory": {"value": "FASTA"},
                "fileSizeBytes": 300,
            },
        ]
        responses = [
            self._file_response(0, first_records),
            self._file_response(1, last_records),
        ]
        mock_get.side_effect = responses

        pages = fetch_pride_file_pages("PXD035844", page_size=2)

        self.assertEqual(
            [record["accession"] for page in pages for record in page.data],
            ["file-1", "file-2", "file-3"],
        )
        self.assertEqual(pages[0].data, first_records)
        self.assertEqual(pages[1].data, last_records)
        self.assertEqual(pages[0].raw_body, responses[0].content)
        self.assertEqual(pages[1].raw_body, responses[1].content)
        self.assertEqual(pages[1].request_url, responses[1].url)
        self.assertEqual(pages[1].headers["total_records"], "3")
        self.assertEqual(
            [call.kwargs["params"]["page"]
             for call in mock_get.call_args_list],
            [0, 1],
        )
        self.assertEqual(
            [call.kwargs["params"]["pageSize"]
             for call in mock_get.call_args_list],
            [2, 2],
        )
    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_preserves_empty_file_list_response(self, mock_get):
        response = self._file_response(0, [], total_records=0)
        mock_get.return_value = response

        pages = fetch_pride_file_pages("PXD035844", page_size=2)

        self.assertEqual(len(pages), 1)
        self.assertEqual(pages[0].data, [])
        self.assertEqual(pages[0].raw_body, response.content)
        self.assertEqual(pages[0].headers["total_records"], "0")
        self.assertEqual(mock_get.call_count, 1)

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_rejects_file_list_that_ends_before_reported_total(self, mock_get):
        record = {
            "accession": "file-1",
            "projectAccessions": ["PXD035844"],
        }
        mock_get.side_effect = [
            self._file_response(0, [record], total_records=2),
            self._file_response(1, [], total_records=2),
        ]

        with self.assertRaisesRegex(ValueError, "Incomplete PRIDE file list"):
            fetch_pride_file_pages("PXD035844", page_size=2)

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_rejects_repeated_files_within_or_across_pages(self, mock_get):
        record = {
            "accession": "file-1",
            "projectAccessions": ["PXD035844"],
        }
        for records_by_page in ([[record, record]], [[record], [record]]):
            with self.subTest(records_by_page=records_by_page):
                mock_get.side_effect = [
                    self._file_response(number, records, total_records=2)
                    for number, records in enumerate(records_by_page)
                ]

                with self.assertRaisesRegex(ValueError, "repeated file"):
                    fetch_pride_file_pages("PXD035844", page_size=2)

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_rejects_files_without_the_requested_project(self, mock_get):
        for projects in (["PXD013928"], [], None):
            with self.subTest(projects=projects):
                record = {
                    "accession": "file-1",
                    "projectAccessions": projects,
                }
                mock_get.return_value = self._file_response(
                    0, [record], total_records=1
                )

                with self.assertRaisesRegex(ValueError, "does not belong"):
                    fetch_pride_file_pages("PXD035844", page_size=2)

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_preserves_files_assigned_to_multiple_projects(self, mock_get):
        record = {
            "accession": "shared-file",
            "projectAccessions": ["PXD013928", "PXD035844"],
        }
        mock_get.return_value = self._file_response(
            0, [record], total_records=1
        )

        pages = fetch_pride_file_pages("PXD035844", page_size=2)

        self.assertEqual(len(pages), 1)
        self.assertEqual(pages[0].data, [record])

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_rejects_missing_or_invalid_file_counts(self, mock_get):
        for count in (None, "not-a-number", "-1"):
            with self.subTest(count=count):
                response = self._file_response(0, [], total_records=0)
                if count is None:
                    del response.headers["total_records"]
                else:
                    response.headers["total_records"] = count
                mock_get.return_value = response

                with self.assertRaises(ValueError):
                    fetch_pride_file_pages("PXD035844", page_size=2)

    @patch("claudio.dataset_catalog.pride.requests.get")
    def test_rejects_file_count_changes_between_pages(self, mock_get):
        mock_get.side_effect = [
            self._file_response(0, [{
                "accession": "file-1",
                "projectAccessions": ["PXD035844"],
            }], total_records=2),
            self._file_response(1, [{
                "accession": "file-2",
                "projectAccessions": ["PXD035844"],
            }], total_records=3),
        ]

        with self.assertRaisesRegex(ValueError, "file count changed"):
            fetch_pride_file_pages("PXD035844", page_size=2)
if __name__ == "__main__":
    unittest.main()