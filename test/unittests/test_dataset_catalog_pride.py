import json
import unittest
from datetime import datetime
from unittest.mock import patch

import requests

from claudio.dataset_catalog.pride import fetch_pride_project


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

if __name__ == "__main__":
    unittest.main()