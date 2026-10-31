import json
import unittest
from unittest.mock import patch

import requests

from claudio.dataset_catalog.proteomecentral import fetch_dataset_pages


def make_response(page_number, records):
    payload = {
        "datasets": records,
        "facets": {},
        "query": {
            "search": "XL-MS",
            "resultType": "resultset",
            "pageNumber": page_number,
            "pageSize": 2,
        },
        "result_set": {
            "datasets_title_list": ["dataset identifier", "title"],
            "n_available_pages": 2,
            "n_available_rows": 3,
            "n_rows_returned": len(records),
            "page_number": page_number,
            "page_size": 2,
        },
        "status": {
            "status": "OK",
            "status_code": 200,
            "error_code": None,
            "description": f"Sent {len(records)} datasets",
        },
    }
    response = requests.Response()
    response.status_code = 200
    response.encoding = "utf-8"
    response.headers["Content-Type"] = "application/json"
    response.url = (
        "https://proteomecentral.proteomexchange.org"
        "/api/proxi/v0.1/datasets"
        f"?search=XL-MS&pageNumber={page_number}&pageSize=2"
    )
    response._content = json.dumps(payload).encode("utf-8")
    return response


class TestDatasetCatalogFetch(unittest.TestCase):
    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_fetches_all_pages_and_preserves_raw_responses(self, mock_get):
        responses = [
            make_response(
                1,
                [
                    ["PXD000001", "First study"],
                    ["PXD000002", "Second study"],
                ],
            ),
            make_response(2, [["PXD000003", "Third study"]]),
        ]
        mock_get.side_effect = responses

        pages = fetch_dataset_pages("XL-MS", page_size=2)

        accessions = [
            record[0]
            for page in pages
            for record in page.data["datasets"]
        ]
        self.assertEqual(
            accessions,
            ["PXD000001", "PXD000002", "PXD000003"],
        )
        self.assertEqual(
            [call.kwargs["params"]["pageNumber"]
             for call in mock_get.call_args_list],
            [1, 2],
        )
        self.assertEqual(pages[0].raw_body, responses[0].content)
        self.assertEqual(pages[1].raw_body, responses[1].content)
        self.assertEqual(pages[1].request_url, responses[1].url)
    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_preserves_empty_search_response(self, mock_get):
        response = make_response(1, [])
        payload = response.json()
        payload["result_set"]["n_available_rows"] = 0
        payload["result_set"]["n_available_pages"] = 0
        response._content = json.dumps(payload).encode("utf-8")
        mock_get.return_value = response

        pages = fetch_dataset_pages("XL-MS", page_size=2)

        self.assertEqual(len(pages), 1)
        self.assertEqual(pages[0].data["datasets"], [])
        self.assertEqual(pages[0].raw_body, response.content)

    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_propagates_http_error_from_later_page(self, mock_get):
        first = make_response(
            1,
            [
                ["PXD000001", "First study"],
                ["PXD000002", "Second study"],
            ],
        )
        second = make_response(2, [["PXD000003", "Third study"]])
        second.status_code = 503
        mock_get.side_effect = [first, second]

        with self.assertRaises(requests.exceptions.HTTPError):
            fetch_dataset_pages("XL-MS", page_size=2)

    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_rejects_incomplete_results(self, mock_get):
        mock_get.side_effect = [
            make_response(1, [["PXD000001", "First study"]]),
            make_response(2, [["PXD000003", "Third study"]]),
        ]

        with self.assertRaisesRegex(ValueError, "Incomplete download"):
            fetch_dataset_pages("XL-MS", page_size=2)

    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_rejects_api_error_in_http_200_response(self, mock_get):
        response = make_response(1, [])
        payload = response.json()
        payload["status"].update(
            status="ERROR",
            status_code=500,
            error_code="SEARCH_FAILED",
            description="Search failed",
        )
        response._content = json.dumps(payload).encode("utf-8")
        mock_get.return_value = response

        with self.assertRaisesRegex(ValueError, "search failed"):
            fetch_dataset_pages("XL-MS", page_size=2)

if __name__ == "__main__":
    unittest.main()