import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import requests
from click.testing import CliRunner

from claudio.dataset_catalog.cli import main


def make_response(url, params, records):
    payload = {
        "datasets": records,
        "facets": {},
        "query": dict(params),
        "result_set": {
            "datasets_title_list": [
                "dataset identifier", "title", "repository", "species"
            ],
            "n_available_pages": 1 if records else 0,
            "n_available_rows": len(records),
            "n_rows_returned": len(records),
            "page_number": 1,
            "page_size": params["pageSize"],
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
    response.url = requests.Request("GET", url, params=params).prepare().url
    response._content = json.dumps(payload).encode("utf-8")
    return response


class TestDatasetCatalogCLI(unittest.TestCase):
    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_discovers_multiple_terms_and_writes_merged_candidates(self, mock_get):
        records_by_term = {
            "XL-MS": [
                ["PXD000001", "First study", "PRIDE", "Homo sapiens"],
            ],
            "crosslinking": [
                ["PXD000001", "First study", "PRIDE", "Homo sapiens"],
                ["PXD000002", "Second study", "MassIVE", "Escherichia coli"],
            ],
        }

        def respond(url, *, params, timeout):
            return make_response(url, params, records_by_term[params["search"]])

        mock_get.side_effect = respond

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            result = CliRunner().invoke(main, [
                "discover",
                "--search-term", "XL-MS",
                "--search-term", "crosslinking",
                "--output", str(output_dir),
            ])

            self.assertEqual(
                result.exit_code, 0, f"{result.output}\n{result.exception}"
            )
            candidates = json.loads(
                (output_dir / "candidates.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                [record["accession"] for record in candidates],
                ["PXD000001", "PXD000002"],
            )
            self.assertEqual(
                candidates[0]["search_terms"], ["XL-MS", "crosslinking"]
            )
            self.assertEqual(
                candidates[1]["search_terms"], ["crosslinking"]
            )

            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["candidate_count"], 2)
            self.assertEqual(
                [search["search_term"] for search in manifest["searches"]],
                ["XL-MS", "crosslinking"],
            )
            for search in manifest["searches"]:
                for page in search["pages"]:
                    self.assertTrue((output_dir / page["path"]).is_file())

        self.assertEqual(
            [call.kwargs["params"]["search"]
             for call in mock_get.call_args_list],
            ["XL-MS", "crosslinking"],
        )
    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_queries_trimmed_duplicate_terms_only_once(self, mock_get):
        def respond(url, *, params, timeout):
            return make_response(url, params, [
                ["PXD000001", "First study", "PRIDE", "Homo sapiens"],
            ])

        mock_get.side_effect = respond

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            result = CliRunner().invoke(main, [
                "discover",
                "--search-term", " XL-MS ",
                "--search-term", "XL-MS",
                "--output", str(output_dir),
            ])

            self.assertEqual(result.exit_code, 0, result.output)
            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                [search["search_term"] for search in manifest["searches"]],
                ["XL-MS"],
            )
            self.assertEqual(
                [call.kwargs["params"]["search"]
                 for call in mock_get.call_args_list],
                ["XL-MS"],
            )

    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_rejects_invalid_arguments_before_querying(self, mock_get):
        cases = [
            (["--search-term", "   "], "--search-term"),
            (["--search-term", "XL-MS", "--page-size", "0"], "--page-size"),
        ]
        for arguments, option in cases:
            with self.subTest(arguments=arguments):
                with TemporaryDirectory() as temp_dir:
                    output_dir = Path(temp_dir) / "run"
                    result = CliRunner().invoke(main, [
                        "discover", "--output", str(output_dir),
                        *arguments,
                    ])

                    self.assertEqual(result.exit_code, 2, result.output)
                    self.assertIn(option, result.output)
                    self.assertFalse(output_dir.exists())

        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_existing_output_is_rejected_before_querying(self, mock_get):
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            old_csv = output_dir / "candidates.csv"
            old_csv.write_text("previous results\n", encoding="utf-8")

            result = CliRunner().invoke(main, [
                "discover",
                "--search-term", "XL-MS",
                "--output", str(output_dir),
            ])

            self.assertEqual(result.exit_code, 2, result.output)
            self.assertIn("--output", result.output)
            self.assertEqual(
                old_csv.read_text(encoding="utf-8"), "previous results\n"
            )
            mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_later_query_failure_does_not_publish_a_complete_snapshot(self, mock_get):
        def respond(url, *, params, timeout):
            if params["search"] == "crosslinking":
                raise requests.exceptions.Timeout("read timed out")
            return make_response(url, params, [
                ["PXD000001", "First study", "PRIDE", "Homo sapiens"],
            ])

        mock_get.side_effect = respond

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            result = CliRunner().invoke(main, [
                "discover",
                "--search-term", "XL-MS",
                "--search-term", "crosslinking",
                "--output", str(output_dir),
            ])

            self.assertEqual(result.exit_code, 1, result.output)
            self.assertIn("read timed out", result.output)
            self.assertFalse((output_dir / "manifest.json").exists())

        self.assertEqual(
            [call.kwargs["params"]["search"]
             for call in mock_get.call_args_list],
            ["XL-MS", "crosslinking"],
        )

    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_zero_hits_are_saved_as_a_successful_empty_search(self, mock_get):
        def respond(url, *, params, timeout):
            return make_response(url, params, [])

        mock_get.side_effect = respond

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            result = CliRunner().invoke(main, [
                "discover",
                "--search-term", "XL-MS",
                "--output", str(output_dir),
            ])

            self.assertEqual(result.exit_code, 0, result.output)
            candidates = json.loads(
                (output_dir / "candidates.json").read_text(encoding="utf-8")
            )
            self.assertEqual(candidates, [])
            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["candidate_count"], 0)
            self.assertEqual(len(manifest["searches"][0]["pages"]), 1)
    @patch("claudio.dataset_catalog.proteomecentral.requests.get")
    def test_uses_builtin_terms_when_none_are_supplied(self, mock_get):
        def respond(url, *, params, timeout):
            return make_response(url, params, [
                ["PXD000001", "First study", "PRIDE", "Homo sapiens"],
            ])

        mock_get.side_effect = respond
        expected_terms = [
            "XL-MS",
            "crosslinking",
            "cross-linking",
            "cross linking",
            "crosslinked",
        ]

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            result = CliRunner().invoke(main, [
                "discover", "--output", str(output_dir),
            ])

            self.assertEqual(result.exit_code, 0, result.output)
            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                [search["search_term"] for search in manifest["searches"]],
                expected_terms,
            )
            candidates = json.loads(
                (output_dir / "candidates.json").read_text(encoding="utf-8")
            )
            self.assertEqual(len(candidates), 1)
            self.assertEqual(candidates[0]["search_terms"], expected_terms)

        self.assertEqual(
            [call.kwargs["params"]["search"]
             for call in mock_get.call_args_list],
            expected_terms,
        )
if __name__ == "__main__":
    unittest.main()