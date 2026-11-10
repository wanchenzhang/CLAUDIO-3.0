from _dataset_catalog_fixtures import install_pridepy_fixture

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
        "datasets": [list(row) + ['<a href="https://doi.org/10.1234/example">Article</a>'] for row in records],
        "facets": {},
        "query": dict(params),
        "result_set": {
            "datasets_title_list": [
                "dataset identifier", "title", "repository", "species", "publications"
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
    def setUp(self):
        install_pridepy_fixture(self)
        from claudio.dataset_catalog.discovery_history import initialize_discovery_history
        from claudio.dataset_catalog.publication_lookup import PublicationLookup
        from test_dataset_catalog_publication_lookup import Clock, Response, Session, crossref_body
        import os

        state = TemporaryDirectory()
        self.addCleanup(state.cleanup)
        initialize_discovery_history(Path(state.name) / "claudio/discovery")
        environment = patch.dict(os.environ, {"XDG_STATE_HOME": state.name})
        environment.start()
        self.addCleanup(environment.stop)

        def lookup(directory, **options):
            clock = Clock()
            return PublicationLookup(directory, session=Session(Response(crossref_body())),
                                     clock=clock, sleeper=clock.sleep, **options)

        publication_http = patch("claudio.dataset_catalog.published_discovery.PublicationLookup",
                                 side_effect=lookup)
        publication_http.start()
        self.addCleanup(publication_http.stop)

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


class TestDatasetCatalogInspectCLI(unittest.TestCase):
    def setUp(self):
        install_pridepy_fixture(self)
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.catalog = self.root / "candidates.json"
        self.output = self.root / "inventory"
        self.runner = CliRunner()

    def write_catalog(self, repositories):
        candidates = [
            {
                "accession": f"PXD{number:06d}",
                "repository": repository,
                "title": f"Study {number}",
                "search_terms": ["XL-MS"],
                "review_status": "pending",
            }
            for number, repository in enumerate(repositories, 1)
        ]
        self.catalog.write_text(json.dumps(candidates), encoding="utf-8")
        return candidates

    def invoke_inspect(self, *extra):
        return self.runner.invoke(main, [
            "inspect", "--catalog", str(self.catalog),
            "--output", str(self.output), "--delay", "0", *extra,
        ])

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_collects_an_inventory_and_reports_saved_progress(self, mock_get):
        candidates = self.write_catalog(["PRIDE"])
        url = "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD000001"
        project = requests.Response()
        project.status_code = 200
        project.url = url
        project._content = json.dumps({"accession": "PXD000001"}).encode()
        files = requests.Response()
        files.status_code = 200
        files.url = url + "/files"
        files.headers["total_records"] = "0"
        files._content = b"[]"
        mock_get.side_effect = [project, files]

        result = self.invoke_inspect()

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("Saved progress: 0/1", result.output)
        self.assertIn("[1/1] PXD000001 (PRIDE): collected; 0 file records", result.output)
        self.assertIn("Collected: 1", result.output)
        manifest = json.loads((self.output / "manifest.json").read_bytes())
        self.assertEqual(manifest["status_counts"]["collected"], 1)
        saved = json.loads((self.output / "projects.json").read_bytes())
        self.assertEqual(saved[0]["review_status"], "pending")
        self.assertEqual(saved[0]["search_terms"], candidates[0]["search_terms"])
        self.assertTrue((self.output / saved[0]["inventory_path"] / "files.csv").is_file())
        self.assertEqual(mock_get.call_count, 2)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_records_project_failures_and_continues(self, mock_get):
        self.write_catalog(["PRIDE", "UnknownRepository"])
        mock_get.side_effect = requests.exceptions.Timeout("Server timed out")

        result = self.invoke_inspect()

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("project_metadata: Server timed out", result.output)
        self.assertIn("Recorded outcomes: 2", result.output)
        self.assertIn("Collected: 0", result.output)
        self.assertIn("Failed: 1", result.output)
        self.assertIn("Adapter pending: 1", result.output)
        saved = json.loads((self.output / "projects.json").read_bytes())
        self.assertEqual([row["inspection_status"] for row in saved],
                         ["failed", "adapter_pending"])
        self.assertTrue(all(row["review_status"] == "pending" for row in saved))

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_existing_output_or_missing_resume_checkpoint(self, mock_get):
        self.write_catalog(["PRIDE"])
        self.output.mkdir()
        marker = self.output / "keep.txt"
        marker.write_text("keep", encoding="utf-8")

        for extra in ((), ("--resume",)):
            with self.subTest(extra=extra):
                result = self.invoke_inspect(*extra)
                self.assertEqual(result.exit_code, 2, result.output)
                self.assertEqual(list(self.output.iterdir()), [marker])
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_invalid_delays_before_creating_a_run(self, mock_get):
        self.write_catalog(["PRIDE"])
        for delay in ("-1", "nan", "inf"):
            with self.subTest(delay=delay):
                result = self.invoke_inspect("--delay", delay)
                self.assertEqual(result.exit_code, 2, result.output)
                self.assertFalse(self.output.exists())
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_resumes_after_interrupt_without_retrying_recorded_failures(self, mock_get):
        self.write_catalog(["PRIDE", "PRIDE"])
        mock_get.side_effect = [requests.exceptions.Timeout("Offline"), KeyboardInterrupt()]
        interrupted = self.invoke_inspect()
        self.assertNotEqual(interrupted.exit_code, 0)
        checkpoint = json.loads((self.output / "checkpoint.json").read_bytes())
        self.assertEqual(len(checkpoint["results"]), 1)

        mock_get.reset_mock()
        mock_get.side_effect = requests.exceptions.Timeout("Still offline")
        resumed = self.invoke_inspect("--resume")

        self.assertEqual(resumed.exit_code, 0, resumed.output)
        self.assertIn("Saved progress: 1/2", resumed.output)
        self.assertIn("Failed: 2", resumed.output)
        self.assertEqual(mock_get.call_count, 1)
        self.assertTrue(mock_get.call_args.args[0].endswith("PXD000002"))

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_rejects_changed_catalog_without_touching_saved_files(self, mock_get):
        candidates = self.write_catalog(["UnknownRepository"])
        first = self.invoke_inspect()
        self.assertEqual(first.exit_code, 0, first.output)
        before = {str(p.relative_to(self.output)): p.read_bytes()
                  for p in self.output.rglob("*") if p.is_file()}
        candidates[0]["title"] = "Changed study"
        self.catalog.write_text(json.dumps(candidates), encoding="utf-8")

        result = self.invoke_inspect("--resume")

        self.assertEqual(result.exit_code, 1, result.output)
        self.assertIn("catalog", result.output.lower())
        after = {str(p.relative_to(self.output)): p.read_bytes()
                 for p in self.output.rglob("*") if p.is_file()}
        self.assertEqual(after, before)
        mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_completed_resume_leaves_reports_unchanged(self, mock_get):
        self.write_catalog(["UnknownRepository"])
        first = self.invoke_inspect()
        self.assertEqual(first.exit_code, 0, first.output)
        before = {str(p.relative_to(self.output)): p.read_bytes()
                  for p in self.output.rglob("*") if p.is_file()}

        result = self.invoke_inspect("--resume")

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertIn("Saved progress: 1/1", result.output)
        after = {str(p.relative_to(self.output)): p.read_bytes()
                 for p in self.output.rglob("*") if p.is_file()}
        self.assertEqual(after, before)
        mock_get.assert_not_called()


class TestDatasetCatalogSummaryCLI(unittest.TestCase):
    def setUp(self):
        install_pridepy_fixture(self)

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_summarizes_saved_outcomes_without_network_requests(self, mock_get):
        from claudio.dataset_catalog.inventory import inspect_catalog

        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog = root / "candidates.json"
            catalog.write_text(json.dumps([{
                "accession": "PXD000001", "repository": "UnknownRepository",
                "review_status": "pending",
            }]), encoding="utf-8")
            run = root / "run"
            inspect_catalog(catalog, run)
            result = CliRunner().invoke(main, [
                "summarize", "--catalog", str(catalog),
                "--run", str(run), "--output", str(root / "summary"),
            ])
            self.assertEqual(result.exit_code, 0, result.output)
            self.assertIn("Candidate datasets: 1", result.output)
            self.assertIn("Recorded attempts: 1", result.output)
            self.assertIn("Adapter pending: 1", result.output)
            self.assertTrue((root / "summary/attempts.json").is_file())
            mock_get.assert_not_called()

    def test_requires_source_runs_and_reports_incomplete_coverage(self):
        from claudio.dataset_catalog.inventory import inspect_catalog

        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate = {"accession": "PXD000001", "repository": "UnknownRepository"}
            catalog = root / "candidates.json"
            catalog.write_text(json.dumps([candidate]), encoding="utf-8")
            run = root / "run"
            inspect_catalog(catalog, run)
            catalog.write_text(json.dumps([candidate, {
                "accession": "PXD000002", "repository": "UnknownRepository",
            }]), encoding="utf-8")
            args = ["summarize", "--catalog", str(catalog), "--output", str(root / "summary")]
            result = CliRunner().invoke(main, args)
            self.assertEqual(result.exit_code, 2, result.output)
            result = CliRunner().invoke(main, args + ["--run", str(run)])
            self.assertEqual(result.exit_code, 1, result.output)
            self.assertFalse((root / "summary").exists())

if __name__ == "__main__":
    unittest.main()