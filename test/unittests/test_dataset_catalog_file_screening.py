from _dataset_catalog_fixtures import install_pridepy_fixture

import csv
import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import requests
from click.testing import CliRunner

from claudio.dataset_catalog.cli import main
from claudio.dataset_catalog.inventory import inspect_catalog, summarize_inventories
from claudio.dataset_catalog.repositories import classify_inventory_file


def file_record(name, category="", **extra):
    return {
        "accession": "file-" + name,
        "fileName": name,
        "fileCategory": {"value": category},
        "projectAccessions": ["PXD000001"],
        "fileSizeBytes": 100,
        "publicFileLocations": [{"value": "https://example.org/data/" + name}],
        **extra,
    }


class TestFileRoleClassification(unittest.TestCase):
    def setUp(self):
        install_pridepy_fixture(self)

    def test_distinguishes_results_databases_spectra_and_supporting_files(self):
        cases = [
            ("pLink2_filtered_crosslinks.csv", "", "identification_result_candidate", "first"),
            ("experiment_xiFDR_Links.csv", "", "identification_result_candidate", "first"),
            ("identifications.mzid", "", "identification_result_candidate", "review"),
            ("export.csv", "RESULT", "identification_result_candidate", "review"),
            ("export.dat", "Search engine output file URI", "identification_result_candidate", "review"),
            ("sample.dat", "Associated raw file URI", "mass_spectra", "defer"),
            ("search.FASTA.gz", "", "search_fasta_candidate", "first"),
            ("sample.raw", "", "mass_spectra", "defer"),
            ("sample.mgf", "", "mass_spectra", "defer"),
            ("sample.mzML.gz", "", "mass_spectra", "defer"),
            ("supplementary_table_1.xlsx", "", "supplementary_candidate", "review"),
            ("methods.pdf", "", "supplementary_candidate", "review"),
            ("table.csv", "OTHER", "uncertain", "review"),
            ("task.plink", "", "uncertain", "review"),
            ("crosslinks.csv.zip", "RESULT", "archive", "review"),
            ("search.fasta.tar.gz", "FASTA", "archive", "review"),
        ]
        for name, category, role, priority in cases:
            with self.subTest(name=name, category=category):
                record = file_record(name, category)
                original = deepcopy(record)
                result = classify_inventory_file(record)
                self.assertEqual(result["file_role"], role)
                self.assertEqual(result["download_priority"], priority)
                self.assertTrue(result["evidence"])
                self.assertEqual(record, original)

    def test_initial_downloads_require_xl_hinted_text_tables_or_fasta(self):
        for name in ("crosslinks.csv.gz", "experiment_xiFDR_links.tsv", "pLink2_results.txt"):
            with self.subTest(name=name):
                result = classify_inventory_file(file_record(name, "RESULT"))
                self.assertEqual(result["download_priority"], "first")
        for name in ("crosslinks.msf", "crosslinks.pdresult", "crosslinks.zhrm",
                     "crosslinks.xml", "crosslinks.xlsx", "crosslinks.html", "crosslinks.svg"):
            with self.subTest(name=name):
                record = file_record(name, "RESULT")
                original = deepcopy(record)
                result = classify_inventory_file(record)
                self.assertEqual(result["file_role"], "identification_result_candidate")
                self.assertEqual(result["download_priority"], "review")
                self.assertTrue(result["xl_specific_hint"])
                self.assertEqual(record, original)

    def test_generic_results_are_not_confirmed_as_crosslink_identifications(self):
        result = classify_inventory_file(file_record("export.csv", "RESULT"))
        self.assertFalse(result["xl_specific_hint"])
        self.assertEqual(result["download_priority"], "review")
        self.assertEqual(result["classification_status"], "candidate")
        self.assertNotIn("fdr", result)
        self.assertNotIn("eligible", result)

    def test_a_crosslink_directory_does_not_make_every_table_a_result(self):
        result = classify_inventory_file(file_record(
            "experimental_design.csv", "OTHER",
            relativePath="crosslink_results/experimental_design.csv",
        ))
        self.assertEqual(result["file_role"], "uncertain")
        self.assertFalse(result["xl_specific_hint"])

    def test_conflicting_file_type_and_repository_category_stay_uncertain(self):
        result = classify_inventory_file(file_record("sample.raw", "RESULT"))
        self.assertEqual(result["file_role"], "uncertain")
        self.assertEqual(result["classification_status"], "conflicting_metadata")
        self.assertEqual(set(result["possible_roles"]),
                         {"mass_spectra", "identification_result_candidate"})
        self.assertEqual(result["download_priority"], "review")

    def test_raw_data_and_peak_lists_remain_distinguishable(self):
        raw = classify_inventory_file(file_record("sample.raw"))
        peak = classify_inventory_file(file_record("sample.mgf"))
        self.assertEqual(raw["spectra_kind"], "raw")
        self.assertEqual(peak["spectra_kind"], "peak_list")

    def test_missing_links_unresolved_paths_and_empty_files_block_download_candidates(self):
        cases = [
            ({"publicFileLocations": []}, "no_public_download_url"),
            ({"pathDecodingStatus": "unresolved"}, "unresolved_path_encoding"),
            ({"fileSizeBytes": 0}, "zero_byte_file"),
        ]
        for fields, reason in cases:
            with self.subTest(reason=reason):
                result = classify_inventory_file(file_record("crosslinks.csv", **fields))
                self.assertEqual(result["file_role"], "identification_result_candidate")
                self.assertEqual(result["download_priority"], "blocked")
                self.assertIn(reason, result["download_blockers"])

    def test_an_unresolved_reported_link_is_not_repaired(self):
        uri = "https://download.iprox.org/IPX0015485000/IPX0015485001/results_1%FDR.plink"
        result = classify_inventory_file(file_record(
            "results_1%FDR.plink", pathDecodingStatus="unresolved",
            publicFileLocations=[{"value": uri}],
        ))
        self.assertEqual(result["download_urls"], [uri])
        self.assertEqual(result["download_priority"], "blocked")


class TestSummaryFileScreening(unittest.TestCase):
    def setUp(self):
        install_pridepy_fixture(self)

    def make_run(self, root):
        candidates = [
            {"accession": "PXD000001", "repository": "PRIDE", "review_status": "pending"},
            {"accession": "PXD000002", "repository": "PRIDE", "review_status": "pending"},
        ]
        catalog = root / "candidates.json"
        catalog.write_text(json.dumps(candidates), encoding="utf-8")
        files = [
            file_record("crosslinks.csv"), file_record("search.fasta.gz"),
            file_record("sample.raw"), file_record("supplementary_table_1.xlsx"),
            file_record("unknown_table.csv"), file_record("results.zip", "RESULT"),
            file_record("blocked_crosslinks.csv", pathDecodingStatus="unresolved"),
        ]

        def respond(url, **kwargs):
            if "PXD000002" in url:
                raise requests.exceptions.HTTPError("404 Not Found")
            response = requests.Response()
            response.status_code = 200
            response.url = requests.Request("GET", url, params=kwargs.get("params")).prepare().url
            response.headers["Content-Type"] = "application/json"
            if url.endswith("/files"):
                response.headers["total_records"] = str(len(files))
                response._content = json.dumps(files).encode("utf-8")
            else:
                response._content = b'{"accession":"PXD000001"}'
            return response

        run = root / "run"
        with patch("claudio.dataset_catalog.repositories.requests.get", side_effect=respond):
            inspect_catalog(catalog, run)
        return catalog, run, files

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_summary_adds_file_reports_without_dropping_uncertain_files_or_failed_projects(self, mock_get):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog, run, files = self.make_run(root)
            before = {str(p.relative_to(run)): p.read_bytes() for p in run.rglob("*") if p.is_file()}
            output = summarize_inventories(catalog, [run], root / "summary")
            with (output / "file_screening.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), len(files))
            self.assertEqual([row["source_file_record_index"] for row in rows],
                             [str(number) for number in range(len(files))])
            with (output / "download_candidates.csv").open(encoding="utf-8", newline="") as handle:
                downloads = list(csv.DictReader(handle))
            self.assertEqual([row["file_name"] for row in downloads],
                             ["crosslinks.csv", "search.fasta.gz"])
            self.assertTrue(all(row["download_priority"] == "first" for row in downloads))
            self.assertEqual(rows[-1]["download_priority"], "blocked")
            self.assertEqual(json.loads(rows[0]["download_urls"]),
                             [files[0]["publicFileLocations"][0]["value"]])
            source_files = Path(rows[0]["source_inventory_path"]) / "files.json"
            self.assertEqual(rows[0]["source_files_sha256"], hashlib.sha256(source_files.read_bytes()).hexdigest())
            projects = json.loads((output / "project_file_screening.json").read_bytes())
            self.assertEqual(len(projects), 2)
            self.assertEqual(projects[0]["screening_status"], "metadata_screened")
            self.assertEqual(projects[1]["screening_status"], "inventory_unavailable")
            self.assertIsNone(projects[1]["listed_file_count"])
            self.assertEqual(projects[1]["review_status"], "pending")
            manifest = json.loads((output / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_screening"]["file_count"], 7)
            self.assertEqual(manifest["file_screening"]["download_candidate_count"], 2)
            self.assertEqual(manifest["file_screening"]["basis"], "metadata_only")
            self.assertEqual(manifest["scientific_review"], "unchanged")
            for name in ("file_screening.csv", "download_candidates.csv", "project_file_screening.json"):
                self.assertEqual(manifest["output_sha256"][name],
                                 hashlib.sha256((output / name).read_bytes()).hexdigest())
            self.assertEqual(before, {str(p.relative_to(run)): p.read_bytes()
                                      for p in run.rglob("*") if p.is_file()})
            mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories.requests.get")
    def test_existing_summary_command_reports_the_screening_counts(self, mock_get):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalog, run, _ = self.make_run(root)
            result = CliRunner().invoke(main, [
                "summarize", "--catalog", str(catalog), "--run", str(run),
                "--output", str(root / "summary"),
            ])
            self.assertEqual(result.exit_code, 0, result.output)
            self.assertIn("File records screened: 7", result.output)
            self.assertIn("Priority download candidates: 2", result.output)
            mock_get.assert_not_called()


if __name__ == "__main__":
    unittest.main()
