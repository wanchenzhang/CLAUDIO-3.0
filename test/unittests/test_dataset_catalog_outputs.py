import csv
import hashlib
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from claudio.dataset_catalog.proteomecentral import write_discovery_outputs
from claudio.dataset_catalog.proteomecentral import CatalogPage


class TestDatasetCatalogOutputs(unittest.TestCase):
    def test_writes_candidates_and_auditable_raw_response(self):
        payload = {
            "datasets": [
                ["PXD000001", "Cross-linking study", "PRIDE", "Homo sapiens"]
            ],
            "facets": {},
            "query": {
                "search": "XL-MS",
                "resultType": "resultset",
                "pageNumber": 1,
                "pageSize": 100,
            },
            "result_set": {
                "datasets_title_list": [
                    "dataset identifier", "title", "repository", "species"
                ],
                "n_available_pages": 1,
                "n_available_rows": 1,
                "n_rows_returned": 1,
                "page_number": 1,
                "page_size": 100,
            },
            "status": {
                "status": "OK",
                "status_code": 200,
                "error_code": None,
                "description": "Sent 1 datasets",
            },
        }
        raw_body = (
            json.dumps(payload, separators=(",", ":")) + "\n\n"
        ).encode("utf-8")
        page = CatalogPage(
            request_url=(
                "https://proteomecentral.proteomexchange.org"
                "/api/proxi/v0.1/datasets?search=XL-MS"
            ),
            retrieved_at_utc="2026-10-31T13:00:00+00:00",
            raw_body=raw_body,
            data=payload,
        )

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            write_discovery_outputs({"XL-MS": [page]}, output_dir)

            self.assertTrue((output_dir / "manifest.json").is_file())
            self.assertTrue((output_dir / "candidates.json").is_file())
            self.assertTrue((output_dir / "candidates.csv").is_file())

            candidates = json.loads(
                (output_dir / "candidates.json").read_text(encoding="utf-8")
            )
            self.assertEqual(len(candidates), 1)
            self.assertEqual(candidates[0]["accession"], "PXD000001")
            self.assertEqual(
                candidates[0]["metadata_variants"][0]["metadata"]["title"],
                "Cross-linking study",
            )

            with (output_dir / "candidates.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["accession"], "PXD000001")
            self.assertEqual(json.loads(rows[0]["search_terms"]), ["XL-MS"])
            self.assertEqual(rows[0]["review_status"], "pending")

            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["candidate_count"], 1)
            search = manifest["searches"][0]
            self.assertEqual(search["search_term"], "XL-MS")
            self.assertEqual(search["reported_total"], 1)
            self.assertEqual(search["downloaded_records"], 1)

            saved_page = search["pages"][0]
            raw_path = Path(saved_page["path"])
            self.assertFalse(raw_path.is_absolute())
            self.assertEqual((output_dir / raw_path).read_bytes(), raw_body)
            self.assertEqual(
                saved_page["sha256"], hashlib.sha256(raw_body).hexdigest()
            )
            self.assertEqual(saved_page["request_url"], page.request_url)
            self.assertEqual(
                saved_page["retrieved_at_utc"], page.retrieved_at_utc
            )
    def _empty_page(self):
        payload = {
            "datasets": [],
            "facets": {},
            "query": {
                "search": "XL-MS",
                "resultType": "resultset",
                "pageNumber": 1,
                "pageSize": 100,
            },
            "result_set": {
                "datasets_title_list": ["dataset identifier", "title"],
                "n_available_pages": 0,
                "n_available_rows": 0,
                "n_rows_returned": 0,
                "page_number": 1,
                "page_size": 100,
            },
            "status": {
                "status": "OK",
                "status_code": 200,
                "error_code": None,
                "description": "Sent 0 datasets",
            },
        }
        return CatalogPage(
            request_url=(
                "https://proteomecentral.proteomexchange.org"
                "/api/proxi/v0.1/datasets?search=XL-MS"
            ),
            retrieved_at_utc="2026-10-31T13:00:00+00:00",
            raw_body=json.dumps(payload).encode("utf-8"),
            data=payload,
        )

    def test_writes_header_and_search_record_when_no_candidates_match(self):
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            write_discovery_outputs(
                {"XL-MS": [self._empty_page()]}, output_dir
            )

            with (output_dir / "candidates.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                reader = csv.DictReader(handle)
                self.assertIn("accession", reader.fieldnames)
                self.assertEqual(list(reader), [])
            candidates = json.loads(
                (output_dir / "candidates.json").read_text(encoding="utf-8")
            )
            self.assertEqual(candidates, [])

            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["candidate_count"], 0)
            search = manifest["searches"][0]
            self.assertEqual(search["reported_total"], 0)
            self.assertEqual(search["downloaded_records"], 0)
            self.assertEqual(len(search["pages"]), 1)

    def test_refuses_to_overwrite_an_existing_snapshot(self):
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            old_csv = output_dir / "candidates.csv"
            old_csv.write_text("previous results\n", encoding="utf-8")

            with self.assertRaises(FileExistsError):
                write_discovery_outputs(
                    {"XL-MS": [self._empty_page()]}, output_dir
                )

            self.assertEqual(
                old_csv.read_text(encoding="utf-8"), "previous results\n"
            )
            self.assertFalse((output_dir / "manifest.json").exists())

    def test_write_failure_does_not_publish_a_complete_manifest(self):
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            with patch(
                "claudio.dataset_catalog.proteomecentral.csv.DictWriter.writeheader",
                side_effect=OSError("No space left on device"),
            ):
                with self.assertRaises(OSError):
                    write_discovery_outputs(
                        {"XL-MS": [self._empty_page()]}, output_dir
                    )

            self.assertTrue((output_dir / "candidates.json").is_file())
            self.assertFalse((output_dir / "manifest.json").exists())

    def test_missing_responses_are_not_treated_as_zero_hits(self):
        for searches in ({}, {"XL-MS": []}):
            with self.subTest(searches=searches):
                with TemporaryDirectory() as temp_dir:
                    output_dir = Path(temp_dir) / "run"

                    with self.assertRaises(ValueError):
                        write_discovery_outputs(searches, output_dir)

                    self.assertFalse(output_dir.exists())

if __name__ == "__main__":
    unittest.main()