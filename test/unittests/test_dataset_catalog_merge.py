import json
import unittest

from claudio.dataset_catalog.discovery import merge_search_results
from claudio.dataset_catalog.proteomecentral import CatalogPage


def make_page(search_term, records):
    payload = {
        "datasets": records,
        "facets": {},
        "query": {
            "search": search_term,
            "resultType": "resultset",
            "pageNumber": 1,
            "pageSize": 100,
        },
        "result_set": {
            "datasets_title_list": [
                "dataset identifier", "title", "repository", "species"
            ],
            "n_available_pages": 1,
            "n_available_rows": len(records),
            "n_rows_returned": len(records),
            "page_number": 1,
            "page_size": 100,
        },
        "status": {
            "status": "OK",
            "status_code": 200,
            "error_code": None,
            "description": f"Sent {len(records)} datasets",
        },
    }
    return CatalogPage(
        request_url=(
            "https://proteomecentral.proteomexchange.org"
            f"/api/proxi/v0.1/datasets?search={search_term}"
        ),
        retrieved_at_utc="2026-10-31T00:00:00+00:00",
        raw_body=json.dumps(payload).encode("utf-8"),
        data=payload,
    )


class TestDatasetCatalogMerge(unittest.TestCase):
    def test_merges_by_accession_and_keeps_each_search_term(self):
        searches = {
            "XL-MS": [
                make_page("XL-MS", [
                    ["PXD000001", "First study", "PRIDE", "Homo sapiens"],
                    ["PXD000002", "Second study", "MassIVE", "Escherichia coli"],
                ])
            ],
            "crosslinking": [
                make_page("crosslinking", [
                    ["PXD000002", "Second study", "MassIVE", "Escherichia coli"],
                    ["PXD000003", "Third study", "PRIDE", "Homo sapiens"],
                ])
            ],
        }

        records = merge_search_results(searches)

        self.assertEqual(
            [record["accession"] for record in records],
            ["PXD000001", "PXD000002", "PXD000003"],
        )
        self.assertEqual(records[0]["search_terms"], ["XL-MS"])
        self.assertEqual(
            records[1]["search_terms"], ["XL-MS", "crosslinking"]
        )
        self.assertEqual(records[2]["search_terms"], ["crosslinking"])
        self.assertEqual(records[1]["title"], "Second study")
        self.assertEqual(records[1]["repository"], "MassIVE")
        self.assertEqual(records[1]["species"], "Escherichia coli")
        self.assertTrue(
            all(record["review_status"] == "pending" for record in records)
        )

    def test_preserves_conflicting_metadata_with_its_search_sources(self):
        searches = {
            "XL-MS": [
                make_page("XL-MS", [
                    ["PXD000001", "Original title", "PRIDE", "Homo sapiens"],
                ])
            ],
            "crosslinking": [
                make_page("crosslinking", [
                    ["PXD000001", "Updated title", "PRIDE", "Homo sapiens"],
                ])
            ],
        }

        records = merge_search_results(searches)

        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertIn("metadata_conflict", record)
        self.assertTrue(record["metadata_conflict"])
        self.assertEqual(record["review_status"], "pending")
        self.assertEqual(
            record["search_terms"], ["XL-MS", "crosslinking"]
        )
        self.assertEqual(
            [variant["metadata"]["title"]
             for variant in record["metadata_variants"]],
            ["Original title", "Updated title"],
        )
        self.assertEqual(
            [variant["search_terms"]
             for variant in record["metadata_variants"]],
            [["XL-MS"], ["crosslinking"]],
        )

    def test_combines_sources_for_identical_metadata(self):
        records = [["PXD000001", "First study", "PRIDE", "Homo sapiens"]]
        searches = {
            "XL-MS": [make_page("XL-MS", records)],
            "crosslinking": [make_page("crosslinking", records)],
        }

        merged = merge_search_results(searches)

        self.assertEqual(len(merged), 1)
        record = merged[0]
        self.assertIn("metadata_conflict", record)
        self.assertFalse(record["metadata_conflict"])
        self.assertEqual(len(record["metadata_variants"]), 1)
        self.assertEqual(
            record["metadata_variants"][0]["search_terms"],
            ["XL-MS", "crosslinking"],
        )
if __name__ == "__main__":
    unittest.main()