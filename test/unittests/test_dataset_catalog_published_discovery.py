
"""Offline end-to-end contracts for published incremental discovery."""

import copy
import csv
import hashlib
import importlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from urllib.parse import parse_qs, unquote, urlsplit

from click.testing import CliRunner
import requests

from claudio.dataset_catalog.discovery_history import DiscoveryHistory, initialize_discovery_history
from claudio.dataset_catalog.proteomecentral import CatalogPage, DATASETS_URL
from claudio.dataset_catalog.publication_lookup import PublicationLookup as RealLookup, replay_publication_lookup
from test_dataset_catalog_publication_lookup import Clock, Response, crossref_body, pubmed_body


DOI = "10.1234/example"


def record(accession="PXD000001", publications=None):
    return {"accession": accession, "title": "An XL-MS dataset", "repository": "PRIDE",
            "species": "Homo sapiens", "publications": publications if publications is not None else
            '<a href="https://doi.org/10.1234/example">Article</a>'}


def page(term, records, page_size=100):
    columns = ["dataset identifier", "title", "repository", "species", "publications", "doi"]
    params = {"search": term, "resultType": "resultset", "pageNumber": 1, "pageSize": page_size}
    data = {
        "datasets": [[row["accession"], row["title"], row["repository"], row["species"],
                      row["publications"], row.get("doi", "")] for row in records],
        "facets": {}, "query": params,
        "result_set": {"datasets_title_list": columns, "n_available_pages": 1 if records else 0,
                       "n_available_rows": len(records), "n_rows_returned": len(records),
                       "page_number": 1, "page_size": page_size},
        "status": {"status": "OK", "status_code": 200, "error_code": None, "description": "Offline fixture"},
    }
    return CatalogPage(requests.Request("GET", DATASETS_URL, params=params).prepare().url,
                       "2026-11-09T12:00:00+00:00", json.dumps(data).encode(), data)


class ProviderSession:
    def __init__(self, calls, responses=None):
        self.calls = calls
        self.responses = responses or {}

    def get(self, url, *, timeout, allow_redirects, stream, headers):
        if timeout != (10, 30) or allow_redirects is not False or stream is not True:
            raise AssertionError("Unexpected publication request policy")
        self.calls.append(url)
        parsed = urlsplit(url)
        if parsed.hostname == "api.crossref.org":
            identifier = unquote(parsed.path.removeprefix("/works/"))
            payload = json.loads(crossref_body())
            payload["message"]["DOI"] = identifier
            default = json.dumps(payload).encode()
        elif parsed.hostname == "eutils.ncbi.nlm.nih.gov":
            identifier = parse_qs(parsed.query)["id"][0]
            default = pubmed_body()
        else:
            raise AssertionError("Publication request left the official endpoints")
        response = self.responses.get(identifier, default)
        if isinstance(response, Exception):
            raise response
        if isinstance(response, dict):
            response = json.dumps(response).encode()
        return Response(response, url=url)

    def close(self):
        pass


class TestPublishedDiscoveryIntegration(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.history = self.root / "history"
        initialize_discovery_history(self.history)
        self.http_calls = []
        self.search_calls = []

    def api(self):
        try:
            return importlib.import_module("claudio.dataset_catalog.published_discovery")
        except ModuleNotFoundError as error:
            if error.name != "claudio.dataset_catalog.published_discovery":
                raise
            self.fail("Published discovery integration is not implemented")

    def factory(self, responses=None):
        def create(directory, **options):
            clock = Clock()
            return RealLookup(directory, session=ProviderSession(self.http_calls, responses),
                              clock=clock, sleeper=clock.sleep, **options)
        return create

    def fetcher(self, records):
        def fetch(term, *, page_size):
            self.search_calls.append(term)
            return [page(term, records, page_size)]
        return fetch

    def run_discovery(self, records, name="run", responses=None, **options):
        api = self.api()
        with patch.object(api, "fetch_dataset_pages", side_effect=self.fetcher(records)), \
             patch.object(api, "PublicationLookup", side_effect=self.factory(responses)):
            return api.run_published_discovery(
                ["XL-MS"], self.root / name, history_dir=self.history, as_of="2026-11-09", **options)

    def read(self, directory, name):
        return json.loads((Path(directory) / name).read_bytes())

    def seen(self):
        with DiscoveryHistory(self.history) as history:
            return history.seen_accessions()

    def test_pipeline_excludes_seen_published_candidates_on_second_run(self):
        first = self.run_discovery([record()], "first")
        second = self.run_discovery([record()], "second")
        self.assertEqual([row["accession"] for row in self.read(first, "candidates.json")], ["PXD000001"])
        self.assertEqual(self.read(second, "candidates.json"), [])
        self.assertEqual(self.read(second, "candidate_audit.json")[0]["discovery_decision"], "excluded_seen")
        self.assertEqual(len(self.http_calls), 1)

    def test_later_publication_can_enter_after_an_earlier_pending_run(self):
        first = self.run_discovery([record(publications="Pending")], "first")
        self.assertEqual(self.read(first, "candidates.json"), [])
        self.assertEqual(self.seen(), set())
        second = self.run_discovery([record()], "second")
        self.assertEqual(len(self.read(second, "candidates.json")), 1)

    def test_bare_pending_context_is_retained_without_publication_acceptance(self):
        output = self.run_discovery([record(publications="Pending")])
        audit = self.read(output, "candidate_audit.json")[0]
        self.assertEqual(audit["publication_status"], "pending")
        marker = next(item for item in audit["publication_review"]["indicators"] if item["kind"] == "pending")
        self.assertEqual(marker["raw_value"], "Pending")
        self.assertIn("candidate.publications", marker["source_fields"])
        self.assertIn("candidate.metadata_variants[0].metadata.publications", marker["source_fields"])
        self.assertEqual(audit["papers"], [])
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.http_calls, [])
        self.assertEqual(self.seen(), set())

    def test_pending_in_an_unrelated_title_is_not_a_publication_indicator(self):
        source = record(publications="")
        source["title"] = "Pending"
        output = self.run_discovery([source])
        audit = self.read(output, "candidate_audit.json")[0]
        self.assertEqual(audit["publication_review"]["indicators"], [])
        self.assertEqual(audit["publication_status"], "unresolved")
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.http_calls, [])

    def test_pending_inside_a_publication_list_retains_its_exact_source_path(self):
        output = self.run_discovery([record(publications=[" Pending "])])
        audit = self.read(output, "candidate_audit.json")[0]
        marker = next(item for item in audit["publication_review"]["indicators"] if item["kind"] == "pending")
        self.assertEqual(marker["raw_value"], " Pending ")
        self.assertIn("candidate.publications[0]", marker["source_fields"])
        self.assertEqual(audit["publication_status"], "pending")
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.http_calls, [])

    def test_missing_history_is_rejected_before_search(self):
        api = self.api()
        with patch.object(api, "fetch_dataset_pages", side_effect=AssertionError("Search must not start")):
            with self.assertRaisesRegex(ValueError, "(?i)initializ|missing"):
                api.run_published_discovery(["XL-MS"], self.root / "run", history_dir=self.root / "missing")
        self.assertFalse((self.root / "missing").exists())
        self.assertFalse((self.root / "run").exists())

    def test_dataset_doi_and_release_metadata_do_not_establish_a_paper(self):
        source = record(publications="")
        source["doi"] = "10.6019/PXD000001"
        output = self.run_discovery([source])
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.http_calls, [])
        self.assertNotEqual(self.read(output, "candidate_audit.json")[0]["publication_status"], "published")

    def test_paper_metadata_is_present_in_json_and_csv(self):
        output = self.run_discovery([record()])
        row = self.read(output, "candidates.json")[0]
        self.assertEqual(row["publication_status"], "published")
        self.assertEqual((row["papers"][0]["doi"], row["papers"][0]["title"], row["papers"][0]["venue"]),
                         (DOI, "An example article", "Example Journal"))
        with (output / "candidates.csv").open(newline="", encoding="utf-8") as handle:
            csv_row = list(csv.DictReader(handle))[0]
        self.assertEqual(json.loads(csv_row["paper_dois_json"]), [DOI])
        self.assertEqual(json.loads(csv_row["paper_pmids_json"]), [])
        self.assertEqual(json.loads(csv_row["papers_json"]), row["papers"])

    def test_pmid_only_article_is_eligible(self):
        output = self.run_discovery([record(publications={"pubmedID": "1234"})],
                                    responses={"1234": pubmed_body(None)})
        papers = self.read(output, "candidates.json")[0]["papers"]
        self.assertEqual(papers[0]["pmid"], "1234")
        self.assertIsNone(papers[0]["doi"])

    def test_independent_doi_and_pmid_links_merge_only_after_identifier_confirmation(self):
        links = '<a href="https://doi.org/10.1234/example">DOI</a> <a href="https://pubmed.ncbi.nlm.nih.gov/1234/">PMID</a>'
        output = self.run_discovery([record(publications=links)])
        papers = self.read(output, "candidates.json")[0]["papers"]
        self.assertEqual(len(papers), 1)
        self.assertEqual((papers[0]["doi"], papers[0]["pmid"]), (DOI, "1234"))
        self.assertEqual(len(papers[0]["verification_evidence"]), 2)

    def test_conflicting_explicit_doi_pmid_pair_is_not_promoted(self):
        output = self.run_discovery([record(publications={"doi": "10.1234/other", "pubmedID": "1234"})])
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.read(output, "candidate_audit.json")[0]["publication_status"], "unresolved")
        self.assertEqual(self.seen(), set())

    def test_preprint_only_reference_is_excluded(self):
        payload = json.loads(crossref_body())
        payload["message"]["type"] = "posted-content"
        output = self.run_discovery([record()], responses={DOI: payload})
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.read(output, "candidate_audit.json")[0]["publication_status"], "ineligible")

    def test_dataset_doi_reference_is_not_treated_as_a_journal_article(self):
        payload = json.loads(crossref_body())
        payload["message"]["type"] = "dataset"
        output = self.run_discovery([record()], responses={DOI: payload})
        self.assertEqual(self.read(output, "candidates.json"), [])

    def test_future_article_is_not_published_at_the_requested_as_of_date(self):
        output = self.run_discovery([record()], responses={DOI: crossref_body(2028)})
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.seen(), set())

    def test_http_failures_remain_unresolved_and_do_not_consume_history(self):
        output = self.run_discovery([record()], responses={DOI: requests.exceptions.Timeout("Offline")},
                                    publication_request_limit=1)
        self.assertEqual(self.read(output, "candidates.json"), [])
        audit = self.read(output, "candidate_audit.json")[0]
        self.assertEqual(audit["publication_status"], "unresolved")
        self.assertEqual(self.seen(), set())

    def test_request_budget_preserves_unresolved_candidates_in_the_audit(self):
        rows = [record(), record("PXD000002", {"doi": "10.1234/other"})]
        output = self.run_discovery(rows, publication_request_limit=1)
        self.assertEqual([row["accession"] for row in self.read(output, "candidates.json")], ["PXD000001"])
        self.assertEqual(self.read(output, "candidate_audit.json")[1]["publication_status"], "unresolved")
        self.assertEqual(len(self.http_calls), 1)

    def test_raw_response_and_parsed_page_must_match(self):
        api = self.api()
        source = page("XL-MS", [record()])
        altered = copy.deepcopy(source.data)
        altered["datasets"][0][1] = "Changed title"
        changed = CatalogPage(source.request_url, source.retrieved_at_utc, source.raw_body, altered)
        with patch.object(api, "fetch_dataset_pages", return_value=[changed]):
            with self.assertRaisesRegex(ValueError, "(?i)raw|response"):
                api.run_published_discovery(["XL-MS"], self.root / "run", history_dir=self.history)
        self.assertEqual(self.seen(), set())

    def test_incomplete_search_pages_cannot_publish_a_complete_run(self):
        api = self.api()
        source = page("XL-MS", [record()])
        data = copy.deepcopy(source.data)
        data["result_set"]["n_available_rows"] = 2
        incomplete = CatalogPage(source.request_url, source.retrieved_at_utc, json.dumps(data).encode(), data)
        with patch.object(api, "fetch_dataset_pages", return_value=[incomplete]):
            with self.assertRaisesRegex(ValueError, "(?i)incomplete|count"):
                api.run_published_discovery(["XL-MS"], self.root / "run", history_dir=self.history)
        self.assertFalse((self.root / "run/manifest.json").exists())

    def test_registry_response_url_must_match_the_requested_search(self):
        api = self.api()
        source = page("XL-MS", [record()])
        wrong = CatalogPage(source.request_url.replace("proteomecentral.proteomexchange.org", "example.org"),
                            source.retrieved_at_utc, source.raw_body, source.data)
        with patch.object(api, "fetch_dataset_pages", return_value=[wrong]):
            with self.assertRaisesRegex(ValueError, "(?i)source|url|registry"):
                api.run_published_discovery(["XL-MS"], self.root / "run", history_dir=self.history)

    def test_empty_search_writes_csv_headers_and_a_complete_audit(self):
        output = self.run_discovery([])
        self.assertEqual(self.read(output, "candidates.json"), [])
        self.assertEqual(self.read(output, "candidate_audit.json"), [])
        with (output / "candidates.csv").open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            self.assertIn("papers_json", reader.fieldnames)
            self.assertEqual(list(reader), [])
        self.assertEqual(self.read(output, "manifest.json")["status"], "complete")
        self.assertEqual(self.http_calls, [])

    def test_different_datasets_with_the_same_paper_are_both_emitted(self):
        output = self.run_discovery([record(), record("PXD000002")])
        self.assertEqual(len(self.read(output, "candidates.json")), 2)
        self.assertEqual(len(self.http_calls), 1)

    def test_search_terms_are_trimmed_and_deduplicated(self):
        api = self.api()
        with patch.object(api, "fetch_dataset_pages", side_effect=self.fetcher([])):
            api.run_published_discovery([" XL-MS ", "XL-MS"], self.root / "run", history_dir=self.history)
        self.assertEqual(self.search_calls, ["XL-MS"])

    def test_original_source_rows_and_raw_publication_fields_are_retained(self):
        rows = [record()]
        before = copy.deepcopy(rows)
        output = self.run_discovery(rows)
        original = self.read(output, "source_candidates.json")[0]
        accepted = self.read(output, "candidates.json")[0]
        self.assertEqual(original["publications"], before[0]["publications"])
        self.assertEqual(accepted["publications"], before[0]["publications"])
        self.assertEqual(rows, before)
        self.assertEqual(accepted["publication_review"]["association_status"], "explicit_source_reference")

    def test_all_metadata_variants_are_checked_and_pending_evidence_is_retained(self):
        api = self.api()
        pages = {"first": [record(publications="Pending")], "second": [record()]}
        with patch.object(api, "fetch_dataset_pages", side_effect=lambda term, page_size: [page(term, pages[term], page_size)]), \
             patch.object(api, "PublicationLookup", side_effect=self.factory()):
            output = api.run_published_discovery(["first", "second"], self.root / "run",
                                                history_dir=self.history, as_of="2026-11-09")
        candidate = self.read(output, "candidates.json")[0]
        self.assertTrue(candidate["metadata_conflict"])
        self.assertTrue(candidate["publication_review"]["indicators"])

    def test_existing_output_is_rejected_before_search(self):
        api = self.api()
        output = self.root / "existing"
        output.mkdir()
        (output / "keep").write_bytes(b"keep")
        with patch.object(api, "fetch_dataset_pages", side_effect=AssertionError("Search must not start")):
            with self.assertRaises(FileExistsError):
                api.run_published_discovery(["XL-MS"], output, history_dir=self.history)
        self.assertEqual((output / "keep").read_bytes(), b"keep")

    def test_later_search_failure_does_not_consume_any_candidates(self):
        api = self.api()
        with patch.object(api, "fetch_dataset_pages", side_effect=[[page("first", [record()])],
                                                                 requests.exceptions.Timeout("Offline")]):
            with self.assertRaises(requests.exceptions.Timeout):
                api.run_published_discovery(["first", "second"], self.root / "run", history_dir=self.history)
        self.assertEqual(self.seen(), set())
        self.assertFalse((self.root / "run/manifest.json").exists())

    def test_csv_write_failure_does_not_consume_the_candidate(self):
        api = self.api()
        with patch.object(api, "_write_candidate_csv", side_effect=OSError("Disk failure")):
            with self.assertRaises(OSError):
                self.run_discovery([record()])
        self.assertEqual(self.seen(), set())
        self.assertFalse((self.root / "run/manifest.json").exists())

    def test_manifest_covers_sources_audit_csv_and_publication_evidence(self):
        output = self.run_discovery([record()])
        manifest = self.read(output, "manifest.json")
        self.assertTrue({"source_candidates.json", "candidate_audit.json", "candidates.csv"} <= set(manifest["output_sha256"]))
        self.assertTrue(any(name.startswith("publication_evidence/bodies/") for name in manifest["output_sha256"]))
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)

    def test_packaged_http_evidence_can_be_requalified_without_network(self):
        output = self.run_discovery([record()])
        evidence = self.read(output, "candidates.json")[0]["papers"][0]["verification_evidence"][0]
        result = replay_publication_lookup(output / evidence["root"], evidence["receipt_path"],
                                          receipt_sha256=evidence["receipt_sha256"], doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["paper"]["doi"], DOI)

    def test_cli_emits_papers_and_excludes_previously_emitted_candidates(self):
        api = self.api()
        from claudio.dataset_catalog.cli import main
        for name in ("first", "second"):
            with patch.object(api, "fetch_dataset_pages", side_effect=self.fetcher([record()])), \
                 patch.object(api, "PublicationLookup", side_effect=self.factory()):
                result = CliRunner().invoke(main, ["discover", "--search-term", "XL-MS",
                                                  "--history", str(self.history), "--output", str(self.root / name)])
            self.assertEqual(result.exit_code, 0, f"{result.output}\n{result.exception}")
        self.assertEqual(len(self.read(self.root / "first", "candidates.json")), 1)
        self.assertEqual(self.read(self.root / "second", "candidates.json"), [])
        self.assertEqual(len(self.http_calls), 1)

    def test_cli_requires_explicit_history_initialization_before_network(self):
        api = self.api()
        from claudio.dataset_catalog.cli import main
        with patch.object(api, "fetch_dataset_pages", side_effect=AssertionError("Search must not start")):
            result = CliRunner().invoke(main, ["discover", "--history", str(self.root / "missing"),
                                              "--output", str(self.root / "run")])
        self.assertEqual(result.exit_code, 1, result.output)
        self.assertIn("initializ", result.output.lower())
        self.assertFalse((self.root / "missing").exists())

    def test_cli_can_initialize_a_fresh_history_but_cannot_reset_it(self):
        from claudio.dataset_catalog.cli import main
        target = self.root / "fresh"
        first = CliRunner().invoke(main, ["init-discovery-history", "--history", str(target)])
        self.assertEqual(first.exit_code, 0, f"{first.output}\n{first.exception}")
        second = CliRunner().invoke(main, ["init-discovery-history", "--history", str(target)])
        self.assertNotEqual(second.exit_code, 0)
        with DiscoveryHistory(target) as history:
            self.assertEqual(history.seen_accessions(), set())

    def test_cli_rejects_invalid_publication_budget_before_search(self):
        from claudio.dataset_catalog.cli import main
        result = CliRunner().invoke(main, ["discover", "--history", str(self.history),
                                          "--publication-request-limit", "0", "--output", str(self.root / "run")])
        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("--publication-request-limit", result.output)


if __name__ == "__main__":
    unittest.main()
