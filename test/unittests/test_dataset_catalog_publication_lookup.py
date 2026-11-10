
"""Offline contracts for bounded, replayable publication lookups."""

import hashlib
import importlib
import json
from pathlib import Path
import tempfile
import unittest

import requests

DOI = "10.1234/example"
PMID = "1234"


def crossref_body(year=2025):
    return json.dumps({
        "status": "ok", "message-type": "work",
        "message": {
            "DOI": DOI, "type": "journal-article",
            "title": ["An example article"],
            "container-title": ["Example Journal"], "author": [],
            "published-online": {"date-parts": [[year, 1, 1]]},
        },
    }).encode("utf-8")


def pubmed_body(doi=DOI):
    doi_element = "" if doi is None else '<ArticleId IdType="doi">' + doi + '</ArticleId>'
    return ("""<PubmedArticleSet><PubmedArticle><MedlineCitation>
      <PMID>1234</PMID><Article><ArticleTitle>An example article</ArticleTitle>
      <Journal><Title>Example Journal</Title><JournalIssue><PubDate>
      <Year>2025</Year><Month>Jan</Month><Day>1</Day>
      </PubDate></JournalIssue></Journal>
      <PublicationTypeList><PublicationType>Journal Article</PublicationType></PublicationTypeList>
      </Article></MedlineCitation><PubmedData><PublicationStatus>epublish</PublicationStatus>
      <ArticleIdList><ArticleId IdType="pubmed">1234</ArticleId>""" + doi_element +
      "</ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>").encode("utf-8")


class Clock:
    def __init__(self):
        self.value = 0.0
        self.waits = []

    def __call__(self):
        return self.value

    def sleep(self, duration):
        self.waits.append(duration)
        self.value += duration


class Response:
    def __init__(self, body=b"", status=200, headers=None, url=None, stream_error=False):
        self.body = body
        self.status_code = status
        self.headers = headers or {"Content-Type": "application/json"}
        self.url = url
        self.stream_error = stream_error
        self.closed = False

    def iter_content(self, chunk_size):
        for start in range(0, len(self.body), chunk_size):
            yield self.body[start:start + chunk_size]
        if self.stream_error:
            raise requests.exceptions.ConnectionError("Interrupted response")

    def close(self):
        self.closed = True


class Session:
    """Replace only external HTTP; real capture and replay run."""
    def __init__(self, *events):
        self.events = list(events)
        self.closed = False

    def get(self, url, *, timeout, allow_redirects, stream, headers):
        if allow_redirects is not False or stream is not True:
            raise AssertionError("HTTP capture must stream without following redirects")
        if timeout != (10, 30):
            raise AssertionError("HTTP timeouts must be bounded")
        if not headers.get("User-Agent"):
            raise AssertionError("HTTP requests must identify the client")
        if not self.events:
            raise AssertionError("Unexpected extra HTTP attempt")
        event = self.events.pop(0)
        if isinstance(event, Exception):
            raise event
        if event.url is None:
            event.url = url
        return event

    def close(self):
        self.closed = True


class TestPublicationLookup(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.clock = Clock()

    def api(self):
        try:
            return importlib.import_module("claudio.dataset_catalog.publication_lookup")
        except ModuleNotFoundError as error:
            if error.name != "claudio.dataset_catalog.publication_lookup":
                raise
            self.fail("Publication HTTP evidence client is not implemented")

    def client(self, session, **options):
        client = self.api().PublicationLookup(
            self.root / "evidence", session=session,
            clock=self.clock, sleeper=self.clock.sleep, **options,
        )
        self.addCleanup(client.close)
        return client

    def receipt(self, result):
        evidence = result["evidence"]
        raw = (self.root / "evidence" / evidence["receipt_path"]).read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), evidence["receipt_sha256"])
        return json.loads(raw)

    def replay(self, result, **options):
        return self.api().replay_publication_lookup(
            self.root / "evidence", result["evidence"]["receipt_path"],
            receipt_sha256=result["evidence"]["receipt_sha256"], **options,
        )

    def test_crossref_lookup_saves_reloadable_raw_evidence(self):
        body = crossref_body()
        response = Response(body)
        client = self.client(Session(response))
        result = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["paper"]["doi"], DOI)
        self.assertEqual(result["verification_scope"], "captured_http_metadata")
        self.assertIs(result["dataset_association_verified"], False)
        receipt = self.receipt(result)
        self.assertEqual(receipt["request_url"], "https://api.crossref.org/works/10.1234%2Fexample")
        attempt = receipt["attempts"][0]
        self.assertEqual((self.root / "evidence" / attempt["body_path"]).read_bytes(), body)
        self.assertEqual(attempt["body_sha256"], hashlib.sha256(body).hexdigest())
        self.assertTrue(attempt["complete"])
        self.assertTrue(response.closed)
        self.assertEqual(self.replay(result, doi=DOI, as_of="2026-11-09"), result)

    def test_pubmed_lookup_checks_requested_doi_and_pmid(self):
        result = self.client(Session(Response(pubmed_body()))).lookup(doi=DOI, pmid=PMID, as_of="2026-11-09")
        self.assertEqual(result["status"], "published")
        self.assertEqual((result["paper"]["doi"], result["paper"]["pmid"]), (DOI, PMID))
        self.assertEqual(self.receipt(result)["request_url"],
                         "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=1234&retmode=xml&rettype=abstract&tool=claudio")

    def test_pmid_only_published_article_does_not_require_a_doi(self):
        result = self.client(Session(Response(pubmed_body(None)))).lookup(pmid=PMID, as_of="2026-11-09")
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["paper"]["pmid"], PMID)
        self.assertIsNone(result["paper"]["doi"])

    def test_cached_raw_response_is_requalified_for_each_as_of_date(self):
        client = self.client(Session(Response(crossref_body())))
        earlier = client.lookup(doi=DOI, as_of="2024-12-31")
        later = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(earlier["status"], "ineligible")
        self.assertEqual(later["status"], "published")
        self.assertEqual(earlier["evidence"], later["evidence"])
        self.assertEqual(client.request_count, 1)

    def test_cached_pmid_does_not_confirm_a_different_doi(self):
        client = self.client(Session(Response(pubmed_body())))
        first = client.lookup(pmid=PMID, as_of="2026-11-09")
        second = client.lookup(pmid=PMID, doi="10.1234/other", as_of="2026-11-09")
        self.assertEqual(first["status"], "published")
        self.assertEqual(second["status"], "unresolved")
        self.assertIsNone(second["paper"])
        self.assertEqual(client.request_count, 1)

    def test_timeout_retry_is_saved_as_a_separate_http_attempt(self):
        client = self.client(Session(requests.exceptions.Timeout("Offline"), Response(crossref_body())))
        result = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "published")
        attempts = self.receipt(result)["attempts"]
        self.assertEqual(len(attempts), 2)
        self.assertEqual(attempts[0]["error"], "request_failed")
        self.assertEqual(attempts[0]["error_type"], "Timeout")
        self.assertEqual(attempts[0]["body_bytes"], 0)
        self.assertEqual(client.request_count, 2)

    def test_request_budget_counts_failed_attempts_and_retries(self):
        client = self.client(Session(Response(b"busy", 503), Response(b"busy", 503)), max_requests=2)
        result = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "request_budget_exhausted")
        self.assertEqual(len(self.receipt(result)["attempts"]), 2)
        again = client.lookup(doi="10.1234/other", as_of="2026-11-09")
        self.assertEqual(again["status"], "unresolved")
        self.assertIsNone(again["evidence"])
        self.assertEqual(client.request_count, 2)

    def test_response_byte_limit_retains_only_a_bounded_prefix(self):
        response = Response(b"x" * 100)
        result = self.client(Session(response), max_response_bytes=20).lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "response_byte_limit")
        attempt = self.receipt(result)["attempts"][0]
        self.assertFalse(attempt["complete"])
        self.assertEqual((self.root / "evidence" / attempt["body_path"]).read_bytes(), b"x" * 20)
        self.assertTrue(response.closed)

    def test_total_byte_limit_applies_across_distinct_lookups(self):
        body = crossref_body()
        client = self.client(Session(Response(body), Response(b"x" * 100)), max_total_bytes=len(body) + 5)
        self.assertEqual(client.lookup(doi=DOI, as_of="2026-11-09")["status"], "published")
        result = client.lookup(doi="10.1234/other", as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "total_byte_limit")
        self.assertEqual(self.receipt(result)["attempts"][0]["body_bytes"], 5)
        self.assertEqual(client.captured_body_bytes, len(body) + 5)

    def test_redirects_and_missing_records_are_unresolved_not_unpublished(self):
        for status, reason in ((302, "redirect_not_followed"), (404, "http_status_error")):
            with self.subTest(status=status), tempfile.TemporaryDirectory() as temporary:
                client = self.api().PublicationLookup(
                    Path(temporary) / "evidence", session=Session(Response(b"not an article", status)),
                    clock=self.clock, sleeper=self.clock.sleep,
                )
                self.addCleanup(client.close)
                result = client.lookup(doi=DOI, as_of="2026-11-09")
                self.assertEqual(result["status"], "unresolved")
                self.assertEqual(result["reason"], reason)
                self.assertEqual(client.request_count, 1)

    def test_unexpected_response_origin_or_identifier_cannot_qualify(self):
        for url in ("https://example.org/work", "https://api.crossref.org/works/10.1234%2Fother"):
            with self.subTest(url=url), tempfile.TemporaryDirectory() as temporary:
                client = self.api().PublicationLookup(
                    Path(temporary) / "evidence", session=Session(Response(crossref_body(), url=url)),
                    clock=self.clock, sleeper=self.clock.sleep,
                )
                self.addCleanup(client.close)
                result = client.lookup(doi=DOI, as_of="2026-11-09")
                self.assertEqual(result["status"], "unresolved")
                self.assertEqual(result["reason"], "unexpected_response_url")
                self.assertIsNone(result["paper"])

    def test_retry_after_is_honored_and_all_responses_are_retained(self):
        response = Response(b"slow down", 429, {"Retry-After": "2"})
        client = self.client(Session(response, Response(crossref_body())))
        result = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "published")
        self.assertIn(2.0, self.clock.waits)
        attempts = self.receipt(result)["attempts"]
        self.assertEqual([row["status_code"] for row in attempts], [429, 200])
        self.assertEqual(attempts[0]["headers"]["retry-after"], "2")
        self.assertTrue(response.closed)

    def test_long_retry_after_stops_without_violating_the_server_delay(self):
        client = self.client(Session(Response(b"busy", 429, {"Retry-After": "120"})))
        result = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "retry_wait_exceeds_limit")
        self.assertEqual(client.request_count, 1)
        self.assertEqual(self.clock.waits, [])

    def test_retry_exhaustion_cannot_be_mistaken_for_an_ineligible_article(self):
        client = self.client(Session(Response(b"busy", 503), Response(b"busy", 503)), max_attempts=2)
        result = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "http_status_error")
        self.assertEqual(len(self.receipt(result)["attempts"]), 2)

    def test_distinct_requests_are_paced(self):
        client = self.client(Session(Response(crossref_body()), Response(crossref_body())))
        client.lookup(doi=DOI, as_of="2026-11-09")
        client.lookup(doi="10.1234/other", as_of="2026-11-09")
        self.assertGreaterEqual(sum(self.clock.waits), 0.4)
        self.assertEqual(client.request_count, 2)

    def test_elapsed_time_limit_prevents_another_http_request(self):
        client = self.client(Session(), max_elapsed_seconds=1)
        self.clock.value = 2
        result = client.lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "time_budget_exhausted")
        self.assertEqual(client.request_count, 0)

    def test_invalid_json_is_saved_but_never_promoted_to_published(self):
        result = self.client(Session(Response(b"<html>maintenance</html>"))).lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        self.assertEqual(result["reason"], "invalid_crossref_json")
        self.assertIsNone(result["paper"])
        self.assertTrue(self.receipt(result)["attempts"][0]["complete"])

    def test_interrupted_stream_preserves_partial_bytes_and_cannot_qualify(self):
        response = Response(b"partial", stream_error=True)
        result = self.client(Session(response), max_attempts=1).lookup(doi=DOI, as_of="2026-11-09")
        self.assertEqual(result["status"], "unresolved")
        attempt = self.receipt(result)["attempts"][0]
        self.assertFalse(attempt["complete"])
        self.assertEqual(attempt["body_bytes"], 7)
        self.assertEqual(attempt["error"], "request_failed")
        self.assertTrue(response.closed)

    def test_replay_rejects_changed_body_bytes(self):
        result = self.client(Session(Response(crossref_body()))).lookup(doi=DOI, as_of="2026-11-09")
        attempt = self.receipt(result)["attempts"][0]
        (self.root / "evidence" / attempt["body_path"]).write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "(?i)checksum|body"):
            self.replay(result, doi=DOI, as_of="2026-11-09")

    def test_replay_rejects_changed_receipt_bytes(self):
        result = self.client(Session(Response(crossref_body()))).lookup(doi=DOI, as_of="2026-11-09")
        path = self.root / "evidence" / result["evidence"]["receipt_path"]
        path.write_bytes(path.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum|receipt"):
            self.replay(result, doi=DOI, as_of="2026-11-09")

    def test_replay_rejects_a_different_requested_identifier(self):
        result = self.client(Session(Response(crossref_body()))).lookup(doi=DOI, as_of="2026-11-09")
        with self.assertRaisesRegex(ValueError, "(?i)identifier|request"):
            self.replay(result, doi="10.1234/other", as_of="2026-11-09")

    def test_replay_rejects_body_paths_outside_the_evidence_directory(self):
        result = self.client(Session(Response(crossref_body()))).lookup(doi=DOI, as_of="2026-11-09")
        path = self.root / "evidence" / result["evidence"]["receipt_path"]
        receipt = self.receipt(result)
        (self.root / "outside.json").write_bytes(crossref_body())
        receipt["attempts"][0]["body_path"] = "../outside.json"
        raw = json.dumps(receipt).encode()
        path.write_bytes(raw)
        with self.assertRaisesRegex(ValueError, "(?i)path|directory"):
            self.api().replay_publication_lookup(
                self.root / "evidence", result["evidence"]["receipt_path"],
                receipt_sha256=hashlib.sha256(raw).hexdigest(), doi=DOI, as_of="2026-11-09",
            )

    def test_invalid_identifiers_and_dates_fail_before_http(self):
        client = self.client(Session())
        for options in ({"doi": "bad"}, {"pmid": "0"}, {}, {"doi": DOI, "as_of": "not-a-date"}):
            with self.subTest(options=options):
                with self.assertRaises(ValueError):
                    client.lookup(**{"as_of": "2026-11-09", **options})
        self.assertEqual(client.request_count, 0)

    def test_invalid_limits_do_not_create_evidence_directories(self):
        for options in ({"max_requests": 0}, {"max_response_bytes": 0},
                        {"max_attempts": 0}, {"max_elapsed_seconds": float("nan")}):
            with self.subTest(options=options):
                target = self.root / "invalid"
                with self.assertRaises(ValueError):
                    self.api().PublicationLookup(target, session=Session(), **options)
                self.assertFalse(target.exists())

    def test_existing_evidence_directory_is_never_overwritten(self):
        target = self.root / "keep"
        target.mkdir()
        marker = target / "original.txt"
        marker.write_bytes(b"original")
        with self.assertRaises(FileExistsError):
            self.api().PublicationLookup(target, session=Session())
        self.assertEqual(marker.read_bytes(), b"original")


if __name__ == "__main__":
    unittest.main()
