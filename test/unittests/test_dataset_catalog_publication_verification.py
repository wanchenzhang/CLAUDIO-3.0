"""Offline Crossref qualification tests using synthetic metadata."""

import copy
import importlib
import unittest


AS_OF = "2026-11-09"
DOI = "10.1000/example.paper"


def record(parts=None, kind="journal-article"):
    parts = [2024, 5, 24] if parts is None else parts
    publication_date = {"date-parts": [list(parts)]}
    return {
        "status": "ok",
        "message-type": "work",
        "message-version": "1.0.0",
        "message": {
            "DOI": DOI,
            "type": kind,
            "title": ["Synthetic publication for offline testing"],
            "container-title": ["Example Research Journal"],
            "publisher": "Example Publisher",
            "source": "Crossref",
            "prefix": "10.1000",
            "member": "1",
            "reference-count": 0,
            "references-count": 0,
            "is-referenced-by-count": 0,
            "URL": "https://doi.org/" + DOI,
            "author": [
                {"given": "Ada", "family": "Lovelace"},
                {"given": "Alan", "family": "Turing"},
            ],
            "created": {
                "date-parts": [[2020, 1, 1]],
                "date-time": "2020-01-01T00:00:00Z",
                "timestamp": 1577836800000,
            },
            "deposited": {
                "date-parts": [[2024, 5, 24]],
                "date-time": "2024-05-24T00:00:00Z",
                "timestamp": 1716508800000,
            },
            "indexed": {
                "date-parts": [[2026, 11, 1]],
                "date-time": "2026-11-01T00:00:00Z",
                "timestamp": 1793491200000,
            },
            "published-online": copy.deepcopy(publication_date),
            "issued": copy.deepcopy(publication_date),
            "relation": {},
        },
    }


class TestCrossrefPublicationQualification(unittest.TestCase):
    def api(self):
        try:
            return importlib.import_module(
                "claudio.dataset_catalog.publication_verification"
            )
        except ModuleNotFoundError as error:
            if error.name == "claudio.dataset_catalog.publication_verification":
                self.fail("Crossref publication verifier is not implemented")
            raise

    def qualify(self, payload, requested_doi=DOI, as_of=AS_OF):
        return self.api().qualify_crossref_publication(
            payload, requested_doi=requested_doi, as_of=as_of
        )

    def test_confirms_published_journal_article_and_preserves_bibliography(self):
        result = self.qualify(record())
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["verification_scope"], "supplied_metadata")
        paper = result["paper"]
        self.assertEqual(paper["doi"], "10.1000/example.paper")
        self.assertIsNone(paper["pmid"])
        self.assertEqual(
            paper["title"], "Synthetic publication for offline testing"
        )
        self.assertEqual(paper["venue"], "Example Research Journal")
        self.assertEqual(paper["year"], 2024)
        self.assertEqual(paper["authors"], [
            {"given": "Ada", "family": "Lovelace"},
            {"given": "Alan", "family": "Turing"},
        ])
        self.assertEqual(paper["publication_date"]["date_parts"], [2024, 5, 24])
        self.assertEqual(paper["publication_date"]["precision"], "day")
        self.assertEqual(paper["url"], "https://doi.org/10.1000/example.paper")

    def test_accepts_formally_published_proceedings_article(self):
        result = self.qualify(record(kind="proceedings-article"))
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["paper"]["work_type"], "proceedings-article")

    def test_does_not_accept_preprint_dataset_or_book_as_article(self):
        for kind in ("posted-content", "dataset", "book"):
            with self.subTest(kind=kind):
                result = self.qualify(record(kind=kind))
                self.assertEqual(result["status"], "ineligible")
                self.assertIsNone(result["paper"])

    def test_requires_matching_requested_doi(self):
        payload = record()
        payload["message"]["DOI"] = "10.1000/unrelated.paper"
        result = self.qualify(payload)
        self.assertEqual(result["status"], "unresolved")
        self.assertIsNone(result["paper"])

    def test_accepts_case_and_resolver_variants_of_same_doi(self):
        result = self.qualify(
            record(),
            requested_doi="https://doi.org/10.1000/EXAMPLE.PAPER",
        )
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["paper"]["doi"], DOI)

    def test_api_error_or_wrong_response_shape_remains_unresolved(self):
        cases = [
            None,
            [],
            {"status": "failed", "message": "temporary failure"},
            {"status": "ok", "message-type": "work-list", "message": {}},
            {"status": "ok", "message-type": "work", "message": {}},
        ]
        for payload in cases:
            with self.subTest(payload=payload):
                result = self.qualify(payload)
                self.assertEqual(result["status"], "unresolved")
                self.assertIsNone(result["paper"])

    def test_requires_article_title_and_publication_venue(self):
        for field in ("title", "container-title"):
            with self.subTest(field=field):
                payload = record()
                payload["message"][field] = []
                result = self.qualify(payload)
                self.assertEqual(result["status"], "unresolved")
                self.assertIsNone(result["paper"])

    def test_registration_and_acceptance_dates_are_not_publication_dates(self):
        payload = record()
        del payload["message"]["issued"]
        del payload["message"]["published-online"]
        payload["message"]["accepted"] = {"date-parts": [[2020, 1, 1]]}
        result = self.qualify(payload)
        self.assertEqual(result["status"], "unresolved")
        self.assertIsNone(result["paper"])

    def test_future_article_is_not_yet_eligible(self):
        result = self.qualify(record(parts=[2027, 1, 1]))
        self.assertEqual(result["status"], "ineligible")
        self.assertIsNone(result["paper"])

    def test_past_online_publication_is_not_blocked_by_future_print_date(self):
        payload = record()
        payload["message"]["published-print"] = {
            "date-parts": [[2027, 1, 1]]
        }
        payload["message"]["issued"] = {"date-parts": [[2027, 1, 1]]}
        result = self.qualify(payload)
        self.assertEqual(result["status"], "published")
        self.assertEqual(
            result["paper"]["publication_date"]["source_field"],
            "published-online",
        )

    def test_partial_dates_preserve_precision_and_defer_ambiguous_dates(self):
        cases = [
            ([2024], "published", "year"),
            ([2026, 10], "published", "month"),
            ([2026], "unresolved", None),
            ([2026, 11], "unresolved", None),
        ]
        for parts, status, precision in cases:
            with self.subTest(parts=parts):
                result = self.qualify(record(parts=parts))
                self.assertEqual(result["status"], status)
                if status == "published":
                    date_info = result["paper"]["publication_date"]
                    self.assertEqual(date_info["precision"], precision)
                    self.assertEqual(date_info["date_parts"], parts)
                else:
                    self.assertIsNone(result["paper"])

    def test_invalid_calendar_dates_do_not_become_published(self):
        for parts in ([2024, 2, 30], [2024, 13, 1], [True], []):
            with self.subTest(parts=parts):
                result = self.qualify(record(parts=parts))
                self.assertEqual(result["status"], "unresolved")
                self.assertIsNone(result["paper"])

    def test_publication_on_as_of_date_is_inclusive(self):
        result = self.qualify(record(parts=[2026, 11, 9]))
        self.assertEqual(result["status"], "published")

    def test_later_publication_can_change_eligibility(self):
        payload = record(parts=[2026, 11, 3])
        earlier = self.qualify(payload, as_of="2026-10-31")
        later = self.qualify(payload, as_of="2026-11-09")
        self.assertEqual(earlier["status"], "ineligible")
        self.assertEqual(later["status"], "published")

    def test_preserves_original_response(self):
        payload = record()
        original = copy.deepcopy(payload)
        self.qualify(payload)
        self.assertEqual(payload, original)

    def test_rejects_invalid_requested_identifier_or_as_of_date(self):
        with self.assertRaises(ValueError):
            self.qualify(record(), requested_doi="not-a-doi")
        with self.assertRaises(ValueError):
            self.qualify(record(), as_of="not-a-date")


if __name__ == "__main__":
    unittest.main()
