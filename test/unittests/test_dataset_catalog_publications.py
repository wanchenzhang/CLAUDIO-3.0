"""Offline tests for explicit dataset-to-paper associations."""

import copy
import importlib
import unittest


class PublicationAPITestCase(unittest.TestCase):
    def api(self):
        try:
            return importlib.import_module(
                "claudio.dataset_catalog.publications"
            )
        except ModuleNotFoundError as error:
            if error.name == "claudio.dataset_catalog.publications":
                self.fail("Publication association parser is not implemented")
            raise


class TestPublicationIdentifiers(PublicationAPITestCase):
    def test_normalizes_doi_identifiers_and_resolver_urls(self):
        cases = [
            (" 10.1000/ABC ", "10.1000/abc"),
            ("doi:10.1000/ABC(2024)", "10.1000/abc(2024)"),
            (
                "HTTPS://DX.DOI.ORG/10.1016/J.MCPRO.2024.100724",
                "10.1016/j.mcpro.2024.100724",
            ),
            (
                "https://doi.org/10.1000%2FABC?utm_source=test",
                "10.1000/abc",
            ),
        ]
        for raw, expected in cases:
            with self.subTest(raw=raw):
                self.assertEqual(self.api().normalize_doi(raw), expected)

    def test_rejects_invalid_or_untrusted_doi_values(self):
        cases = [
            None, "", 123, True,
            "10.1000/", "10.1000/a b",
            "https://example.org/10.1000/article",
            "https://doi.org.example.org/10.1000/article",
        ]
        for raw in cases:
            with self.subTest(raw=raw):
                self.assertIsNone(self.api().normalize_doi(raw))

    def test_normalizes_positive_pmids_and_official_urls(self):
        cases = [
            (38266916, "38266916"),
            (" 38266916 ", "38266916"),
            ("00038266916", "38266916"),
            (
                "https://pubmed.ncbi.nlm.nih.gov/38266916/",
                "38266916",
            ),
            (
                "http://www.ncbi.nlm.nih.gov/pubmed/38266916",
                "38266916",
            ),
        ]
        for raw, expected in cases:
            with self.subTest(raw=raw):
                self.assertEqual(self.api().normalize_pmid(raw), expected)

    def test_rejects_zero_boolean_float_and_untrusted_pmids(self):
        cases = [
            None, "", 0, "0", -1, "-1", True, False,
            38266916.0, "38 266916", "1.5",
            "https://example.org/pubmed/38266916",
            "https://pubmed.ncbi.nlm.nih.gov.example.org/38266916/",
        ]
        for raw in cases:
            with self.subTest(raw=raw):
                self.assertIsNone(self.api().normalize_pmid(raw))


class TestPublicationLinks(PublicationAPITestCase):
    def extract(self, candidate, repository_metadata=None):
        return self.api().extract_publication_links(
            candidate, repository_metadata=repository_metadata
        )

    def ids(self, result):
        return {
            (reference["doi"], reference["pmid"])
            for reference in result["references"]
        }

    def test_html_extracts_article_ids_without_promoting_dataset_doi(self):
        candidate = {
            "accession": "PXD046392",
            "doi": "10.6019/PXD046392",
            "publications": (
                '<a href="https://dx.doi.org/'
                '10.1016/J.MCPRO.2024.100724">Article DOI</a>; '
                '<a href="https://www.ncbi.nlm.nih.gov/pubmed/'
                '38266916">Article citation</a>'
            ),
        }
        result = self.extract(candidate)
        self.assertEqual(result["dataset_accession"], "PXD046392")
        self.assertEqual(result["publication_status"], "unverified")
        self.assertEqual(self.ids(result), {
            ("10.1016/j.mcpro.2024.100724", None),
            (None, "38266916"),
        })
        self.assertNotIn(("10.6019/pxd046392", None), self.ids(result))
        for reference in result["references"]:
            self.assertIn("candidate.publications", reference["source_fields"])

    def test_pending_placeholder_does_not_create_a_paper(self):
        result = self.extract({
            "accession": "PXD000001",
            "publications": "Dataset with its publication pending",
        })
        self.assertEqual(result["references"], [])
        self.assertTrue(any(
            item["kind"] == "pending" for item in result["indicators"]
        ))
        self.assertEqual(result["publication_status"], "unverified")

    def test_no_manuscript_cv_term_is_not_a_pmid(self):
        result = self.extract({
            "accession": "PXD000001",
            "publications": [{
                "cvLabel": "MS",
                "accession": "MS:1002853",
                "name": "Dataset with no associated published manuscript",
            }],
        })
        self.assertEqual(result["references"], [])
        self.assertTrue(any(
            item["kind"] == "no_publication"
            for item in result["indicators"]
        ))

    def test_structured_reference_keeps_doi_when_pmid_is_zero(self):
        result = self.extract(
            {"accession": "PXD000001"},
            {
                "accession": "PXD000001",
                "doi": "10.6019/PXD000001",
                "references": [{
                    "doi": "10.1000/PAPER",
                    "pubmedID": 0,
                }],
            },
        )
        self.assertEqual(self.ids(result), {("10.1000/paper", None)})
        self.assertIn(
            "repository.references[0].doi",
            result["references"][0]["source_fields"],
        )

    def test_structured_reference_preserves_asserted_identifier_pair(self):
        result = self.extract(
            {"accession": "PXD000001"},
            {
                "accession": "PXD000001",
                "references": [{
                    "doi": "10.1000/PAPER",
                    "pubmedID": 38266916,
                }],
            },
        )
        self.assertEqual(
            self.ids(result), {("10.1000/paper", "38266916")}
        )
        self.assertEqual(result["publication_status"], "unverified")

    def test_proxi_pubmed_reference_requires_matching_dataset_binding(self):
        result = self.extract(
            {"accession": "PXD000001"},
            {
                "accession": [
                    {"accession": "MS:1001919", "value": "PXD000001"},
                    {"accession": "MS:1002487", "value": "MSV000000001"},
                ],
                "publications": [{
                    "accession": "MS:1000879",
                    "name": "PubMed identifier",
                    "value": "38266916",
                }],
            },
        )
        self.assertEqual(self.ids(result), {(None, "38266916")})
        self.assertIn(
            "repository.publications[0]",
            result["references"][0]["source_fields"],
        )

    def test_multiple_links_do_not_merge_identifiers_by_adjacency(self):
        result = self.extract({
            "accession": "PXD000001",
            "publications": (
                '<a href="https://doi.org/10.1000/first">Same title</a> '
                '<a href="https://pubmed.ncbi.nlm.nih.gov/38266916/">'
                'Same title</a> '
                '<a href="https://doi.org/10.1000/second">Second</a>'
            ),
        })
        self.assertEqual(self.ids(result), {
            ("10.1000/first", None),
            (None, "38266916"),
            ("10.1000/second", None),
        })
        self.assertEqual(len(result["references"]), 3)

    def test_duplicate_identifiers_retain_metadata_variant_sources(self):
        link = '<a href="https://doi.org/10.1000/PAPER">Paper</a>'
        result = self.extract({
            "accession": "PXD000001",
            "publications": link,
            "metadata_variants": [{
                "metadata": {"publications": link},
                "search_terms": ["crosslinking"],
            }],
        })
        self.assertEqual(self.ids(result), {("10.1000/paper", None)})
        self.assertEqual(len(result["references"]), 1)
        paths = result["references"][0]["source_fields"]
        self.assertIn("candidate.publications", paths)
        self.assertIn(
            "candidate.metadata_variants[0].metadata.publications", paths
        )

    def test_description_urls_and_release_date_do_not_create_a_paper(self):
        result = self.extract({
            "accession": "PXD000001",
            "title": "Study mentioning 10.1000/unlinked",
            "projectDescription":
                "See https://doi.org/10.1000/unlinked",
            "doi": "10.6019/PXD000001",
            "publicationDate": "2020-01-01",
            "status": "PUBLIC",
        })
        self.assertEqual(result["references"], [])
        self.assertEqual(result["publication_status"], "unverified")

    def test_preprint_doi_is_only_an_unverified_association(self):
        result = self.extract({
            "accession": "PXD000001",
            "publications":
                "https://doi.org/10.1101/2026.01.01.123456",
        })
        self.assertEqual(
            self.ids(result), {("10.1101/2026.01.01.123456", None)}
        )
        self.assertEqual(result["publication_status"], "unverified")

    def test_rejects_repository_metadata_for_another_dataset(self):
        with self.assertRaisesRegex(ValueError, "(?i)accession|dataset"):
            self.extract(
                {"accession": "PXD000001"},
                {
                    "accession": "PXD000002",
                    "references": [{"doi": "10.1000/unrelated"}],
                },
            )

    def test_rejects_repository_metadata_without_accession_binding(self):
        with self.assertRaisesRegex(ValueError, "(?i)accession|dataset"):
            self.extract(
                {"accession": "PXD000001"},
                {"references": [{"doi": "10.1000/unbound"}]},
            )

    def test_rejects_missing_or_invalid_candidate_accession(self):
        for accession in (None, "", "   ", 123):
            with self.subTest(accession=accession):
                with self.assertRaises(ValueError):
                    self.extract({"accession": accession})

    def test_preserves_original_candidate_and_repository_metadata(self):
        candidate = {
            "accession": "PXD000001",
            "publications": "Dataset with its publication pending",
            "review_status": "pending",
            "metadata_variants": [],
        }
        repository = {
            "accession": "PXD000001",
            "references": [{"doi": "10.1000/PAPER", "pubmedID": 0}],
        }
        original_candidate = copy.deepcopy(candidate)
        original_repository = copy.deepcopy(repository)
        result = self.extract(candidate, repository)
        self.assertEqual(self.ids(result), {("10.1000/paper", None)})
        self.assertEqual(candidate, original_candidate)
        self.assertEqual(repository, original_repository)


if __name__ == "__main__":
    unittest.main()
