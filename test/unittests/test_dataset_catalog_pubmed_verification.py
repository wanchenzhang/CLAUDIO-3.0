"""Offline PubMed qualification tests using synthetic XML."""

import copy
import importlib
import unittest
from xml.etree import ElementTree as ET


PMID = "38266916"
DOI = "10.1000/example.paper"
AS_OF = "2026-11-09"

XML = """<PubmedArticleSet>
<PubmedArticle>
 <MedlineCitation Status="MEDLINE" Owner="NLM">
  <PMID Version="1">38266916</PMID>
  <Article PubModel="Electronic-Print">
   <Journal>
    <ISSN IssnType="Electronic">1234-5678</ISSN>
    <JournalIssue CitedMedium="Internet">
     <Volume>1</Volume><Issue>2</Issue>
     <PubDate><Year>2024</Year><Month>May</Month><Day>24</Day></PubDate>
    </JournalIssue>
    <Title>Example Research Journal</Title>
    <ISOAbbreviation>Example Res J</ISOAbbreviation>
   </Journal>
   <ArticleTitle>Synthetic <i>publication</i> for offline testing</ArticleTitle>
   <ELocationID EIdType="doi" ValidYN="Y">10.1000/example.paper</ELocationID>
   <AuthorList CompleteYN="Y">
    <Author ValidYN="Y"><LastName>Lovelace</LastName><ForeName>Ada</ForeName></Author>
    <Author ValidYN="Y"><CollectiveName>Example Consortium</CollectiveName></Author>
   </AuthorList>
   <Language>eng</Language>
   <PublicationTypeList>
    <PublicationType UI="D016428">Journal Article</PublicationType>
   </PublicationTypeList>
   <ArticleDate DateType="Electronic">
    <Year>2024</Year><Month>05</Month><Day>24</Day>
   </ArticleDate>
  </Article>
 </MedlineCitation>
 <PubmedData>
  <History>
   <PubMedPubDate PubStatus="received">
    <Year>2020</Year><Month>1</Month><Day>1</Day>
   </PubMedPubDate>
   <PubMedPubDate PubStatus="pubmed">
    <Year>2026</Year><Month>11</Month><Day>1</Day>
   </PubMedPubDate>
  </History>
  <PublicationStatus>ppublish</PublicationStatus>
  <ArticleIdList>
   <ArticleId IdType="pubmed">38266916</ArticleId>
   <ArticleId IdType="doi">10.1000/example.paper</ArticleId>
  </ArticleIdList>
 </PubmedData>
</PubmedArticle>
</PubmedArticleSet>"""

DOCTYPE = (
    '<?xml version="1.0" encoding="UTF-8"?>\n'
    '<!DOCTYPE PubmedArticleSet PUBLIC '
    '"-//NLM//DTD PubMedArticle, 1st January 2026//EN" '
    '"https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_260101.dtd">\n'
)


def fixture():
    return ET.fromstring(XML)


def body(root):
    return DOCTYPE.encode("utf-8") + ET.tostring(root, encoding="utf-8")


def article(root):
    return root.find("PubmedArticle/MedlineCitation/Article")


def replace_date(node, parts):
    for child in list(node):
        node.remove(child)
    for name, value in zip(("Year", "Month", "Day"), parts):
        ET.SubElement(node, name).text = str(value)


def remove_dois(root):
    item = article(root)
    for node in list(item.findall("ELocationID")):
        if node.get("EIdType") == "doi":
            item.remove(node)
    ids = root.find("PubmedArticle/PubmedData/ArticleIdList")
    for node in list(ids):
        if node.get("IdType") == "doi":
            ids.remove(node)


class TestPubMedPublicationQualification(unittest.TestCase):
    def api(self):
        try:
            return importlib.import_module(
                "claudio.dataset_catalog.pubmed_verification"
            )
        except ModuleNotFoundError as error:
            if error.name == "claudio.dataset_catalog.pubmed_verification":
                self.fail("PubMed publication verifier is not implemented")
            raise

    def qualify(self, root, **options):
        payload = body(root) if isinstance(root, ET.Element) else root
        return self.api().qualify_pubmed_publication(
            payload,
            requested_pmid=options.get("requested_pmid", PMID),
            requested_doi=options.get("requested_doi"),
            as_of=options.get("as_of", AS_OF),
        )

    def test_accepts_published_pmid_only_article_without_doi(self):
        root = fixture()
        remove_dois(root)
        result = self.qualify(root)
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["verification_scope"], "supplied_metadata")
        self.assertEqual(result["paper"]["pmid"], "38266916")
        self.assertIsNone(result["paper"]["doi"])

    def test_preserves_title_authors_venue_and_date(self):
        result = self.qualify(fixture())
        paper = result["paper"]
        self.assertEqual(result["status"], "published")
        self.assertEqual(
            paper["title"], "Synthetic publication for offline testing"
        )
        self.assertEqual(paper["venue"], "Example Research Journal")
        self.assertEqual(paper["year"], 2024)
        self.assertEqual(paper["authors"], [
            {"given": "Ada", "family": "Lovelace"},
            {"name": "Example Consortium"},
        ])
        self.assertEqual(paper["publication_date"]["date_parts"], [2024, 5, 24])
        self.assertEqual(paper["url"], "https://pubmed.ncbi.nlm.nih.gov/38266916/")

    def test_confirms_explicit_doi_pmid_pair(self):
        result = self.qualify(
            fixture(), requested_doi="https://doi.org/10.1000/EXAMPLE.PAPER"
        )
        self.assertEqual(result["status"], "published")
        self.assertEqual(result["paper"]["doi"], DOI)

    def test_does_not_claim_unconfirmed_doi_pmid_equivalence(self):
        for missing in (False, True):
            root = fixture()
            if missing:
                remove_dois(root)
            result = self.qualify(root, requested_doi="10.1000/other.paper")
            self.assertEqual(result["status"], "unresolved")
            self.assertIsNone(result["paper"])

    def test_rejects_different_requested_pmid(self):
        result = self.qualify(fixture(), requested_pmid="12345678")
        self.assertEqual(result["status"], "unresolved")

    def test_rejects_conflicting_pubmed_identifiers(self):
        root = fixture()
        root.find(
            "PubmedArticle/PubmedData/ArticleIdList/ArticleId[@IdType='pubmed']"
        ).text = "12345678"
        self.assertEqual(self.qualify(root)["status"], "unresolved")

    def test_rejects_conflicting_doi_metadata(self):
        root = fixture()
        article(root).find("ELocationID").text = "10.1000/other.paper"
        self.assertEqual(self.qualify(root)["status"], "unresolved")

    def test_preprint_is_not_promoted_by_published_status(self):
        root = fixture()
        article(root).find("PublicationTypeList/PublicationType").text = "Preprint"
        result = self.qualify(root)
        self.assertEqual(result["status"], "ineligible")
        self.assertIsNone(result["paper"])

    def test_received_or_accepted_status_is_not_formal_publication(self):
        for status in ("received", "accepted"):
            with self.subTest(status=status):
                root = fixture()
                root.find("PubmedArticle/PubmedData/PublicationStatus").text = status
                self.assertEqual(self.qualify(root)["status"], "ineligible")

    def test_ahead_of_print_is_distinct_from_preprint(self):
        root = fixture()
        root.find("PubmedArticle/PubmedData/PublicationStatus").text = "aheadofprint"
        self.assertEqual(self.qualify(root)["status"], "published")

    def test_indexing_and_received_dates_are_not_publication_dates(self):
        root = fixture()
        item = article(root)
        item.remove(item.find("ArticleDate"))
        issue = item.find("Journal/JournalIssue")
        issue.remove(issue.find("PubDate"))
        self.assertEqual(self.qualify(root)["status"], "unresolved")

    def test_past_electronic_publication_survives_future_print_date(self):
        root = fixture()
        replace_date(
            article(root).find("Journal/JournalIssue/PubDate"), [2027, 1, 1]
        )
        self.assertEqual(self.qualify(root)["status"], "published")

    def test_future_and_partial_dates_are_not_guessed(self):
        cases = [
            ([2027, 1, 1], "ineligible", None),
            ([2026], "unresolved", None),
            ([2024], "published", "year"),
            ([2026, 10], "published", "month"),
            ([2026, 11, 9], "published", "day"),
        ]
        for parts, status, precision in cases:
            with self.subTest(parts=parts):
                root = fixture()
                item = article(root)
                replace_date(item.find("ArticleDate"), parts)
                replace_date(item.find("Journal/JournalIssue/PubDate"), parts)
                result = self.qualify(root)
                self.assertEqual(result["status"], status)
                if precision:
                    self.assertEqual(
                        result["paper"]["publication_date"]["precision"],
                        precision,
                    )

    def test_invalid_calendar_date_is_unresolved(self):
        root = fixture()
        replace_date(article(root).find("ArticleDate"), [2024, 2, 30])
        self.assertEqual(self.qualify(root)["status"], "unresolved")

    def test_invalid_xml_or_entity_declarations_are_unresolved(self):
        cases = [
            None,
            b"<broken",
            b"<html><body>Service unavailable</body></html>",
            (
                b'<!DOCTYPE PubmedArticleSet [<!ENTITY x "unexpected">]>'
                b"<PubmedArticleSet>&x;</PubmedArticleSet>"
            ),
        ]
        for payload in cases:
            with self.subTest(payload=payload):
                self.assertEqual(self.qualify(payload)["status"], "unresolved")

    def test_selects_only_the_requested_article_in_a_batch(self):
        root = fixture()
        other = copy.deepcopy(root.find("PubmedArticle"))
        other.find("MedlineCitation/PMID").text = "12345678"
        other.find(
            "PubmedData/ArticleIdList/ArticleId[@IdType='pubmed']"
        ).text = "12345678"
        other.find("MedlineCitation/Article/ArticleTitle").text = "Unrelated"
        root.append(other)
        result = self.qualify(root)
        self.assertEqual(result["status"], "published")
        self.assertEqual(
            result["paper"]["title"], "Synthetic publication for offline testing"
        )

    def test_duplicate_requested_articles_are_ambiguous(self):
        root = fixture()
        root.append(copy.deepcopy(root.find("PubmedArticle")))
        self.assertEqual(self.qualify(root)["status"], "unresolved")

    def test_rejects_invalid_request_arguments(self):
        with self.assertRaises(ValueError):
            self.qualify(fixture(), requested_pmid="0")
        with self.assertRaises(ValueError):
            self.qualify(fixture(), requested_doi="not-a-doi")
        with self.assertRaises(ValueError):
            self.qualify(fixture(), as_of="not-a-date")


if __name__ == "__main__":
    unittest.main()
