"""Qualify supplied PubMed XML without making network requests."""

import re
from xml.etree import ElementTree as ET

from .publications import normalize_doi, normalize_pmid
from .publication_verification import _as_of_date, _date_interval


_PUBLISHED_STATUSES = {"ppublish", "epublish", "aheadofprint", "ecollection"}
_MONTH_NAMES = (
    "january", "february", "march", "april", "may", "june",
    "july", "august", "september", "october", "november", "december",
)
_MONTHS = {
    key: number
    for number, name in enumerate(_MONTH_NAMES, 1)
    for key in (name, name[:3])
}
_MONTHS["sept"] = 9


def _text(node):
    return "".join(node.itertext()).strip() if node is not None else ""


def _xml_date(node, source_field):
    year = _text(node.find("Year"))
    month = _text(node.find("Month"))
    day = _text(node.find("Day"))
    raw_medline = _text(node.find("MedlineDate"))
    if not year and raw_medline:
        years = set(re.findall(r"\b[0-9]{4}\b", raw_medline))
        if len(years) != 1:
            raise ValueError("MedlineDate has no unique year interval")
        year = next(iter(years))
    if not re.fullmatch(r"[0-9]{4}", year):
        raise ValueError("Publication year is missing or invalid")
    parts = [int(year)]
    if month:
        if re.fullmatch(r"[0-9]+", month):
            parts.append(int(month))
        elif month.lower().rstrip(".") in _MONTHS:
            parts.append(_MONTHS[month.lower().rstrip(".")])
        else:
            raise ValueError("Publication month is unresolved")
    if day:
        if not month or not re.fullmatch(r"[0-9]+", day):
            raise ValueError("Publication day is unresolved")
        parts.append(int(day))
    result = _date_interval({"date-parts": [parts]}, source_field)
    if raw_medline:
        result["raw_medline_date"] = raw_medline
    return result


def qualify_pubmed_publication(
    xml_body, *, requested_pmid, requested_doi=None, as_of
):
    """Interpret supplied PubMed XML; HTTP evidence must be verified separately."""
    pmid = normalize_pmid(requested_pmid)
    if pmid is None:
        raise ValueError("requested_pmid is invalid")
    doi_request = normalize_doi(requested_doi)
    if requested_doi is not None and doi_request is None:
        raise ValueError("requested_doi is invalid")
    cutoff = _as_of_date(as_of)

    def outcome(status, reason, paper=None):
        return {
            "provider": "pubmed",
            "requested_pmid": pmid,
            "requested_doi": doi_request,
            "verification_scope": "supplied_metadata",
            "status": status,
            "reason": reason,
            "paper": paper,
        }

    if not isinstance(xml_body, (bytes, str)):
        return outcome("unresolved", "invalid_xml_response")
    try:
        raw = xml_body.encode("utf-8") if isinstance(xml_body, str) else xml_body
        if len(raw) > 2 * 1024 * 1024:
            return outcome("unresolved", "xml_response_too_large")
        if re.search(br"<!ENTITY\b|<!DOCTYPE[^>]*\[", raw, re.IGNORECASE):
            return outcome("unresolved", "xml_entity_declaration_not_allowed")
        root = ET.fromstring(raw)
    except (ET.ParseError, UnicodeError, ValueError):
        return outcome("unresolved", "invalid_xml_response")
    if root.tag != "PubmedArticleSet":
        return outcome("unresolved", "unexpected_xml_root")

    matches = [
        item for item in root.findall("PubmedArticle")
        if normalize_pmid(_text(item.find("MedlineCitation/PMID"))) == pmid
    ]
    if len(matches) != 1:
        return outcome("unresolved", "requested_article_not_unique")
    selected = matches[0]
    item = selected.find("MedlineCitation/Article")
    data = selected.find("PubmedData")
    if item is None or data is None:
        return outcome("unresolved", "article_metadata_missing")

    identifier_nodes = data.findall("ArticleIdList/ArticleId")
    article_pmids = [
        normalize_pmid(_text(node)) for node in identifier_nodes
        if node.get("IdType") == "pubmed"
    ]
    if any(value != pmid for value in article_pmids):
        return outcome("unresolved", "conflicting_pubmed_identifiers")

    doi_values = {
        value
        for node in (
            [node for node in identifier_nodes if node.get("IdType") == "doi"]
            + [
                node for node in item.findall("ELocationID")
                if node.get("EIdType") == "doi" and node.get("ValidYN") != "N"
            ]
        )
        if (value := normalize_doi(_text(node))) is not None
    }
    if len(doi_values) > 1:
        return outcome("unresolved", "conflicting_doi_identifiers")
    doi = next(iter(doi_values), None)
    if doi_request is not None and doi != doi_request:
        return outcome("unresolved", "doi_pmid_equivalence_not_confirmed")

    publication_types = [
        _text(node) for node in item.findall("PublicationTypeList/PublicationType")
    ]
    if any(value.lower() in {"preprint", "preprints"} for value in publication_types):
        return outcome("ineligible", "preprint_not_formal_article")

    publication_status = _text(data.find("PublicationStatus")).lower()
    if publication_status in {"received", "accepted"}:
        return outcome("ineligible", "article_not_yet_published")
    if publication_status not in _PUBLISHED_STATUSES:
        return outcome("unresolved", "publication_status_not_confirmed")

    title = _text(item.find("ArticleTitle"))
    venue = _text(item.find("Journal/Title"))
    if not title or not venue:
        return outcome("unresolved", "missing_title_or_venue")

    dates = []
    try:
        for index, node in enumerate(item.findall("ArticleDate")):
            if node.get("DateType") in {"Electronic", "Print"}:
                dates.append(_xml_date(node, f"ArticleDate[{index}]"))
        node = item.find("Journal/JournalIssue/PubDate")
        if node is not None:
            dates.append(_xml_date(node, "Journal/JournalIssue/PubDate"))
        for index, node in enumerate(data.findall("History/PubMedPubDate")):
            if node.get("PubStatus") in {"epublish", "ppublish", "aheadofprint"}:
                dates.append(_xml_date(node, f"History/PubMedPubDate[{index}]"))
    except (ValueError, TypeError, OverflowError):
        return outcome("unresolved", "invalid_publication_date")
    if not dates:
        return outcome("unresolved", "publication_date_missing")

    confirmed = [value for value in dates if value["latest"] <= cutoff.isoformat()]
    if not confirmed:
        if all(value["earliest"] > cutoff.isoformat() for value in dates):
            return outcome("ineligible", "publication_date_in_future")
        return outcome("unresolved", "publication_date_precision_insufficient")

    authors = []
    for node in item.findall("AuthorList/Author"):
        collective = _text(node.find("CollectiveName"))
        if collective:
            authors.append({"name": collective})
        else:
            family = _text(node.find("LastName"))
            given = _text(node.find("ForeName")) or _text(node.find("Initials"))
            if family or given:
                authors.append({"given": given, "family": family})

    date_info = confirmed[0]
    paper = {
        "doi": doi,
        "pmid": pmid,
        "title": title,
        "authors": authors,
        "venue": venue,
        "year": date_info["date_parts"][0],
        "publication_date": date_info,
        "publication_dates": dates,
        "work_type": "journal-article",
        "publication_types": publication_types,
        "source_publication_status": publication_status,
        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    }
    return outcome("published", "published_article_metadata_confirmed", paper)
