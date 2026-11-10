"""Qualify supplied bibliographic metadata without making HTTP requests."""

import calendar
import copy
import re
from datetime import date

from .publications import normalize_doi


_ARTICLE_TYPES = {"journal-article", "proceedings-article"}
_PUBLICATION_DATE_FIELDS = (
    "published-online", "published-print", "published", "issued",
)


def _as_of_date(value):
    if type(value) is date:
        return value
    if not isinstance(value, str) or re.fullmatch(
        r"[0-9]{4}-[0-9]{2}-[0-9]{2}", value
    ) is None:
        raise ValueError("as_of must be an ISO calendar date")
    return date.fromisoformat(value)


def _date_interval(value, source_field):
    if not isinstance(value, dict):
        raise ValueError("Publication date must contain date-parts")
    arrays = value.get("date-parts")
    if not isinstance(arrays, list) or len(arrays) != 1:
        raise ValueError("Publication date has no unique date-parts array")
    parts = arrays[0]
    if (
        not isinstance(parts, list)
        or not 1 <= len(parts) <= 3
        or any(type(part) is not int for part in parts)
    ):
        raise ValueError("Invalid publication date parts")
    year = parts[0]
    month = parts[1] if len(parts) >= 2 else 1
    day = parts[2] if len(parts) == 3 else 1
    earliest = date(year, month, day)
    if len(parts) == 1:
        latest = date(year, 12, 31)
    elif len(parts) == 2:
        latest = date(year, month, calendar.monthrange(year, month)[1])
    else:
        latest = earliest
    return {
        "date_parts": list(parts),
        "precision": {1: "year", 2: "month", 3: "day"}[len(parts)],
        "earliest": earliest.isoformat(),
        "latest": latest.isoformat(),
        "source_field": source_field,
    }


def _first_text(value):
    if not isinstance(value, list):
        return None
    return next(
        (item.strip() for item in value if isinstance(item, str) and item.strip()),
        None,
    )


def qualify_crossref_publication(payload, *, requested_doi, as_of):
    """Interpret a supplied Crossref work response as of a fixed date.

    The caller must separately validate the HTTP origin, status and evidence
    checksum. This pure function does not establish a dataset association.
    """
    requested = normalize_doi(requested_doi)
    if requested is None:
        raise ValueError("requested_doi is invalid")
    cutoff = _as_of_date(as_of)

    def outcome(status, reason, paper=None):
        return {
            "provider": "crossref",
            "requested_doi": requested,
            "verification_scope": "supplied_metadata",
            "status": status,
            "reason": reason,
            "paper": paper,
        }

    if (
        not isinstance(payload, dict)
        or payload.get("status") != "ok"
        or payload.get("message-type") != "work"
        or not isinstance(payload.get("message"), dict)
    ):
        return outcome("unresolved", "invalid_crossref_response")

    work = payload["message"]
    doi = normalize_doi(work.get("DOI"))
    if doi != requested:
        return outcome("unresolved", "identifier_not_confirmed")

    work_type = work.get("type")
    if not isinstance(work_type, str) or not work_type:
        return outcome("unresolved", "missing_work_type")
    if work_type not in _ARTICLE_TYPES:
        return outcome("ineligible", "non_article_work_type")

    title = _first_text(work.get("title"))
    venue = _first_text(work.get("container-title"))
    if title is None or venue is None:
        return outcome("unresolved", "missing_title_or_venue")

    dates = []
    try:
        for field in _PUBLICATION_DATE_FIELDS:
            if field in work and work[field] is not None:
                dates.append(_date_interval(work[field], field))
    except (ValueError, TypeError, OverflowError):
        return outcome("unresolved", "invalid_publication_date")

    if not dates:
        return outcome("unresolved", "publication_date_missing")

    confirmed = [
        item for item in dates
        if date.fromisoformat(item["latest"]) <= cutoff
    ]
    if not confirmed:
        if all(date.fromisoformat(item["earliest"]) > cutoff for item in dates):
            return outcome("ineligible", "publication_date_in_future")
        return outcome("unresolved", "publication_date_precision_insufficient")

    authors = work.get("author", [])
    if not isinstance(authors, list) or any(
        not isinstance(author, dict) for author in authors
    ):
        return outcome("unresolved", "invalid_author_metadata")

    selected = confirmed[0]
    paper = {
        "doi": doi,
        "pmid": None,
        "title": title,
        "authors": copy.deepcopy(authors),
        "venue": venue,
        "year": selected["date_parts"][0],
        "publication_date": selected,
        "publication_dates": dates,
        "work_type": work_type,
        "url": "https://doi.org/" + doi,
    }
    return outcome("published", "published_article_metadata_confirmed", paper)
