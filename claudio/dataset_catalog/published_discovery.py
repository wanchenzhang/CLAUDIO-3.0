
"""Published-only incremental discovery with explicit paper and source evidence."""

from collections import Counter
from copy import deepcopy
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tempfile
from urllib.parse import parse_qs, urlsplit

from . import _write_json
from .discovery_history import DiscoveryHistory, default_history_directory
from .proteomecentral import CSV_COLUMNS, DATASETS_URL, fetch_dataset_pages, merge_search_results
from .publications import extract_publication_links
from .publication_lookup import PublicationLookup
from .publication_verification import _as_of_date


PUBLISHED_CSV_COLUMNS = CSV_COLUMNS + [
    "publication_status", "paper_dois_json", "paper_pmids_json", "papers_json",
]


def _canonical_accession(value):
    from .discovery_history import _accession
    return _accession(value)


def _extend_unique(destination, values):
    for value in values:
        if value not in destination:
            destination.append(deepcopy(value))


def _retain_pending_context(candidate, indicators):
    """Retain an exact bare Pending placeholder as context, never as a paper."""
    retained = deepcopy(indicators)

    def visit(value, path):
        if isinstance(value, str) and value.strip().casefold() == "pending":
            existing = next((item for item in retained if item.get("kind") == "pending"
                             and item.get("raw_value") == value), None)
            if existing is None:
                existing = {"kind": "pending", "raw_value": value, "source_fields": []}
                retained.append(existing)
            _extend_unique(existing["source_fields"], [path])
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, path + "[" + str(index) + "]")

    sources = [("candidate", candidate)]
    for index, variant in enumerate(candidate.get("metadata_variants") or []):
        if isinstance(variant, dict) and isinstance(variant.get("metadata"), dict):
            sources.append(("candidate.metadata_variants[" + str(index) + "].metadata", variant["metadata"]))
    for path, metadata in sources:
        for field in ("publications", "references"):
            visit(metadata.get(field), path + "." + field)
    return retained


def _add_paper(papers, paper):
    for previous in papers:
        same = any(previous.get(key) and previous.get(key) == paper.get(key) for key in ("doi", "pmid"))
        conflict = any(previous.get(key) and paper.get(key) and previous[key] != paper[key] for key in ("doi", "pmid"))
        if same and not conflict:
            evidence = deepcopy(previous["verification_evidence"])
            fields = deepcopy(previous["association_source_fields"])
            if paper.get("pmid") and not previous.get("pmid"):
                previous.update(deepcopy(paper))
            _extend_unique(evidence, paper["verification_evidence"])
            _extend_unique(fields, paper["association_source_fields"])
            previous["verification_evidence"] = evidence
            previous["association_source_fields"] = fields
            return
    papers.append(deepcopy(paper))


def qualify_dataset_candidate(candidate, lookup, *, as_of):
    """Qualify explicit source references; never infer papers from a title or dataset DOI."""
    as_of = _as_of_date(as_of).isoformat()
    accession = _canonical_accession(candidate.get("accession"))
    result = deepcopy(candidate)
    result.update(accession=accession, publication_status="unresolved", papers=[])
    review = {"as_of": as_of, "association_status": "not_established",
              "evidence_root": "publication_evidence", "indicators": [], "references": [], "checks": []}
    result["publication_review"] = review
    try:
        extracted = extract_publication_links(candidate)
    except ValueError as error:
        review.update(reason="invalid_source_publication_metadata", error=str(error))
        return result
    review["indicators"] = _retain_pending_context(candidate, extracted["indicators"])
    review["references"] = deepcopy(extracted["references"])
    if not extracted["references"]:
        kinds = {item["kind"] for item in review["indicators"]}
        result["publication_status"] = ("pending" if "pending" in kinds else
                                         "unpublished" if "no_publication" in kinds else "unresolved")
        review["reason"] = "no_explicit_article_identifier"
        return result
    statuses = []
    for reference in extracted["references"]:
        checked = lookup.lookup(doi=reference.get("doi"), pmid=reference.get("pmid"), as_of=as_of)
        statuses.append(checked["status"])
        review["checks"].append({"reference": deepcopy(reference), "qualification": deepcopy(checked)})
        if checked["status"] != "published":
            continue
        if not isinstance(checked.get("paper"), dict) or not isinstance(checked.get("evidence"), dict):
            raise ValueError("Published article qualification is missing captured HTTP evidence")
        paper = deepcopy(checked["paper"])
        paper["association_source_fields"] = deepcopy(reference["source_fields"])
        paper["verification_evidence"] = [{
            "root": "publication_evidence", **deepcopy(checked["evidence"]),
            "provider": checked["provider"], "requested_doi": reference.get("doi"),
            "requested_pmid": reference.get("pmid"),
        }]
        _add_paper(result["papers"], paper)
    if result["papers"]:
        result["publication_status"] = "published"
        review.update(association_status="explicit_source_reference", reason="published_article_confirmed")
    elif "unresolved" in statuses:
        result["publication_status"] = "unresolved"
        review["reason"] = "article_qualification_unresolved"
    else:
        result["publication_status"] = "ineligible"
        review["reason"] = "no_qualifying_published_article"
    return result


def _save_search(staging, term, pages, search_number, page_size):
    if not isinstance(pages, list) or not pages:
        raise ValueError("Each search must have complete response pages")
    expected_counts = None
    downloaded = 0
    sources = []
    for number, page in enumerate(pages, 1):
        try:
            raw_data = json.loads(page.raw_body)
            if raw_data != page.data:
                raise ValueError("Parsed registry response differs from its raw body")
            uri = urlsplit(page.request_url)
            if (uri.scheme != "https" or uri.hostname != "proteomecentral.proteomexchange.org"
                    or uri.port not in (None, 443) or uri.username or uri.password or uri.fragment
                    or uri.path != "/api/proxi/v0.1/datasets"):
                raise ValueError("Registry source URL is not the official dataset search")
            query = parse_qs(uri.query)
            expected_query = {"search": term, "resultType": "resultset",
                              "pageNumber": str(number), "pageSize": str(page_size)}
            if any(query.get(key) != [value] for key, value in expected_query.items()):
                raise ValueError("Registry source URL does not match the requested search")
            if raw_data["status"]["status"] != "OK":
                raise ValueError("Registry search did not report success")
            result_set = raw_data["result_set"]
            counts = result_set["n_available_rows"], result_set["n_available_pages"]
            if any(type(value) is not int or value < 0 for value in counts):
                raise ValueError("Invalid registry result counts")
            if expected_counts is None:
                expected_counts = counts
            elif counts != expected_counts:
                raise ValueError("Registry result counts changed between pages")
            rows = raw_data["datasets"]
            columns = result_set["datasets_title_list"]
            if (not isinstance(rows, list) or not isinstance(columns, list)
                    or any(not isinstance(column, str) for column in columns)
                    or len(columns) != len(set(columns)) or "dataset identifier" not in columns
                    or any(not isinstance(row, list) or len(row) != len(columns) for row in rows)):
                raise ValueError("Malformed registry dataset rows")
            if (type(result_set["page_number"]) is not int or result_set["page_number"] != number
                    or type(result_set["n_rows_returned"]) is not int or result_set["n_rows_returned"] != len(rows)):
                raise ValueError("Registry page identity or row count mismatch")
            downloaded += len(rows)
        except (KeyError, TypeError, UnicodeError, AttributeError) as error:
            raise ValueError("Invalid registry response structure") from error
        name = f"raw/search_{search_number:03d}_page_{number:04d}.json"
        (staging / name).write_bytes(page.raw_body)
        sources.append({"path": name, "request_url": page.request_url,
                        "retrieved_at_utc": page.retrieved_at_utc, "query": deepcopy(raw_data["query"]),
                        "sha256": hashlib.sha256(page.raw_body).hexdigest()})
    if len(pages) != max(1, expected_counts[1]) or downloaded != expected_counts[0]:
        raise ValueError("Incomplete registry search page or record counts")
    return {"search_term": term, "reported_total": expected_counts[0],
            "downloaded_records": downloaded, "pages": sources}


def _write_candidate_csv(path, candidates):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PUBLISHED_CSV_COLUMNS)
        writer.writeheader()
        for candidate in candidates:
            row = {column: candidate.get(column, "") for column in CSV_COLUMNS}
            row["search_terms"] = json.dumps(candidate["search_terms"], ensure_ascii=False)
            row["metadata_variant_count"] = len(candidate["metadata_variants"])
            row["publication_status"] = candidate["publication_status"]
            for name, key in (("paper_dois_json", "doi"), ("paper_pmids_json", "pmid")):
                row[name] = json.dumps(list(dict.fromkeys(
                    paper[key] for paper in candidate["papers"] if paper.get(key))), ensure_ascii=False)
            row["papers_json"] = json.dumps(candidate["papers"], ensure_ascii=False, allow_nan=False)
            writer.writerow(row)


def run_published_discovery(search_terms, output_dir, *, history_dir=None, page_size=100,
                            publication_request_limit=1000, as_of=None, progress=None):
    """Publish unseen datasets with confirmed papers and retain a complete review audit.

    The registry's explicit publication references establish association. Missing
    references remain unresolved for later review; no title search or release-date
    inference is used. This function does not retrieve repository file inventories.
    """
    if not isinstance(search_terms, (list, tuple)) or not search_terms:
        raise ValueError("At least one search term is required")
    if any(not isinstance(term, str) or not term.strip() for term in search_terms):
        raise ValueError("Search terms must be nonempty strings")
    terms = list(dict.fromkeys(term.strip() for term in search_terms))
    if type(page_size) is not int or page_size < 1:
        raise ValueError("page_size must be a positive integer")
    if type(publication_request_limit) is not int or not 1 <= publication_request_limit <= 10000:
        raise ValueError("publication_request_limit must be between 1 and 10000")
    as_of = _as_of_date(as_of if as_of is not None else datetime.now(timezone.utc).date()).isoformat()
    if _as_of_date(as_of) > datetime.now(timezone.utc).date():
        raise ValueError("Publication qualification date must not be in the future")
    output = Path(output_dir).expanduser().resolve()
    if output.exists():
        raise FileExistsError("Choose a new discovery output directory")
    history_path = history_dir if history_dir is not None else default_history_directory()
    with DiscoveryHistory(history_path) as history:
        seen = history.seen_accessions()
        run = history.prepare_run(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=".claudio_discovery_" + run["run_id"] + "_", dir=output.parent))
        (staging / "raw").mkdir()
        _write_json(staging / "discovery_context.json", {
            "run_id": run["run_id"], "namespace_id": run["namespace_id"],
            "as_of": as_of, "publication_request_limit": publication_request_limit,
            "association_scope": "explicit_registry_publication_references",
        })
        searches = {}
        sources = []
        for number, term in enumerate(terms, 1):
            pages = fetch_dataset_pages(term, page_size=page_size)
            sources.append(_save_search(staging, term, pages, number, page_size))
            searches[term] = pages
            if progress:
                progress(f"{term}: {sources[-1]['downloaded_records']} records")
        originals = merge_search_results(searches)
        for candidate in originals:
            candidate["accession"] = _canonical_accession(candidate["accession"])
        if len({candidate["accession"] for candidate in originals}) != len(originals):
            raise ValueError("Registry candidates contain duplicate canonical accessions")
        _write_json(staging / "source_candidates.json", originals)
        reviewed = []
        with PublicationLookup(staging / "publication_evidence", max_requests=publication_request_limit) as lookup:
            for index, candidate in enumerate(originals, 1):
                if progress and (index == 1 or index % 10 == 0 or index == len(originals)):
                    progress(f"Publication review: {index}/{len(originals)}; HTTP attempts: {lookup.request_count}")
                if candidate["accession"] in seen:
                    annotated = deepcopy(candidate)
                    annotated.update(publication_status="unverified", papers=[],
                                     discovery_decision="excluded_seen",
                                     publication_review={"reason": "already_emitted_in_history",
                                                         "current_publications_rechecked": False})
                else:
                    annotated = qualify_dataset_candidate(candidate, lookup, as_of=as_of)
                    annotated["discovery_decision"] = ("eligible" if annotated["publication_status"] == "published"
                                                        else "excluded_publication_status")
                reviewed.append(annotated)
            transport = {"request_limit": publication_request_limit, "http_attempts": lookup.request_count,
                         "captured_body_bytes": lookup.captured_body_bytes}

        def write_outputs(destination, selected):
            # Renames retain the captured bytes and leave the staging directory
            # empty. Failed output directories remain available for inspection.
            for path in sorted(staging.iterdir()):
                path.rename(destination / path.name)
            emitted = {candidate["accession"] for candidate in selected}
            audit = deepcopy(reviewed)
            for candidate in audit:
                if candidate["publication_status"] == "published":
                    candidate["discovery_decision"] = ("emitted" if candidate["accession"] in emitted else "excluded_seen")
            _write_json(destination / "candidate_audit.json", audit)
            _write_json(destination / "candidates.json", selected)
            _write_candidate_csv(destination / "candidates.csv", selected)
            return {
                "source": "ProteomeCentral", "searches": sources, "source_candidate_count": len(originals),
                "publication_as_of": as_of, "publication_transport": transport,
                "association_scope": "explicit_registry_publication_references",
                "publication_status_counts": dict(Counter(candidate["publication_status"] for candidate in audit)),
                "discovery_decision_counts": dict(Counter(candidate["discovery_decision"] for candidate in audit)),
            }

        history.publish(run["run_id"], reviewed, write_outputs)
        staging.rmdir()
        return output
