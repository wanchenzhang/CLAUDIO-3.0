"""Search ProteomeCentral and save candidate catalogs with their sources."""

import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List
from urllib.parse import parse_qs, urlsplit
import requests

from . import _write_json


DATASETS_URL = "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets"
CSV_COLUMNS = [
    "accession", "title", "repository", "species", "search_terms",
    "review_status", "metadata_conflict", "metadata_variant_count",
]


@dataclass(frozen=True)
class CatalogPage:
    request_url: str
    retrieved_at_utc: str
    raw_body: bytes
    data: dict

def validate_search_page(page, term, *, page_number=None, page_size=None):
    """Validate the registry URL, page identity, and dataset rows."""
    if not isinstance(page.request_url, str):
        raise ValueError("Missing registry source URL")

    uri = urlsplit(page.request_url)
    if (
        uri.scheme != "https"
        or uri.hostname != "proteomecentral.proteomexchange.org"
        or uri.port not in (None, 443)
        or uri.username or uri.password or uri.fragment
        or uri.path != "/api/proxi/v0.1/datasets"
    ):
        raise ValueError("Invalid registry source URL")

    try:
        data = page.data
        result = data["result_set"]
        number = result["page_number"]
        size = result["page_size"]
        if (
            type(number) is not int or number < 1
            or type(size) is not int or size < 1
            or (page_number is not None and number != page_number)
            or (page_size is not None and size != page_size)
        ):
            raise ValueError("Registry page identity mismatch")

        query = parse_qs(uri.query)
        expected = {
            "search": term,
            "resultType": "resultset",
            "pageNumber": str(number),
            "pageSize": str(size),
        }
        if any(query.get(key) != [value] for key, value in expected.items()):
            raise ValueError("Registry source URL does not match the requested search")

        if data["status"]["status"] != "OK":
            raise ValueError("Registry search did not report success")

        columns = result["datasets_title_list"]
        rows = data["datasets"]
        if (
            not isinstance(columns, list)
            or any(not isinstance(column, str) for column in columns)
            or len(columns) != len(set(columns))
            or "dataset identifier" not in columns
            or not isinstance(rows, list)
            or any(not isinstance(row, list) or len(row) != len(columns) for row in rows)
        ):
            raise ValueError("Malformed registry dataset rows")

        returned = result["n_rows_returned"]
        if type(returned) is not int or returned != len(rows):
            raise ValueError("Registry row count mismatch")
    except (KeyError, TypeError, AttributeError) as error:
        raise ValueError("Invalid registry response structure") from error
def fetch_dataset_pages(search_term: str, *, page_size: int = 100) -> List[CatalogPage]:
    """Fetch all matching pages, retaining their original response bodies."""
    pages = []
    expected_counts = None

    while True:
        page_number = len(pages) + 1
        params = {
            "search": search_term,
            "resultType": "resultset",
            "pageNumber": page_number,
            "pageSize": page_size,
        }
        response = requests.get(DATASETS_URL, params=params, timeout=(15, 120))
        response.raise_for_status()
        retrieved_at = datetime.now(timezone.utc).isoformat()
        data = json.loads(response.content)

        if data["status"]["status"] != "OK":
            raise ValueError(f"ProteomeCentral search failed: {data['status']}")

        result_set = data["result_set"]
        if result_set["page_number"] != page_number:
            raise ValueError("ProteomeCentral returned an unexpected page")

        counts = (result_set["n_available_rows"], result_set["n_available_pages"])
        if expected_counts is None:
            expected_counts = counts
        elif counts != expected_counts:
            raise ValueError("ProteomeCentral result counts changed during search")

        records = data["datasets"]
        if len(records) != result_set["n_rows_returned"]:
            raise ValueError("ProteomeCentral returned an inconsistent row count")
        if not records and counts[0] > 0:
            raise ValueError("ProteomeCentral returned an unexpected empty page")

        pages.append(CatalogPage(
            request_url=response.url,
            retrieved_at_utc=retrieved_at,
            raw_body=response.content,
            data=data,
        ))
        if page_number >= counts[1]:
            break

    received = sum(len(page.data["datasets"]) for page in pages)
    if received != expected_counts[0]:
        raise ValueError(
            f"Incomplete download: received {received} of {expected_counts[0]} records"
        )
    return pages


def merge_search_results(searches):
    """Merge candidates and retain metadata variants with their search sources."""
    candidates = {}
    for search_term, pages in searches.items():
        for page in pages:
            columns = page.data["result_set"]["datasets_title_list"]
            for values in page.data["datasets"]:
                if len(values) != len(columns):
                    raise ValueError("Dataset row does not match its column headers")

                metadata = dict(zip(columns, values))
                accession = metadata.pop("dataset identifier")
                if accession not in candidates:
                    candidates[accession] = {
                        **metadata,
                        "accession": accession,
                        "search_terms": [],
                        "review_status": "pending",
                        "metadata_conflict": False,
                        "metadata_variants": [],
                    }

                candidate = candidates[accession]
                variants = candidate["metadata_variants"]
                variant = next(
                    (item for item in variants if item["metadata"] == metadata), None
                )
                if variant is None:
                    variant = {"metadata": metadata, "search_terms": []}
                    variants.append(variant)

                for record in (candidate, variant):
                    if search_term not in record["search_terms"]:
                        record["search_terms"].append(search_term)
                candidate["metadata_conflict"] = len(variants) > 1

    return [candidates[accession] for accession in sorted(candidates)]

def write_search_pages(output_dir, search_term, pages, search_number):
    search = {
        "search_term": search_term,
        "reported_total": pages[0].data["result_set"]["n_available_rows"],
        "downloaded_records": sum(len(page.data["datasets"]) for page in pages),
        "pages": [],
    }
    for page_number, page in enumerate(pages, start=1):
        filename = f"search_{search_number:03d}_page_{page_number:04d}.json"
        relative_path = Path("raw") / filename
        (output_dir / relative_path).write_bytes(page.raw_body)
        search["pages"].append({
            "path": relative_path.as_posix(),
            "request_url": page.request_url,
            "retrieved_at_utc": page.retrieved_at_utc,
            "query": page.data["query"],
            "sha256": hashlib.sha256(page.raw_body).hexdigest(),
        })
    return search
def write_discovery_outputs(searches, output_dir):
    """Write a discovery snapshot into a new directory.

    Parameters
    ----------
    searches : dict
        Search terms mapped to their retrieved CatalogPage objects.
    output_dir : str or pathlib.Path
        New directory for the snapshot.

    Returns
    -------
    output_dir : pathlib.Path
    """
    if not searches or any(not pages for pages in searches.values()):
        raise ValueError("Each search must include at least one response page")

    candidates = merge_search_results(searches)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "raw").mkdir()
    manifest = {
        "schema_version": 1,
        "source": "ProteomeCentral",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "candidate_count": len(candidates),
        "searches": [],
    }

    for search_number, (search_term, pages) in enumerate(searches.items(), start=1):
        manifest["searches"].append(
            write_search_pages(output_dir, search_term, pages, search_number)
        )

    _write_json(output_dir / "candidates.json", candidates)
    with (output_dir / "candidates.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for candidate in candidates:
            row = {column: candidate.get(column, "") for column in CSV_COLUMNS}
            row["search_terms"] = json.dumps(candidate["search_terms"], ensure_ascii=False)
            row["metadata_variant_count"] = len(candidate["metadata_variants"])
            writer.writerow(row)

    _write_json(output_dir / "manifest.json", manifest)
    return output_dir
