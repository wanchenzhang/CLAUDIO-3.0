"""Search ProteomeCentral and save candidate catalogs with their sources."""

import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List

import requests

from . import _write_json


DATASETS_URL = (
    "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets"
)


@dataclass(frozen=True)
class CatalogPage:
    request_url: str
    retrieved_at_utc: str
    raw_body: bytes
    data: dict


def fetch_dataset_pages(
    search_term: str, *, page_size: int = 100
) -> List[CatalogPage]:
    """Fetch all matching pages, retaining their original response bodies."""
    pages = []
    page_number = 1
    expected_total = None
    expected_pages = None

    while True:
        params = {
            "search": search_term,
            "resultType": "resultset",
            "pageNumber": page_number,
            "pageSize": page_size,
        }
        response = requests.get(
            DATASETS_URL, params=params, timeout=(15, 120)
        )
        response.raise_for_status()
        retrieved_at = datetime.now(timezone.utc).isoformat()
        data = json.loads(response.content)

        if data["status"]["status"] != "OK":
            raise ValueError(
                f"ProteomeCentral search failed: {data['status']}"
            )

        result_set = data["result_set"]
        if result_set["page_number"] != page_number:
            raise ValueError(
                "ProteomeCentral returned an unexpected page"
            )

        total = result_set["n_available_rows"]
        page_count = result_set["n_available_pages"]
        if expected_total is None:
            expected_total = total
            expected_pages = page_count
        elif total != expected_total or page_count != expected_pages:
            raise ValueError(
                "ProteomeCentral result counts changed during search"
            )

        records = data["datasets"]
        if len(records) != result_set["n_rows_returned"]:
            raise ValueError(
                "ProteomeCentral returned an inconsistent row count"
            )
        if not records and expected_total > 0:
            raise ValueError(
                "ProteomeCentral returned an unexpected empty page"
            )

        pages.append(
            CatalogPage(
                request_url=response.url,
                retrieved_at_utc=retrieved_at,
                raw_body=response.content,
                data=data,
            )
        )

        if page_number >= expected_pages:
            break
        page_number += 1

    received = sum(len(page.data["datasets"]) for page in pages)
    if received != expected_total:
        raise ValueError(
            f"Incomplete download: received {received} "
            f"of {expected_total} records"
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
                    raise ValueError(
                        "Dataset row does not match its column headers"
                    )

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
                if search_term not in candidate["search_terms"]:
                    candidate["search_terms"].append(search_term)

                variants = candidate["metadata_variants"]
                for variant in variants:
                    if variant["metadata"] == metadata:
                        if search_term not in variant["search_terms"]:
                            variant["search_terms"].append(search_term)
                        break
                else:
                    variants.append({
                        "metadata": metadata,
                        "search_terms": [search_term],
                    })

                candidate["metadata_conflict"] = len(variants) > 1

    return [candidates[accession] for accession in sorted(candidates)]


CSV_COLUMNS = [
    "accession",
    "title",
    "repository",
    "species",
    "search_terms",
    "review_status",
    "metadata_conflict",
    "metadata_variant_count",
]




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

    for search_number, (search_term, pages) in enumerate(
        searches.items(), start=1
    ):
        search = {
            "search_term": search_term,
            "reported_total": pages[0].data["result_set"]["n_available_rows"],
            "downloaded_records": sum(
                len(page.data["datasets"]) for page in pages
            ),
            "pages": [],
        }

        for page_number, page in enumerate(pages, start=1):
            relative_path = Path("raw") / (
                f"search_{search_number:03d}_page_{page_number:04d}.json"
            )
            (output_dir / relative_path).write_bytes(page.raw_body)
            search["pages"].append({
                "path": relative_path.as_posix(),
                "request_url": page.request_url,
                "retrieved_at_utc": page.retrieved_at_utc,
                "query": page.data["query"],
                "sha256": hashlib.sha256(page.raw_body).hexdigest(),
            })

        manifest["searches"].append(search)

    _write_json(output_dir / "candidates.json", candidates)

    with (output_dir / "candidates.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for candidate in candidates:
            row = {
                column: candidate.get(column, "")
                for column in CSV_COLUMNS
            }
            row["search_terms"] = json.dumps(
                candidate["search_terms"], ensure_ascii=False
            )
            row["metadata_variant_count"] = len(
                candidate["metadata_variants"]
            )
            writer.writerow(row)

    _write_json(output_dir / "manifest.json", manifest)
    return output_dir
