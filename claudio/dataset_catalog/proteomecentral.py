"""Retrieve dataset search results from ProteomeCentral."""

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import List

import requests


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