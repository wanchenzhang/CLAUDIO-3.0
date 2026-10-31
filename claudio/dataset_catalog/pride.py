"""Retrieve project metadata and file listings from PRIDE."""

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import List, Mapping, Union

import requests


PRIDE_API_URL = "https://www.ebi.ac.uk/pride/ws/archive/v3"


@dataclass(frozen=True)
class PrideResponse:
    request_url: str
    retrieved_at_utc: str
    raw_body: bytes
    headers: Mapping[str, str]
    data: Union[dict, list]


def fetch_pride_project(accession: str) -> PrideResponse:
    """Fetch a PRIDE project and retain its original response.

    Parameters
    ----------
    accession : str
        ProteomeXchange project accession.

    Returns
    -------
    project : PrideResponse
    """
    response = requests.get(
        f"{PRIDE_API_URL}/projects/{accession}", timeout=(15, 120)
    )
    response.raise_for_status()
    retrieved_at = datetime.now(timezone.utc).isoformat()
    data = json.loads(response.content)

    if not isinstance(data, dict):
        raise ValueError("PRIDE project response must be a JSON object")
    if data.get("accession") != accession:
        raise ValueError(
            f"PRIDE project accession mismatch: requested {accession}, "
            f"received {data.get('accession')}"
        )

    return PrideResponse(
        request_url=response.url,
        retrieved_at_utc=retrieved_at,
        raw_body=response.content,
        headers=response.headers.copy(),
        data=data,
    )
def fetch_pride_file_pages(
    accession: str, *, page_size: int = 100
) -> List[PrideResponse]:
    """Fetch every indexed file page for a PRIDE project.

    Parameters
    ----------
    accession : str
        ProteomeXchange project accession.
    page_size : int
        Number of file records requested per page.

    Returns
    -------
    pages : list of PrideResponse
    """
    if (
        not isinstance(page_size, int)
        or isinstance(page_size, bool)
        or page_size < 1
    ):
        raise ValueError("page_size must be a positive integer")

    pages = []
    page_number = 0
    expected_total = None
    seen_files = set()

    while True:
        response = requests.get(
            f"{PRIDE_API_URL}/projects/{accession}/files",
            params={"pageSize": page_size, "page": page_number},
            timeout=(15, 120),
        )
        response.raise_for_status()
        retrieved_at = datetime.now(timezone.utc).isoformat()
        records = json.loads(response.content)

        if not isinstance(records, list):
            raise ValueError("PRIDE file response must be a JSON array")

        try:
            total = int(response.headers["total_records"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(
                "PRIDE file response has no valid total_records header"
            ) from error
        if total < 0:
            raise ValueError("PRIDE file count must not be negative")

        if expected_total is None:
            expected_total = total
        elif total != expected_total:
            raise ValueError("PRIDE file count changed during retrieval")

        for record in records:
            if not isinstance(record, dict):
                raise ValueError("PRIDE file record must be a JSON object")

            file_id = record.get("accession")
            if not isinstance(file_id, str) or not file_id.strip():
                raise ValueError("PRIDE file record has no valid accession")

            projects = record.get("projectAccessions")
            if not isinstance(projects, list) or accession not in projects:
                raise ValueError(
                    f"PRIDE file does not belong to {accession}: {file_id}"
                )
            if file_id in seen_files:
                raise ValueError(f"PRIDE returned a repeated file: {file_id}")
            seen_files.add(file_id)

        received = len(seen_files)
        if received > expected_total:
            raise ValueError("PRIDE returned more files than its reported total")
        if not records and received < expected_total:
            raise ValueError(
                f"Incomplete PRIDE file list: received {received} "
                f"of {expected_total} files"
            )

        pages.append(
            PrideResponse(
                request_url=response.url,
                retrieved_at_utc=retrieved_at,
                raw_body=response.content,
                headers=response.headers.copy(),
                data=records,
            )
        )

        if received == expected_total:
            return pages
        page_number += 1