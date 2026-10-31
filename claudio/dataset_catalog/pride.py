"""Retrieve and save PRIDE project metadata and file inventories."""

import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Mapping, Union

import requests

from . import _write_json


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


FILE_COLUMNS = [
    "project_accession",
    "file_accession",
    "file_name",
    "file_category",
    "file_size_bytes",
    "checksum",
    "project_accessions",
    "public_file_locations",
]


def _validated_files(accession, file_pages):
    if not file_pages:
        raise ValueError("At least one file-list response is required")

    files = []
    expected_total = None
    seen_files = set()

    for page in file_pages:
        if not isinstance(page.data, list):
            raise ValueError("PRIDE file response must be a JSON array")
        if json.loads(page.raw_body) != page.data:
            raise ValueError("Parsed file metadata differs from its raw response")

        try:
            total = int(page.headers["total_records"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("Invalid PRIDE total_records header") from error
        if total < 0:
            raise ValueError("PRIDE file count must not be negative")
        if expected_total is None:
            expected_total = total
        elif total != expected_total:
            raise ValueError("PRIDE file counts differ between responses")

        for record in page.data:
            if not isinstance(record, dict):
                raise ValueError("PRIDE file record must be a JSON object")
            file_id = record.get("accession")
            if not isinstance(file_id, str) or not file_id.strip():
                raise ValueError("PRIDE file record has no valid accession")
            projects = record.get("projectAccessions")
            if not isinstance(projects, list) or accession not in projects:
                raise ValueError(f"PRIDE file does not belong to {accession}")
            if file_id in seen_files:
                raise ValueError(f"Repeated PRIDE file: {file_id}")
            seen_files.add(file_id)
            files.append(record)

    if len(files) != expected_total:
        raise ValueError(
            f"PRIDE file count mismatch: received {len(files)}, "
            f"expected {expected_total}"
        )
    return files, expected_total


def _save_response(output_dir, relative_path, response):
    (output_dir / relative_path).write_bytes(response.raw_body)
    return {
        "path": relative_path.as_posix(),
        "request_url": response.request_url,
        "retrieved_at_utc": response.retrieved_at_utc,
        "headers": dict(response.headers),
        "sha256": hashlib.sha256(response.raw_body).hexdigest(),
    }


def _file_row(accession, record):
    category = record.get("fileCategory")
    return {
        "project_accession": accession,
        "file_accession": record["accession"],
        "file_name": record.get("fileName", ""),
        "file_category": (
            category.get("value", "") if isinstance(category, dict) else ""
        ),
        "file_size_bytes": record.get("fileSizeBytes", ""),
        "checksum": record.get("checksum", ""),
        "project_accessions": json.dumps(
            record["projectAccessions"], ensure_ascii=False
        ),
        "public_file_locations": json.dumps(
            record.get("publicFileLocations") or [], ensure_ascii=False
        ),
    }


def write_pride_outputs(accession, project, file_pages, output_dir):
    """Write a complete PRIDE project snapshot into a new directory.

    Parameters
    ----------
    accession : str
        ProteomeXchange project accession.
    project : PrideResponse
        Retrieved project metadata.
    file_pages : list of PrideResponse
        Complete retrieved file-list responses.
    output_dir : str or pathlib.Path
        New directory for the snapshot.

    Returns
    -------
    output_dir : pathlib.Path
    """
    if (
        not isinstance(project.data, dict)
        or project.data.get("accession") != accession
    ):
        raise ValueError("PRIDE project accession mismatch")
    if json.loads(project.raw_body) != project.data:
        raise ValueError("Parsed project metadata differs from its raw response")

    files, reported_total = _validated_files(accession, file_pages)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "raw").mkdir()

    project_source = _save_response(
        output_dir, Path("raw") / "project.json", project
    )
    file_sources = [
        _save_response(
            output_dir,
            Path("raw") / f"files_page_{number:04d}.json",
            page,
        )
        for number, page in enumerate(file_pages)
    ]

    _write_json(output_dir / "project.json", project.data)
    _write_json(output_dir / "files.json", files)

    with (output_dir / "files.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=FILE_COLUMNS)
        writer.writeheader()
        for record in files:
            writer.writerow(_file_row(accession, record))

    manifest = {
        "schema_version": 1,
        "source": "PRIDE",
        "project_accession": accession,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "screening_status": "pending",
        "file_count": len(files),
        "reported_file_count": reported_total,
        "project_response": project_source,
        "file_responses": file_sources,
    }
    _write_json(output_dir / "manifest.json", manifest)
    return output_dir
