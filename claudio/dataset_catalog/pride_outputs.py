"""Write project metadata and file inventories retrieved from PRIDE."""

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from .outputs import _write_json


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
            category = record.get("fileCategory")
            writer.writerow({
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
            })

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