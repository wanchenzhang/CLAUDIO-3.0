"""Inspect candidates, save progress, and resume interrupted inventories."""

import csv
import hashlib
import json
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import requests

from . import _write_json
from .pride import (
    FILE_COLUMNS, PrideResponse, _file_row, _validated_files,
    fetch_pride_project, fetch_pride_file_pages, write_pride_outputs,
)


def inspect_candidate(candidate, output_dir):
    """Inspect one candidate and return its technical collection status.

    Parameters
    ----------
    candidate : dict
        Candidate record from the discovery catalog.
    output_dir : str or pathlib.Path
        New directory for this project's inventory.

    Returns
    -------
    result : dict
        Candidate metadata with inspection status and error information.
    """
    result = deepcopy(candidate)
    result.update({
        "inspection_status": "adapter_pending",
        "file_count": None,
        "inventory_path": "",
        "error_stage": "",
        "error_type": "",
        "error_message": "",
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "finished_at_utc": "",
    })

    if candidate.get("repository") != "PRIDE":
        result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        return result

    stage = "candidate"
    try:
        accession = candidate.get("accession")
        if not isinstance(accession, str) or not accession.strip():
            raise ValueError("Candidate has no valid project accession")

        stage = "project_metadata"
        project = fetch_pride_project(accession)

        stage = "file_list"
        pages = fetch_pride_file_pages(accession)

        stage = "output"
        destination = write_pride_outputs(
            accession, project, pages, output_dir
        )
        result.update({
            "inspection_status": "collected",
            "file_count": sum(len(page.data) for page in pages),
            "inventory_path": str(destination),
        })
    except (requests.exceptions.RequestException, ValueError, OSError) as error:
        result.update({
            "inspection_status": "failed",
            "error_stage": stage,
            "error_type": type(error).__name__,
            "error_message": str(error),
        })

    result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    return result


def _read_catalog(catalog_path):
    body = Path(catalog_path).read_bytes()
    candidates = json.loads(body)
    if not isinstance(candidates, list):
        raise ValueError("Candidate catalog must be a JSON array")
    json.dumps(candidates, allow_nan=False)

    seen = set()
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise ValueError("Each candidate must be a JSON object")

        accession = candidate.get("accession")
        if not isinstance(accession, str) or not accession.strip():
            raise ValueError("Each candidate must have a valid accession")
        if accession in seen:
            raise ValueError(f"Repeated candidate accession: {accession}")
        seen.add(accession)

    return body, candidates


def inspect_catalog(catalog_path, output_dir, *, resume=False):
    """Inspect catalog candidates and save progress after each result.

    Parameters
    ----------
    catalog_path : str or pathlib.Path
        Discovery candidates.json file.
    output_dir : str or pathlib.Path
        Directory for the batch inventory.
    resume : bool, optional
        Request recovery from an existing checkpoint.

    Returns
    -------
    results : list of dict
        One inspection result per candidate, in input order.
    """
    body, candidates = _read_catalog(catalog_path)
    output_dir = Path(output_dir).resolve()
    if resume:
        checkpoint = _load_checkpoint(output_dir, body, candidates)
        # Completed reports need no rewrite. A missing manifest can be
        # regenerated if execution stopped while publishing the reports.
        if len(checkpoint["results"]) == len(candidates):
            if not (output_dir / "manifest.json").is_file():
                finish_inventory(output_dir, checkpoint)
            return checkpoint["results"]
        _archive_reports(output_dir)
    else:
        output_dir.mkdir(parents=True, exist_ok=False)
        (output_dir / "projects").mkdir()
        (output_dir / "input_candidates.json").write_bytes(body)

        started_at = datetime.now(timezone.utc).isoformat()
        checkpoint = {
            "schema_version": 1,
            "status": "running",
            "source_catalog": {
                "path": "input_candidates.json",
                "original_path": str(Path(catalog_path).resolve()),
                "sha256": hashlib.sha256(body).hexdigest(),
            },
            "candidate_count": len(candidates),
            "started_at_utc": started_at,
            "updated_at_utc": started_at,
            "results": [],
        }
    checkpoint["status"] = "running"
    checkpoint["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_checkpoint(output_dir, checkpoint)

    completed = len(checkpoint["results"])
    for number, candidate in enumerate(candidates[completed:], start=completed + 1):
        project_path = Path("projects") / f"candidate_{number:06d}"
        attempt = 1
        relative_path = project_path / f"attempt_{attempt:04d}"
        while _saved_path(output_dir, relative_path.as_posix()).exists():
            attempt += 1
            relative_path = project_path / f"attempt_{attempt:04d}"
        result = inspect_candidate(
            candidate, _saved_path(output_dir, relative_path.as_posix())
        )
        result["candidate_index"] = number
        result["attempt_number"] = (
            attempt if candidate.get("repository") == "PRIDE" else 0
        )
        result["attempt_path"] = (
            relative_path.as_posix() if result["attempt_number"] else ""
        )

        if result["inventory_path"]:
            result["inventory_path"] = (
                Path(result["inventory_path"])
                .relative_to(output_dir)
                .as_posix()
            )

        checkpoint["results"].append(result)
        checkpoint["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        write_checkpoint(output_dir, checkpoint)

    finish_inventory(output_dir, checkpoint)
    return checkpoint["results"]


def _saved_path(root, relative_path):
    """Resolve a saved path within its inventory directory."""
    if not isinstance(relative_path, str) or not relative_path:
        raise ValueError("Missing saved inventory path")
    path = Path(relative_path)
    if path.is_absolute():
        raise ValueError("Saved inventory paths must be relative")
    resolved = (root / path).resolve()
    resolved.relative_to(root.resolve())
    return resolved


def _check_inventory(output_dir, result):
    """Check that a collected project's saved sources are still usable."""
    root = _saved_path(output_dir, result["inventory_path"])
    manifest = json.loads((root / "manifest.json").read_bytes())
    accession = result["accession"]
    if (
        manifest["schema_version"] != 1
        or manifest["status"] != "complete"
        or manifest["project_accession"] != accession
    ):
        raise ValueError(f"Incomplete or mismatched inventory for {accession}")

    sources = [manifest["project_response"], *manifest["file_responses"]]
    responses = []
    for source in sources:
        raw = _saved_path(root, source["path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != source["sha256"]:
            raise ValueError(f"Saved response checksum mismatch for {accession}")
        responses.append(PrideResponse(
            request_url=source["request_url"],
            retrieved_at_utc=source["retrieved_at_utc"],
            raw_body=raw,
            headers=requests.structures.CaseInsensitiveDict(source["headers"]),
            data=json.loads(raw),
        ))

    project = responses[0].data
    if not isinstance(project, dict) or project.get("accession") != accession:
        raise ValueError(f"Saved project identity mismatch for {accession}")
    files, total = _validated_files(accession, responses[1:])
    if (
        json.loads((root / "project.json").read_bytes()) != project
        or json.loads((root / "files.json").read_bytes()) != files
        or manifest["file_count"] != total
        or manifest["reported_file_count"] != total
        or result["file_count"] != total
    ):
        raise ValueError(f"Saved project files differ from their sources: {accession}")

    expected_rows = [
        {key: "" if value is None else str(value)
         for key, value in _file_row(accession, record).items()}
        for record in files
    ]
    with (root / "files.csv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != FILE_COLUMNS or list(reader) != expected_rows:
            raise ValueError(f"Saved file summary differs from its sources: {accession}")


def _load_checkpoint(output_dir, catalog_body, candidates):
    """Validate saved progress before writing files or requesting data."""
    try:
        checkpoint = json.loads((output_dir / "checkpoint.json").read_bytes())
        if checkpoint["schema_version"] != 1:
            raise ValueError("Unsupported inventory checkpoint version")
        if checkpoint["status"] not in ("running", "finished"):
            raise ValueError("Invalid inventory checkpoint status")
        source = checkpoint["source_catalog"]
        saved_catalog = _saved_path(output_dir, source["path"]).read_bytes()
        if (
            hashlib.sha256(catalog_body).hexdigest() != source["sha256"]
            or hashlib.sha256(saved_catalog).hexdigest() != source["sha256"]
            or checkpoint["candidate_count"] != len(candidates)
        ):
            raise ValueError("Input catalog differs from the saved catalog")

        results = checkpoint["results"]
        if not isinstance(results, list) or len(results) > len(candidates):
            raise ValueError("Invalid checkpoint results")
        if checkpoint["status"] == "finished" and len(results) != len(candidates):
            raise ValueError("Finished checkpoint has missing candidate results")
        for number, (candidate, result) in enumerate(zip(candidates, results), start=1):
            if (
                not isinstance(result, dict)
                or result["candidate_index"] != number
                or any(key not in result or result[key] != value
                       for key, value in candidate.items())
            ):
                raise ValueError("Checkpoint results do not match the catalog")
            if result["inspection_status"] == "collected":
                _check_inventory(output_dir, result)
            elif result["inspection_status"] in ("failed", "adapter_pending"):
                if result["file_count"] is not None or result["inventory_path"]:
                    raise ValueError("Uncollected candidate has a saved inventory")
            else:
                raise ValueError("Invalid checkpoint inspection status")
        return checkpoint
    except (OSError, KeyError, TypeError, ValueError, AttributeError) as error:
        raise ValueError(f"Cannot resume inventory: {error}") from error


def _archive_reports(output_dir):
    """Keep the previous batch reports before continuing an interrupted run."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    history = _saved_path(output_dir, f"history/resume_{stamp}")
    history.mkdir(parents=True, exist_ok=False)
    for name in ("checkpoint.json", "manifest.json", "projects.json", "projects.csv"):
        path = output_dir / name
        if path.is_file():
            (history / name).write_bytes(path.read_bytes())
    (output_dir / "manifest.json").unlink(missing_ok=True)


PROJECT_COLUMNS = [
    "accession",
    "repository",
    "title",
    "search_terms",
    "review_status",
    "inspection_status",
    "file_count",
    "inventory_path",
    "attempt_path",
    "error_stage",
    "error_type",
    "error_message",
    "started_at_utc",
    "finished_at_utc",
]


def write_checkpoint(output_dir, checkpoint):
    """Replace the checkpoint after the new JSON has been fully written."""
    output_dir = Path(output_dir)
    temporary_path = output_dir / ".checkpoint.json.tmp"
    _write_json(temporary_path, checkpoint)
    temporary_path.replace(output_dir / "checkpoint.json")


def finish_inventory(output_dir, checkpoint):
    """Write final reports and publish the finished batch manifest."""
    output_dir = Path(output_dir)
    results = checkpoint["results"]
    if len(results) != checkpoint["candidate_count"]:
        raise ValueError(
            "Cannot finish an inventory with missing candidate results"
        )

    _write_json(output_dir / "projects.json", results)

    with (output_dir / "projects.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=PROJECT_COLUMNS)
        writer.writeheader()
        for result in results:
            row = {
                column: result.get(column, "")
                for column in PROJECT_COLUMNS
            }
            row["search_terms"] = json.dumps(
                result.get("search_terms", []), ensure_ascii=False
            )
            writer.writerow(row)

    checkpoint["status"] = "finished"
    checkpoint["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
    write_checkpoint(output_dir, checkpoint)

    counts = Counter(row["inspection_status"] for row in results)
    manifest = {
        "schema_version": 1,
        "status": "finished",
        "source_catalog": checkpoint["source_catalog"],
        "candidate_count": checkpoint["candidate_count"],
        "result_count": len(results),
        "status_counts": {
            status: counts[status]
            for status in ("collected", "failed", "adapter_pending")
        },
        "started_at_utc": checkpoint["started_at_utc"],
        "finished_at_utc": checkpoint["updated_at_utc"],
        "checkpoint": "checkpoint.json",
        "results": "projects.json",
        "summary_csv": "projects.csv",
    }
    _write_json(output_dir / "manifest.json", manifest)
