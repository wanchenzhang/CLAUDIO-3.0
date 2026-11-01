"""Inspect candidates, save progress, and resume interrupted inventories."""

import csv
import hashlib
import io
import json
import math
import re
import time
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

import requests

from . import _write_json
from .repositories import (
    FILE_COLUMNS, RepositoryResponse, _file_row, _repository_files, _json_response_data,
    FILE_SCREENING_RULE_VERSION, classify_inventory_file, download_pilot_file,
    _massive_identifiers, _parse_massive_file_page,
    _jpost_identifiers, _select_jpost_project, _iprox_identifiers, _select_proxi_project,
    _panorama_project, fetch_repository_inventory, write_repository_outputs,
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

    repository = candidate.get("repository")
    if repository not in ("PRIDE", "MassIVE", "jPOST", "iProX", "PanoramaPublic"):
        result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        return result

    stage = "candidate"
    try:
        accession = candidate.get("accession")
        if not isinstance(accession, str) or not accession.strip():
            raise ValueError("Candidate has no valid project accession")

        def set_stage(value):
            nonlocal stage
            stage = value

        project, pages = fetch_repository_inventory(candidate, on_stage=set_stage)

        stage = "output"
        destination = write_repository_outputs(
            repository, accession, project, pages, output_dir
        )
        manifest = json.loads((destination / "manifest.json").read_bytes())
        result.update({
            "inspection_status": "collected",
            "file_count": manifest["file_count"],
            "inventory_path": str(destination),
        })
    except Exception as error:
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


def inspect_catalog(
    catalog_path, output_dir, *, resume=False, delay=0, on_progress=None
):
    """Inspect catalog candidates and save progress after each result.

    Parameters
    ----------
    catalog_path : str or pathlib.Path
        Discovery candidates.json file.
    output_dir : str or pathlib.Path
        Directory for the batch inventory.
    resume : bool, optional
        Request recovery from an existing checkpoint. Recorded failures are
        retained; recovery continues from the first unrecorded candidate.
    delay : float, optional
        Seconds to wait between candidates processed in this call.
    on_progress : callable, optional
        Called with (completed, total, result) after saving progress. The
        initial notification uses result=None; result notifications are copies.

    Returns
    -------
    results : list of dict
        One inspection result per candidate, in input order.
    """
    if (isinstance(delay, bool) or not isinstance(delay, (int, float))
            or not math.isfinite(delay) or delay < 0):
        raise ValueError("delay must be a finite, non-negative number")
    if on_progress is not None and not callable(on_progress):
        raise ValueError("on_progress must be callable")

    body, candidates = _read_catalog(catalog_path)
    output_dir = Path(output_dir).resolve()
    if resume:
        checkpoint = _load_checkpoint(output_dir, body, candidates)
        # Completed reports need no rewrite. A missing manifest can be
        # regenerated if execution stopped while publishing the reports.
        if len(checkpoint["results"]) == len(candidates):
            if not (output_dir / "manifest.json").is_file():
                finish_inventory(output_dir, checkpoint)
            if on_progress is not None:
                on_progress(len(candidates), len(candidates), None)
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
    if on_progress is not None:
        on_progress(completed, len(candidates), None)
    for number, candidate in enumerate(candidates[completed:], start=completed + 1):
        if delay and number > completed + 1:
            time.sleep(delay)
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
            0 if result["inspection_status"] == "adapter_pending" else attempt
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
        if on_progress is not None:
            on_progress(number, len(candidates), deepcopy(result))

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
        or manifest["source"] != result["repository"]
    ):
        raise ValueError(f"Incomplete or mismatched inventory for {accession}")

    sources = [manifest["project_response"], *manifest["file_responses"]]
    responses = []
    repository = manifest["source"]
    for number, source in enumerate(sources):
        raw = _saved_path(root, source["path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != source["sha256"]:
            raise ValueError(f"Saved response checksum mismatch for {accession}")
        data = (
            _parse_massive_file_page(raw)
            if repository == "MassIVE" and number > 0 and manifest.get("backend") != "pridepy"
            else _json_response_data(raw, source["headers"])
        )
        responses.append(RepositoryResponse(
            request_url=source["request_url"],
            retrieved_at_utc=source["retrieved_at_utc"],
            raw_body=raw,
            headers=requests.structures.CaseInsensitiveDict(source["headers"]),
            data=data,
        ))

    files, reported_total = _repository_files(repository, accession, responses[0], responses[1:])
    project = responses[0].data
    if manifest.get("backend") == "pridepy":
        export = responses[1].data
        if (manifest.get("backend_version") != export["backend_version"]
                or manifest.get("repository_accession") != export["repository_accession"]
                or manifest.get("file_inventory_scope") != "provider_returned_files"
                or manifest.get("provider_record_count") != len(export["records"])
                or manifest.get("duplicate_record_count") != len(export["records"]) - len(files)):
            raise ValueError(f"Saved pridepy export metadata differs from its source: {accession}")
    else:
        if repository == "jPOST":
            link_count = len(_select_jpost_project(project, accession)["datasetFiles"])
            if (
                manifest.get("reported_link_count", link_count) != link_count
                or manifest.get("duplicate_link_count", link_count - len(files)) != link_count - len(files)
            ):
                raise ValueError(f"Saved duplicate link counts differ from their source: {accession}")
        if repository in ("jPOST", "iProX"):
            unresolved = sum(record.get("pathDecodingStatus") == "unresolved" for record in files)
            if manifest.get("unresolved_path_count", unresolved) != unresolved:
                raise ValueError(f"Saved unresolved path count differs from its source: {accession}")
        if repository == "MassIVE" and (
            manifest.get("repository_accession") != _massive_identifiers(project)["msv"]
        ):
            raise ValueError(f"Saved native project identity mismatch for {accession}")
        if repository == "jPOST" and (
            manifest.get("repository_accession") != _jpost_identifiers(
                _select_jpost_project(project, accession)
            )["jpst"]
            or manifest.get("file_inventory_scope") != "repository_reported_links"
        ):
            raise ValueError(f"Saved native project identity mismatch for {accession}")
        if repository == "iProX" and (
            manifest.get("repository_accession") != _iprox_identifiers(
                _select_proxi_project(project, accession, repository)
            )["ipx"]
            or manifest.get("file_inventory_scope") != "repository_reported_links"
        ):
            raise ValueError(f"Saved native project identity mismatch for {accession}")
        if repository == "PanoramaPublic":
            context = _panorama_project(project, accession)
            if (
                manifest.get("repository_accession") != context["native_id"]
                or manifest.get("repository_project_url") != context["project_url"]
                or manifest.get("file_inventory_scope") != "recursive_webdav_listing"
            ):
                raise ValueError(f"Saved native project identity mismatch for {accession}")
    if (
        json.loads((root / "project.json").read_bytes()) != project
        or json.loads((root / "files.json").read_bytes()) != files
        or manifest["file_count"] != len(files)
        or manifest["reported_file_count"] != reported_total
        or result["file_count"] != len(files)
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


def _summary_source_run(run_dir):
    """Validate a finished run against its checkpoint and saved raw sources."""
    names = ("input_candidates.json", "checkpoint.json", "projects.json",
             "projects.csv", "manifest.json")
    try:
        bodies = {name: (run_dir / name).read_bytes() for name in names}
        catalog_body, candidates = _read_catalog(run_dir / "input_candidates.json")
        if catalog_body != bodies["input_candidates.json"]:
            raise ValueError("Source catalog changed while being read")
        checkpoint = _load_checkpoint(run_dir, catalog_body, candidates)
        if checkpoint != json.loads(bodies["checkpoint.json"]):
            raise ValueError("Source checkpoint changed while being read")
        if checkpoint["status"] != "finished":
            raise ValueError("Source run has not finished")
        results = checkpoint["results"]
        if json.loads(bodies["projects.json"]) != results:
            raise ValueError("Project report differs from its checkpoint")

        reader = csv.DictReader(io.StringIO(bodies["projects.csv"].decode("utf-8")))
        expected = []
        for result in results:
            row = {
                column: "" if result.get(column, "") is None
                else str(result.get(column, ""))
                for column in PROJECT_COLUMNS
            }
            row["search_terms"] = json.dumps(result.get("search_terms", []), ensure_ascii=False)
            expected.append(row)
        if reader.fieldnames != PROJECT_COLUMNS or list(reader) != expected:
            raise ValueError("Project CSV differs from its checkpoint")

        manifest = json.loads(bodies["manifest.json"])
        counts = Counter(row["inspection_status"] for row in results)
        status_counts = {status: counts[status]
                         for status in ("collected", "failed", "adapter_pending")}
        if (
            manifest["schema_version"] != 1 or manifest["status"] != "finished"
            or manifest["candidate_count"] != len(candidates)
            or manifest["result_count"] != len(results)
            or manifest["source_catalog"] != checkpoint["source_catalog"]
            or manifest["status_counts"] != status_counts
            or manifest["checkpoint"] != "checkpoint.json"
            or manifest["results"] != "projects.json"
            or manifest["summary_csv"] != "projects.csv"
        ):
            raise ValueError("Source manifest differs from its checkpoint")
        source = {
            "path": run_dir.as_posix(),
            "candidate_count": len(candidates),
            "sha256": {name: hashlib.sha256(body).hexdigest()
                       for name, body in bodies.items()},
        }
        return candidates, results, source
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
        raise ValueError(f"Cannot summarize source run {run_dir}: {error}") from error


def summarize_inventories(catalog_path, run_dirs, output_dir):
    """Combine finished inventories, retaining all attempts and existing successes.

    Parameters
    ----------
    catalog_path : str or pathlib.Path
        Full discovery catalog. Source runs may contain subsets of this catalog.
    run_dirs : iterable of str or pathlib.Path
        Finished runs in chronological order. For each candidate, the last
        successful result in this order is preferred. If none succeeded, the
        last recorded outcome is used.
    output_dir : str or pathlib.Path
        New directory for the summary. All sources are checked before writing.

    Returns
    -------
    output_dir : pathlib.Path
        Summary with projects.json, projects.csv, attempts.json and manifest.json.
        Selected projects use absolute inventory and attempt paths. Rows in
        attempts.json keep the original relative paths under their source_run.
    """
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise ValueError("Choose a new summary output directory")
    catalog_body, candidates = _read_catalog(catalog_path)
    by_accession = {row["accession"]: row for row in candidates}
    runs = [Path(path).resolve() for path in run_dirs]
    if not runs or len(set(runs)) != len(runs):
        raise ValueError("Provide at least one source run, without duplicate runs")
    for run in runs:
        if output_dir == run or run in output_dir.parents:
            raise ValueError("Summary output must be outside its source runs")

    attempts = []
    sources = []
    selected = {}
    for run_index, run in enumerate(runs, 1):
        run_candidates, results, source = _summary_source_run(run)
        for candidate in run_candidates:
            if by_accession.get(candidate["accession"]) != candidate:
                raise ValueError(f"Source candidate differs from the full catalog: {candidate['accession']}")
        source["run_index"] = run_index
        sources.append(source)
        for result_index, result in enumerate(results, 1):
            attempt = deepcopy(result)
            attempt.update({
                "source_run": run.as_posix(),
                "source_run_index": run_index,
                "source_result_index": result_index,
                "summary_attempt_index": len(attempts) + 1,
                "selected": False,
            })
            attempts.append(attempt)
            accession = result["accession"]
            previous = selected.get(accession)
            if (previous is None or result["inspection_status"] == "collected"
                    or attempts[previous]["inspection_status"] != "collected"):
                selected[accession] = len(attempts) - 1

    missing = set(by_accession) - set(selected)
    if missing:
        raise ValueError(f"Source runs have no recorded outcome for {len(missing)} candidate(s)")
    attempt_counts = Counter(row["accession"] for row in attempts)
    projects = []
    for number, candidate in enumerate(candidates, 1):
        attempt = attempts[selected[candidate["accession"]]]
        attempt["selected"] = True
        result = deepcopy(attempt)
        result.pop("selected")
        result["selected_attempt_index"] = result.pop("summary_attempt_index")
        result["source_candidate_index"] = result["candidate_index"]
        result["candidate_index"] = number
        result["attempt_count"] = attempt_counts[candidate["accession"]]
        run = Path(result["source_run"])
        for field in ("inventory_path", "attempt_path"):
            if result.get(field):
                result[field] = _saved_path(run, result[field]).as_posix()
        projects.append(result)

    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "input_candidates.json").write_bytes(catalog_body)
    _write_json(output_dir / "projects.json", projects)
    _write_json(output_dir / "attempts.json", attempts)
    columns = PROJECT_COLUMNS + [
        "source_run", "source_candidate_index", "attempt_count", "selected_attempt_index",
    ]
    with (output_dir / "projects.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for result in projects:
            row = {column: result.get(column, "") for column in columns}
            row["search_terms"] = json.dumps(result.get("search_terms", []), ensure_ascii=False)
            writer.writerow(row)

    screening = _write_file_screening(projects, output_dir)
    counts = Counter(row["inspection_status"] for row in projects)
    manifest = {
        "schema_version": 1,
        "kind": "inventory_summary",
        "status": "complete",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_policy": "last_collected_in_run_order_else_last_outcome",
        "scientific_review": "unchanged",
        "file_screening": screening,
        "candidate_count": len(candidates),
        "result_count": len(projects),
        "attempt_count": len(attempts),
        "status_counts": {status: counts[status]
                          for status in ("collected", "failed", "adapter_pending")},
        "source_catalog": {
            "path": "input_candidates.json",
            "original_path": str(Path(catalog_path).resolve()),
            "sha256": hashlib.sha256(catalog_body).hexdigest(),
        },
        "source_runs": sources,
        "output_sha256": {
            name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
            for name in (
                "input_candidates.json", "projects.json", "projects.csv", "attempts.json",
                "file_screening.csv", "download_candidates.csv", "project_file_screening.json",
            )
        },
    }
    _write_json(output_dir / "manifest.json", manifest)
    return output_dir


FILE_SCREENING_COLUMNS = [
    "project_accession", "repository", "file_accession", "file_name", "relative_path",
    "reported_category", "file_size_bytes", "file_role", "classification_status",
    "possible_roles", "xl_specific_hint", "spectra_kind", "download_priority",
    "download_blockers", "evidence", "download_urls", "source_inventory_path",
    "source_file_record_index", "source_files_sha256",
]


def _write_file_screening(projects, output_dir):
    """Stream metadata screening rows; keep every project and every listed file."""
    role_counts = Counter()
    priority_counts = Counter()
    project_reports = []
    with (output_dir / "file_screening.csv").open("w", encoding="utf-8", newline="") as all_handle:
        with (output_dir / "download_candidates.csv").open("w", encoding="utf-8", newline="") as download_handle:
            all_writer = csv.DictWriter(all_handle, fieldnames=FILE_SCREENING_COLUMNS)
            download_writer = csv.DictWriter(download_handle, fieldnames=FILE_SCREENING_COLUMNS)
            all_writer.writeheader()
            download_writer.writeheader()
            for project in projects:
                report = {
                    "accession": project["accession"], "repository": project.get("repository", ""),
                    "inspection_status": project["inspection_status"],
                    "review_status": project.get("review_status", "pending"),
                    "screening_status": "inventory_unavailable", "listed_file_count": None,
                    "source_inventory_path": project.get("inventory_path", ""),
                    "role_counts": {}, "download_priority_counts": {},
                    "review_notes": ["Inventory unavailable; scientific eligibility remains unassessed"],
                }
                if project["inspection_status"] != "collected":
                    project_reports.append(report)
                    continue
                inventory = Path(project["inventory_path"])
                body = (inventory / "files.json").read_bytes()
                files = json.loads(body)
                if not isinstance(files, list) or len(files) != project["file_count"]:
                    raise ValueError(f"File inventory count changed for {project['accession']}")
                digest = hashlib.sha256(body).hexdigest()
                per_role = Counter()
                per_priority = Counter()
                for index, record in enumerate(files):
                    classification = classify_inventory_file(record)
                    category = record.get("fileCategory", "")
                    if isinstance(category, dict):
                        category = category.get("value", "")
                    row = {
                        "project_accession": project["accession"],
                        "repository": project.get("repository", ""),
                        "file_accession": record["accession"],
                        "file_name": record.get("fileName", ""),
                        "relative_path": record.get("relativePath", ""),
                        "reported_category": category,
                        "file_size_bytes": record.get("fileSizeBytes"),
                        **classification,
                        "source_inventory_path": inventory.as_posix(),
                        "source_file_record_index": index,
                        "source_files_sha256": digest,
                    }
                    for field in ("possible_roles", "download_blockers", "evidence", "download_urls"):
                        row[field] = json.dumps(row[field], ensure_ascii=False)
                    all_writer.writerow(row)
                    if classification["download_priority"] == "first":
                        download_writer.writerow(row)
                    per_role[classification["file_role"]] += 1
                    per_priority[classification["download_priority"]] += 1
                notes = ["Metadata screening only; file contents and scientific eligibility remain unassessed"]
                if not per_role["identification_result_candidate"]:
                    notes.append("No identification result candidate identified from metadata; review required")
                if not per_role["search_fasta_candidate"]:
                    notes.append("No search FASTA candidate identified from metadata; review required")
                report.update({
                    "screening_status": "metadata_screened", "listed_file_count": len(files),
                    "source_files_sha256": digest,
                    "role_counts": dict(per_role), "download_priority_counts": dict(per_priority),
                    "review_notes": notes,
                })
                role_counts.update(per_role)
                priority_counts.update(per_priority)
                project_reports.append(report)
    _write_json(output_dir / "project_file_screening.json", project_reports)
    return {
        "rule_version": FILE_SCREENING_RULE_VERSION,
        "basis": "metadata_only", "downloads_performed": False,
        "file_count": sum(role_counts.values()),
        "download_candidate_count": priority_counts["first"],
        "role_counts": dict(role_counts), "download_priority_counts": dict(priority_counts),
        "files": "file_screening.csv", "download_candidates": "download_candidates.csv",
        "projects": "project_file_screening.json",
    }


def _read_content_plan(plan_path, output_dir):
    """Bind every planned transfer to the saved candidate and inventory records."""
    body = Path(plan_path).read_bytes()
    plan = json.loads(body)
    try:
        if plan["schema_version"] != 1 or plan["kind"] != "xlms_file_content_pilot":
            raise ValueError("Unsupported content plan")
        files = plan["files"]
        if not isinstance(files, list) or not files or plan["counts"]["files"] != len(files):
            raise ValueError("Content plan file count is invalid")
        limits = plan["proposed_download_limits"]
        for key in ("per_file_bytes", "total_payload_bytes"):
            if type(limits[key]) is not int or limits[key] <= 0:
                raise ValueError("Download limits must be positive integer byte counts")
        with Path(plan["source"]["server_candidate_csv"]).open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != plan["source"]["candidate_record_count"]:
            raise ValueError("Candidate table record count changed")
        ids, destinations, cache = set(), set(), {}
        for item in files:
            identifier = item["pilot_file_id"]
            if not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z0-9_-]+", identifier) or identifier in ids:
                raise ValueError("Pilot file IDs must be unique simple identifiers")
            ids.add(identifier)
            relative = Path(item["local_relative_path"])
            if not relative.parts or relative.parts[0] != "files":
                raise ValueError("Pilot payload paths must be within files/")
            destination = _saved_path(output_dir, item["local_relative_path"])
            if destination in destinations:
                raise ValueError("Pilot payload destinations must be unique")
            destinations.add(destination)
            number = item["source_candidate_record_number"]
            if type(number) is not int or not 1 <= number <= len(rows):
                raise ValueError("Invalid source candidate record number")
            row = rows[number - 1]
            digest = hashlib.sha256(json.dumps(row, ensure_ascii=False, sort_keys=True,
                                               separators=(",", ":")).encode("utf-8")).hexdigest()
            if row != item["source_candidate_record"] or digest != item["source_candidate_record_sha256"]:
                raise ValueError("Planned source candidate record changed")
            urls = json.loads(row["download_urls"])
            size = int(row["file_size_bytes"]) if row["file_size_bytes"] else None
            fields = (("project_accession", "project_accession"), ("repository", "repository"),
                      ("file_accession", "file_accession"), ("reported_file_name", "file_name"),
                      ("source_inventory_path", "source_inventory_path"))
            if (any(item[left] != row[right] for left, right in fields)
                    or item["download_urls"] != urls or item["reported_size_bytes"] != size
                    or row["download_priority"] != "first" or json.loads(row["download_blockers"])):
                raise ValueError("Planned transfer differs from its source candidate")
            if not urls or any(not isinstance(url, str) or re.match(r"^(?:https?|ftp)://", url) is None for url in urls):
                raise ValueError("Pilot transfers require an HTTP, HTTPS or FTP URL")
            if size is not None and (type(size) is not int or size <= 0):
                raise ValueError("Reported file size must be positive or unknown")
            expected = item.get("expected_download_file_sha256")
            if expected is not None and (not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected)):
                raise ValueError("Expected file SHA256 is invalid")
            inventory = Path(item["source_inventory_path"]).resolve()
            if inventory == output_dir or inventory in output_dir.parents:
                raise ValueError("Pilot output must be outside source inventories")
            if inventory not in cache:
                raw = (inventory / "files.json").read_bytes()
                cache[inventory] = (hashlib.sha256(raw).hexdigest(), json.loads(raw))
            inventory_hash, records = cache[inventory]
            index = item["source_file_record_index"]
            if (inventory_hash != item["inventory_files_json_sha256"]
                    or inventory_hash != row["source_files_sha256"]
                    or type(index) is not int or not 0 <= index < len(records)
                    or str(index) != row["source_file_record_index"]):
                raise ValueError("Source inventory checksum or record index changed")
            record = records[index]
            reported_urls = [entry.get("value") for entry in record.get("publicFileLocations", [])
                             if isinstance(entry, dict)]
            if (record["accession"] != item["file_accession"]
                    or record.get("fileSizeBytes") != size
                    or record.get("fileName") != item["reported_file_name"]
                    or item["project_accession"] not in record["projectAccessions"]
                    or any(url not in reported_urls for url in urls)
                    or record.get("pathDecodingStatus") == "unresolved"):
                raise ValueError("Source file identity or reported download link changed")
            if item["pilot_role"] not in (
                "result_table_candidate", "fasta_candidate", "possible_library_control", "possible_instructions_control",
            ):
                raise ValueError("Unsupported pilot content role")
    except (KeyError, TypeError, AttributeError, IndexError) as error:
        raise ValueError(f"Invalid content plan or source record: {error}") from error
    return body, plan


def _save_content_checkpoint(output_dir, checkpoint):
    temporary = output_dir / ".content_checkpoint.json.tmp"
    _write_json(temporary, checkpoint)
    temporary.replace(output_dir / "content_checkpoint.json")


def _validate_pilot_result(result, item, output_dir):
    if result["pilot_file_id"] != item["pilot_file_id"]:
        raise ValueError("Saved content results do not match the plan order")
    for key in ("charged_payload_bytes", "confirmed_payload_bytes"):
        if type(result[key]) is not int or result[key] < 0:
            raise ValueError("Invalid saved payload accounting")
    if result["charged_payload_bytes"] < result["confirmed_payload_bytes"]:
        raise ValueError("Saved payload budget understates confirmed bytes")
    if result.get("accounting_mode") == "planned_file_sizes":
        expected = item["reported_size_bytes"]
        charged = result["charged_payload_bytes"]
        if result["download_status"] == "downloaded":
            if charged != expected or result["confirmed_payload_bytes"] != expected:
                raise ValueError("Saved pridepy file size differs from its plan")
        elif result["download_status"] in ("skipped_unknown_size", "skipped_file_limit", "budget_exhausted"):
            if charged or result["confirmed_payload_bytes"]:
                raise ValueError("A skipped file must not consume the planned budget")
        elif charged not in (0, expected) or result["confirmed_payload_bytes"]:
            raise ValueError("Saved failed pridepy transfer has invalid plan accounting")
    if result["download_status"] == "downloaded":
        path = _saved_path(output_dir, item["local_relative_path"])
        raw = path.read_bytes()
        if (str(path) != result["local_path"] or len(raw) != result["confirmed_payload_bytes"]
                or hashlib.sha256(raw).hexdigest() != result["file_sha256"]):
            raise ValueError("Saved downloaded file differs from its recorded checksum")
    elif result["download_status"] not in ("failed", "budget_exhausted", "skipped_unknown_size", "skipped_file_limit"):
        raise ValueError("Invalid saved download status")


def _finish_content_pilot(output_dir, checkpoint):
    results = checkpoint["results"]
    _write_json(output_dir / "results.json", results)
    columns = ["pilot_file_id", "project_accession", "repository", "reported_file_name", "pilot_role", "download_status",
               "content_status", "confirmed_payload_bytes", "charged_payload_bytes", "file_sha256", "local_path",
               "error_type", "error_message"]
    with (output_dir / "results.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for result in results:
            row = {key: result.get(key, "") for key in columns}
            row["content_status"] = result["content_probe"]["status"]
            writer.writerow(row)
    checkpoint["status"] = "finished"
    _save_content_checkpoint(output_dir, checkpoint)
    _write_json(output_dir / "manifest.json", {
        "schema_version": 1, "kind": "xlms_file_content_pilot", "status": "finished",
        "scientific_review": "unchanged", "file_count": len(results),
        "plan_sha256": checkpoint["plan_sha256"], "limits": checkpoint["limits"],
        "charged_payload_bytes": checkpoint["charged_payload_bytes"],
        "confirmed_completed_payload_bytes": sum(row["confirmed_payload_bytes"] for row in results),
        "accounting_mode": checkpoint.get("accounting_mode", "bounded_payload_reads"),
        "budget_accounting": (
            "Known file sizes are reserved before pridepy calls; failed/interrupted calls retain their reservation. "
            "These are plan limits, not hard network-byte limits. Unknown sizes are skipped."
            if checkpoint.get("accounting_mode") == "planned_file_sizes" else
            "Interrupted or uncertain reads retain their reservation; confirmed reads charge actual bytes."
        ),
        "download_status_counts": dict(Counter(row["download_status"] for row in results)),
        "content_status_counts": dict(Counter(row["content_probe"]["status"] for row in results)),
        "interrupted_attempt_count": len(checkpoint["interrupted_attempts"]),
        "output_sha256": {name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
                          for name in ("input_plan.json", "content_checkpoint.json", "results.json", "results.csv")},
    })


def inspect_content_plan(plan_path, output_dir, *, resume=False, delay=0, on_progress=None):
    """Execute a bounded pilot plan, checkpointing every file outcome.

    Resume reuses verified outcomes and retries an interrupted active transfer.
    Its original reservation remains charged when consumption cannot be known.
    Failed completed outcomes are retained rather than silently retried.
    """
    if isinstance(delay, bool) or not isinstance(delay, (int, float)) or not math.isfinite(delay) or delay < 0:
        raise ValueError("delay must be finite and non-negative")
    if on_progress is not None and not callable(on_progress):
        raise ValueError("on_progress must be callable")
    output_dir = Path(output_dir).resolve()
    body, plan = _read_content_plan(plan_path, output_dir)
    plan_hash = hashlib.sha256(body).hexdigest()
    files, limits = plan["files"], plan["proposed_download_limits"]
    if resume:
        try:
            checkpoint = json.loads((output_dir / "content_checkpoint.json").read_bytes())
            if (checkpoint["schema_version"] not in (1, 2) or checkpoint["plan_sha256"] != plan_hash
                    or hashlib.sha256((output_dir / "input_plan.json").read_bytes()).hexdigest() != plan_hash
                    or checkpoint["limits"] != limits or checkpoint["status"] not in ("running", "finished")
                    or not isinstance(checkpoint["results"], list) or len(checkpoint["results"]) > len(files)):
                raise ValueError("Saved content checkpoint differs from the plan")
            if checkpoint["schema_version"] == 2 and checkpoint.get("accounting_mode") != "planned_file_sizes":
                raise ValueError("Saved pridepy plan accounting mode changed")
            if (checkpoint["schema_version"] == 1 and checkpoint["status"] != "finished"
                    and len(checkpoint["results"]) != len(files)):
                raise ValueError("An unfinished legacy download run cannot switch download backends; use a new plan/run")
            for result, item in zip(checkpoint["results"], files):
                _validate_pilot_result(result, item, output_dir)
                saved = json.loads((output_dir / "outcomes" / (item["pilot_file_id"] + ".json")).read_bytes())
                if saved != result:
                    raise ValueError("Saved per-file outcome differs from its checkpoint")
            if any(type(row["charged_payload_bytes"]) is not int or row["charged_payload_bytes"] < 0
                   or row["pilot_file_id"] not in {item["pilot_file_id"] for item in files}
                   for row in checkpoint["interrupted_attempts"]):
                raise ValueError("Invalid interrupted transfer accounting")
            charged = sum(row["charged_payload_bytes"] for row in checkpoint["results"])
            charged += sum(row["charged_payload_bytes"] for row in checkpoint["interrupted_attempts"])
            active = checkpoint["active"]
            if active:
                if (len(checkpoint["results"]) >= len(files)
                        or active["pilot_file_id"] != files[len(checkpoint["results"])]["pilot_file_id"]
                        or type(active["reservation_bytes"]) is not int or active["reservation_bytes"] <= 0):
                    raise ValueError("Invalid active transfer reservation")
                charged += active["reservation_bytes"]
            if (type(checkpoint["charged_payload_bytes"]) is not int
                    or charged != checkpoint["charged_payload_bytes"]
                    or not 0 <= charged <= limits["total_payload_bytes"]):
                raise ValueError("Saved transfer budget is invalid")
            if checkpoint["status"] == "finished" and (active or len(checkpoint["results"]) != len(files)):
                raise ValueError("Finished pilot has missing outcomes")
        except (OSError, KeyError, TypeError, ValueError) as error:
            raise ValueError(f"Cannot resume content pilot: {error}") from error
        if len(checkpoint["results"]) == len(files):
            if not (output_dir / "manifest.json").is_file():
                _finish_content_pilot(output_dir, checkpoint)
            else:
                try:
                    manifest = json.loads((output_dir / "manifest.json").read_bytes())
                    names = {"input_plan.json", "content_checkpoint.json", "results.json", "results.csv"}
                    if (manifest["kind"] != "xlms_file_content_pilot" or manifest["status"] != "finished"
                            or manifest["plan_sha256"] != plan_hash
                            or set(manifest["output_sha256"]) != names
                            or json.loads((output_dir / "results.json").read_bytes()) != checkpoint["results"]
                            or any(hashlib.sha256((output_dir / name).read_bytes()).hexdigest() != digest
                                   for name, digest in manifest["output_sha256"].items())):
                        raise ValueError("Completed pilot reports changed")
                except (KeyError, TypeError, OSError, ValueError) as error:
                    raise ValueError(f"Cannot resume content pilot: {error}") from error
            if on_progress:
                on_progress(len(files), len(files), None)
            return checkpoint["results"]
        if active:
            item = files[len(checkpoint["results"])]
            outcome_path = output_dir / "outcomes" / (item["pilot_file_id"] + ".json")
            if outcome_path.is_file():
                result = json.loads(outcome_path.read_bytes())
                _validate_pilot_result(result, item, output_dir)
                if result["charged_payload_bytes"] > active["reservation_bytes"]:
                    raise ValueError("Recovered outcome exceeds its reservation")
                checkpoint["charged_payload_bytes"] -= active["reservation_bytes"] - result["charged_payload_bytes"]
                checkpoint["results"].append(result)
            else:
                event = {"pilot_file_id": item["pilot_file_id"],
                         "charged_payload_bytes": active["reservation_bytes"],
                         "reason": "Interrupted transfer with unknown consumption", "preserved_paths": []}
                checkpoint["interrupted_attempts"].append(event)
                destination = _saved_path(output_dir, item["local_relative_path"])
                for path in (destination, destination.with_name(destination.name + ".part")):
                    if path.exists():
                        archive = output_dir / "interrupted" / str(len(checkpoint["interrupted_attempts"])) / path.name
                        archive.parent.mkdir(parents=True, exist_ok=True)
                        path.replace(archive)
                        event["preserved_paths"].append(archive.relative_to(output_dir).as_posix())
            checkpoint["active"] = None
            _save_content_checkpoint(output_dir, checkpoint)
    else:
        output_dir.mkdir(parents=True, exist_ok=False)
        (output_dir / "outcomes").mkdir()
        (output_dir / "input_plan.json").write_bytes(body)
        checkpoint = {
            "schema_version": 2, "status": "running", "plan_sha256": plan_hash,
            "accounting_mode": "planned_file_sizes",
            "limits": limits, "charged_payload_bytes": 0, "active": None,
            "results": [], "interrupted_attempts": [],
        }
        _save_content_checkpoint(output_dir, checkpoint)
    completed = len(checkpoint["results"])
    if on_progress:
        on_progress(completed, len(files), None)
    for number, item in enumerate(files[completed:], completed + 1):
        if delay and number > completed + 1:
            time.sleep(delay)
        remaining = limits["total_payload_bytes"] - checkpoint["charged_payload_bytes"]
        expected = item["reported_size_bytes"]
        reservation = expected if type(expected) is int and 0 < expected <= min(remaining, limits["per_file_bytes"]) else 0
        destination = _saved_path(output_dir, item["local_relative_path"])
        if destination.exists() or destination.with_name(destination.name + ".part").exists():
            raise ValueError("Unrecorded payload already exists; refusing to overwrite it")
        if reservation:
            checkpoint["active"] = {"pilot_file_id": item["pilot_file_id"], "reservation_bytes": reservation}
            checkpoint["charged_payload_bytes"] += reservation
            _save_content_checkpoint(output_dir, checkpoint)
        result = download_pilot_file(item, destination, limits["per_file_bytes"], reservation)
        if reservation:
            checkpoint["charged_payload_bytes"] -= reservation - result["charged_payload_bytes"]
        outcome = output_dir / "outcomes" / (item["pilot_file_id"] + ".json")
        _write_json(outcome.with_suffix(".json.tmp"), result)
        outcome.with_suffix(".json.tmp").replace(outcome)
        checkpoint["results"].append(result)
        checkpoint["active"] = None
        _save_content_checkpoint(output_dir, checkpoint)
        if on_progress:
            on_progress(number, len(files), deepcopy(result))
    _finish_content_pilot(output_dir, checkpoint)
    return checkpoint["results"]


def prepare_content_retry(previous_run, output_plan):
    """Create a failure-only retry plan within the original cumulative budget.

    Finished parent outcomes, downloaded payloads and source bindings are
    verified without requesting files. The parent run is kept unchanged.
    """
    previous_run = Path(previous_run).resolve()
    output_plan = Path(output_plan).resolve()
    if output_plan.exists() or previous_run in output_plan.parents:
        raise ValueError("Choose a new retry plan outside the parent run")
    plan_path = previous_run / "input_plan.json"
    try:
        parent_body = plan_path.read_bytes()
        original = json.loads(parent_body)
        manifest_body = (previous_run / "manifest.json").read_bytes()
        manifest = json.loads(manifest_body)
        checkpoint = json.loads((previous_run / "content_checkpoint.json").read_bytes())
        if (manifest["kind"] != "xlms_file_content_pilot" or manifest["status"] != "finished"
                or checkpoint["status"] != "finished" or checkpoint["active"] is not None
                or len(checkpoint["results"]) != len(original["files"])):
            raise ValueError("The parent content pilot must be finished")
        # Finished resume validates the parent files/reports and makes no requests.
        results = inspect_content_plan(plan_path, previous_run, resume=True)
        spent = checkpoint["charged_payload_bytes"]
        remaining = original["proposed_download_limits"]["total_payload_bytes"] - spent
        if remaining <= 0:
            raise ValueError("The original payload budget is exhausted")
        failed = {row["pilot_file_id"] for row in results
                  if row["download_status"] in ("failed", "budget_exhausted")}
        if not failed:
            raise ValueError("The parent run has no failed downloads to retry")
        retry = deepcopy(original)
        retry["files"] = [item for item in retry["files"] if item["pilot_file_id"] in failed]
        retry["proposed_download_limits"]["total_payload_bytes"] = remaining
        root_budget = original.get("retry_context", {}).get(
            "root_total_payload_bytes", original["proposed_download_limits"]["total_payload_bytes"]
        )
        retry["retry_context"] = {
            "parent_run": previous_run.as_posix(),
            "parent_plan_sha256": hashlib.sha256(parent_body).hexdigest(),
            "parent_manifest_sha256": hashlib.sha256(manifest_body).hexdigest(),
            "parent_charged_payload_bytes": spent,
            "root_total_payload_bytes": root_budget,
            "cumulative_charged_before_retry": root_budget - remaining,
        }
        retry["selection_policy"] = {
            "mode": "retry_failed_downloads_only",
            "repository_coverage": sorted({item["repository"] for item in retry["files"]}),
            "source_bindings": "Original per-file candidate records and inventory checksums remain unchanged",
            "scientific_review": "unchanged",
        }
        retry["counts"] = {
            "projects": len({item["project_accession"] for item in retry["files"]}),
            "files": len(retry["files"]),
            "roles": dict(Counter(item["pilot_role"] for item in retry["files"])),
            "reported_size_total_bytes": sum(item["reported_size_bytes"] or 0 for item in retry["files"]),
            "files_with_unknown_size": sum(item["reported_size_bytes"] is None for item in retry["files"]),
        }
        retry["status"], retry["downloads_performed"] = "planned", False
        text = json.dumps(retry, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    except (KeyError, TypeError, OSError) as error:
        raise ValueError(f"Cannot prepare content retry: {error}") from error
    output_plan.parent.mkdir(parents=True, exist_ok=True)
    with output_plan.open("x", encoding="utf-8") as handle:
        handle.write(text)
    return output_plan
