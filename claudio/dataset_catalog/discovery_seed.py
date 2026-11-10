
"""Explicit migration of original raw discovery snapshots into emission history."""

from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import tempfile
from urllib.parse import parse_qs, urlsplit

from . import _write_json
from .discovery_history import DiscoveryHistory, default_history_directory
from .proteomecentral import CatalogPage, merge_search_results
from .publication_lookup import PublicationLookup
from .published_discovery import _canonical_accession, _save_search, qualify_dataset_candidate


def _body(directory, name, limit):
    if not isinstance(name, str) or not name or "\\" in name:
        raise ValueError("Invalid legacy snapshot path")
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(
            "Legacy snapshot paths must stay inside the original snapshot"
        )
    path = directory / relative
    if (
        path.is_symlink()
        or not path.resolve().is_relative_to(directory)
        or not path.is_file()
    ):
        raise ValueError("Missing or unsafe legacy snapshot path")
    if path.stat().st_size > limit:
        raise ValueError(
            "Legacy snapshot file exceeds the validation limit"
        )
    return path.read_bytes()


def _timestamp(value, label):
    if not isinstance(value, str):
        raise ValueError("Missing legacy " + label + " date")
    try:
        parsed = datetime.fromisoformat(
            value.replace("Z", "+00:00")
        )
    except ValueError as error:
        raise ValueError(
            "Invalid legacy " + label + " date"
        ) from error
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(
            "Legacy " + label + " date must include a timezone"
        )
    return parsed.astimezone(timezone.utc)


def _load_original(directory):
    bodies = {
        "manifest.json": _body(
            directory, "manifest.json", 2 * 1024 * 1024
        ),
        "candidates.json": _body(
            directory, "candidates.json", 64 * 1024 * 1024
        ),
    }
    try:
        manifest = json.loads(bodies["manifest.json"])
        candidates = json.loads(bodies["candidates.json"])
    except (ValueError, UnicodeError) as error:
        raise ValueError(
            "Invalid original discovery JSON"
        ) from error

    if (
        not isinstance(manifest, dict)
        or manifest.get("schema_version") != 1
        or manifest.get("source") != "ProteomeCentral"
        or manifest.get("status") != "complete"
        or not isinstance(manifest.get("searches"), list)
        or not manifest["searches"]
    ):
        raise ValueError(
            "A complete original ProteomeCentral discovery snapshot is required"
        )

    created = _timestamp(
        manifest.get("created_at_utc"), "run"
    )
    if created > datetime.now(timezone.utc):
        raise ValueError(
            "Legacy source run date is in the future"
        )

    searches = {}
    sizes = {}
    all_paths = set(bodies)
    total_bytes = sum(len(raw) for raw in bodies.values())

    for search in manifest["searches"]:
        if not isinstance(search, dict):
            raise ValueError(
                "Invalid original discovery search receipt"
            )
        term = search.get("search_term")
        if (
            not isinstance(term, str)
            or not term.strip()
            or term in searches
        ):
            raise ValueError(
                "Missing or repeated original search term"
            )
        entries = search.get("pages")
        if not isinstance(entries, list) or not entries:
            raise ValueError(
                "Incomplete original discovery source pages"
            )

        pages = []
        for entry in entries:
            if not isinstance(entry, dict):
                raise ValueError(
                    "Invalid original source page receipt"
                )
            name = entry.get("path")
            if not isinstance(name, str) or name in all_paths:
                raise ValueError(
                    "Invalid or repeated original source path"
                )
            all_paths.add(name)

            remaining = 256 * 1024 * 1024 - total_bytes
            raw = _body(
                directory,
                name,
                min(32 * 1024 * 1024, remaining),
            )
            total_bytes += len(raw)
            digest = entry.get("sha256")
            if (
                not isinstance(digest, str)
                or re.fullmatch("[0-9a-f]{64}", digest) is None
                or hashlib.sha256(raw).hexdigest() != digest
            ):
                raise ValueError(
                    "Original discovery raw response checksum mismatch"
                )

            try:
                data = json.loads(raw)
            except (ValueError, UnicodeError) as error:
                raise ValueError(
                    "Invalid original discovery response JSON"
                ) from error

            if (
                not isinstance(data, dict)
                or not isinstance(data.get("result_set"), dict)
                or not isinstance(data.get("datasets"), list)
                or entry.get("query") != data.get("query")
            ):
                raise ValueError(
                    "Original page query receipt differs from its raw response"
                )

            retrieved = _timestamp(
                entry.get("retrieved_at_utc"), "retrieval"
            )
            if retrieved > created:
                raise ValueError(
                    "Original response retrieval date is after the discovery run date"
                )

            url = entry.get("request_url")
            if not isinstance(url, str):
                raise ValueError(
                    "Original page has no source URL"
                )
            try:
                values = parse_qs(
                    urlsplit(url).query
                ).get("pageSize", [])
                size = int(values[0]) if len(values) == 1 else 0
            except (TypeError, ValueError, IndexError) as error:
                raise ValueError(
                    "Invalid original search page size"
                ) from error

            if size < 1 or (
                term in sizes and sizes[term] != size
            ):
                raise ValueError(
                    "Inconsistent original search page size"
                )

            sizes[term] = size
            bodies[name] = raw
            pages.append(CatalogPage(
                url, entry["retrieved_at_utc"], raw, data
            ))

        if (
            search.get("reported_total")
            != pages[0].data.get("result_set", {}).get(
                "n_available_rows"
            )
            or search.get("downloaded_records")
            != sum(
                len(page.data.get("datasets", []))
                for page in pages
            )
        ):
            raise ValueError(
                "Original discovery search count receipt mismatch"
            )
        searches[term] = pages

    if sum(len(raw) for raw in bodies.values()) > 256 * 1024 * 1024:
        raise ValueError(
            "Original discovery source package exceeds the validation limit"
        )

    try:
        recomputed = merge_search_results(searches)
    except (KeyError, TypeError, AttributeError) as error:
        raise ValueError(
            "Invalid original discovery rows"
        ) from error

    if (
        not isinstance(candidates, list)
        or candidates != recomputed
    ):
        raise ValueError(
            "Original candidates.json does not equal the recomputed raw search merge"
        )
    if (
        type(manifest.get("candidate_count")) is not int
        or manifest["candidate_count"] != len(candidates)
    ):
        raise ValueError(
            "Original discovery candidate count mismatch"
        )

    identifiers = [
        _canonical_accession(candidate.get("accession"))
        for candidate in candidates
    ]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError(
            "Duplicate canonical accessions in original discovery"
        )

    return {
        "manifest": manifest,
        "candidates": candidates,
        "searches": searches,
        "page_sizes": sizes,
        "created": created,
        "bodies": bodies,
    }


def seed_legacy_discovery(
    source_run,
    output_dir,
    *,
    history_dir=None,
    publication_request_limit=1000,
    progress=None,
):
    """Register only old candidates published by the old UTC run date.

    Qualification has calendar-date granularity. Current metadata can
    verify an old explicit article identifier, but cannot add a new
    dataset-paper link to the original snapshot. Seed outputs are
    historical registration receipts, not a new processing queue.
    """
    output = Path(output_dir).expanduser().resolve()
    source = Path(source_run).expanduser().resolve()

    if output.exists():
        raise FileExistsError(
            "Choose a new legacy seeding output directory"
        )
    if (
        output == source
        or source in output.parents
        or output in source.parents
    ):
        raise ValueError(
            "Legacy source and seed output directories must be separate"
        )
    if (
        type(publication_request_limit) is not int
        or not 1 <= publication_request_limit <= 10000
    ):
        raise ValueError(
            "publication_request_limit must be between 1 and 10000"
        )

    history_path = (
        history_dir
        if history_dir is not None
        else default_history_directory()
    )

    with DiscoveryHistory(history_path) as history:
        seen = history.seen_accessions()
        original = _load_original(source)
        as_of = original["created"].date().isoformat()
        run = history.prepare_run(output)

        output.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(
            prefix=".claudio_legacy_seed_" + run["run_id"] + "_",
            dir=output.parent,
        ))

        retained = staging / "legacy_source"
        for name, raw in original["bodies"].items():
            path = retained / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)

        validation = staging / "source_validation"
        (validation / "raw").mkdir(parents=True)
        validated = []
        for number, (term, pages) in enumerate(
            original["searches"].items(), 1
        ):
            validated.append(_save_search(
                validation,
                term,
                pages,
                number,
                original["page_sizes"][term],
            ))
        _write_json(
            validation / "searches.json", validated
        )

        reviewed = []
        with PublicationLookup(
            staging / "publication_evidence",
            max_requests=publication_request_limit,
        ) as lookup:
            for index, candidate in enumerate(
                original["candidates"], 1
            ):
                accession = _canonical_accession(
                    candidate["accession"]
                )

                if progress and (
                    index == 1
                    or index % 10 == 0
                    or index == len(original["candidates"])
                ):
                    progress(
                        f"Historical publication review: "
                        f"{index}/{len(original['candidates'])}; "
                        f"HTTP attempts: {lookup.request_count}"
                    )

                if accession in seen:
                    annotated = deepcopy(candidate)
                    annotated.update(
                        accession=accession,
                        publication_status="unverified",
                        papers=[],
                        seed_decision="excluded_seen",
                        publication_review={
                            "reason": "already_emitted_in_history",
                            "current_publications_rechecked": False,
                        },
                    )
                else:
                    annotated = qualify_dataset_candidate(
                        candidate, lookup, as_of=as_of
                    )
                    annotated["seed_decision"] = (
                        "eligible"
                        if annotated["publication_status"] == "published"
                        else "not_confirmed_published_by_legacy_date"
                    )

                reviewed.append(annotated)

            transport = {
                "request_limit": publication_request_limit,
                "http_attempts": lookup.request_count,
                "captured_body_bytes": lookup.captured_body_bytes,
            }

        for name, raw in original["bodies"].items():
            if _body(source, name, len(raw) + 1) != raw:
                raise ValueError(
                    "Original discovery source changed during seeding"
                )

        source_manifest_sha = hashlib.sha256(
            original["bodies"]["manifest.json"]
        ).hexdigest()
        source_candidates_sha = hashlib.sha256(
            original["bodies"]["candidates.json"]
        ).hexdigest()

        def write_outputs(destination, selected):
            for path in sorted(staging.iterdir()):
                path.rename(destination / path.name)

            emitted = {
                candidate["accession"]
                for candidate in selected
            }
            audit = deepcopy(reviewed)
            for candidate in audit:
                if candidate["publication_status"] == "published":
                    candidate["seed_decision"] = (
                        "seeded"
                        if candidate["accession"] in emitted
                        else "excluded_seen"
                    )

            report = {
                "kind": "legacy_discovery_history_seed_report",
                "run_id": run["run_id"],
                "namespace_id": run["namespace_id"],
                "source_directory": str(source),
                "source_manifest_sha256": source_manifest_sha,
                "source_candidates_sha256": source_candidates_sha,
                "source_created_at_utc":
                    original["manifest"]["created_at_utc"],
                "publication_as_of": as_of,
                "date_precision": "calendar_day",
                "source_candidates": len(reviewed),
                "seeded_count": len(selected),
                "seeded_accessions": [
                    candidate["accession"]
                    for candidate in selected
                ],
                "not_seeded_accessions": [
                    candidate["accession"]
                    for candidate in audit
                    if candidate["accession"] not in emitted
                ],
                "publication_status_counts": dict(Counter(
                    candidate["publication_status"]
                    for candidate in audit
                )),
                "seed_decision_counts": dict(Counter(
                    candidate["seed_decision"]
                    for candidate in audit
                )),
                "publication_transport": transport,
                "candidate_output_purpose":
                    "historical_history_registration_only",
            }

            _write_json(
                destination / "seed_audit.json", audit
            )
            _write_json(
                destination / "candidates.json", selected
            )
            _write_json(
                destination / "seed_report.json", report
            )

            return {
                "source": "ProteomeCentral",
                "operation": "legacy_history_seed",
                "candidate_output_purpose":
                    "historical_history_registration_only",
                "source_manifest_sha256": source_manifest_sha,
                "source_candidates_sha256": source_candidates_sha,
                "source_created_at_utc":
                    original["manifest"]["created_at_utc"],
                "publication_as_of": as_of,
                "source_candidate_count": len(reviewed),
                "input_exported": False,
                "assessment_performed": False,
            }

        history.publish(
            run["run_id"], reviewed, write_outputs
        )
        staging.rmdir()
        return output
