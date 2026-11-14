"""Load, verify, and package dataset-scoped publication metadata offline."""

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

from .proteomecentral import CatalogPage, merge_search_results, validate_search_page
from .publications import extract_publication_links
from .publication_lookup import replay_publication_lookup
from .publication_verification import _as_of_date

MAX_PACKAGE_BYTES = 16 * 1024 * 1024
SOURCE_MANIFEST = "publication_sources/discovery_manifest.json"
HTTP_ROOT = "publication_sources/http"
REGISTRY_ROOT = "publication_sources/registry/"


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _accession(value):
    _require(isinstance(value, str), "A dataset accession is required")
    value = value.strip().upper()
    _require(
        re.fullmatch(r"PXD[0-9]+", value) is not None,
        "Invalid dataset accession",
    )
    return value


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _encode(value):
    return (
        json.dumps(
            value, indent=2, ensure_ascii=False, allow_nan=False
        ) + "\n"
    ).encode("utf-8")


def _read(root, name, limit=64 * 1024 * 1024):
    _require(
        isinstance(name, str) and name and "\\" not in name,
        "Invalid metadata source path",
    )
    path = Path(name)
    _require(
        not path.is_absolute() and ".." not in path.parts,
        "Metadata source path escapes its directory",
    )
    path = root / path
    _require(
        not path.is_symlink()
        and path.resolve().is_relative_to(root)
        and path.is_file(),
        "Metadata source is missing or outside its directory",
    )
    _require(
        path.stat().st_size <= limit,
        "Metadata source exceeds its size limit",
    )
    return path.read_bytes()


def _json(raw):
    try:
        return json.loads(raw)
    except (ValueError, UnicodeError) as error:
        raise ValueError("Invalid publication metadata JSON") from error


def _manifest(root):
    raw = _read(root, "manifest.json", 8 * 1024 * 1024)
    value = _json(raw)
    _require(
        isinstance(value, dict)
        and value.get("status") == "complete"
        and isinstance(value.get("output_sha256"), dict),
        "A complete checksummed manifest is required",
    )
    return value, raw


def _verified(root, name, manifest):
    expected = manifest["output_sha256"].get(name)
    _require(
        isinstance(expected, str)
        and re.fullmatch("[0-9a-f]{64}", expected) is not None,
        "Metadata source is not pinned by its manifest",
    )
    raw = _read(root, name)
    _require(
        _sha(raw) == expected,
        "Publication metadata source checksum mismatch",
    )
    return raw


def _discovery_policy(manifest, accession):
    _require(
        manifest.get("schema_version") == 2
        and manifest.get("kind") == "published_incremental_discovery"
        and manifest.get("publication_policy") == "confirmed_published_only"
        and manifest.get("source") == "ProteomeCentral"
        and manifest.get("operation") != "legacy_history_seed",
        "Metadata requires a verified published discovery output",
    )
    emitted = manifest.get("emitted_accessions")
    _require(
        isinstance(emitted, list) and accession in emitted,
        "Dataset is not an emitted candidate in the discovery",
    )
    _require(
        type(manifest.get("candidate_count")) is int
        and manifest["candidate_count"] == len(emitted)
        and all(isinstance(value, str) for value in emitted)
        and len(set(emitted)) == len(emitted),
        "Invalid discovery candidate receipt",
    )
    as_of = _as_of_date(manifest.get("publication_as_of"))
    try:
        created = datetime.fromisoformat(
            manifest["created_at_utc"].replace("Z", "+00:00")
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("Invalid discovery creation date") from error
    _require(
        created.tzinfo is not None,
        "Discovery creation date requires a timezone",
    )
    _require(
        as_of <= created.astimezone(timezone.utc).date()
        and as_of <= datetime.now(timezone.utc).date(),
        "Future publication qualification date",
    )
    return as_of.isoformat()


def _page(raw, descriptor, term):
    page = CatalogPage(
        descriptor.get("request_url"),
        descriptor.get("retrieved_at_utc"),
        raw,
        _json(raw),
    )
    validate_search_page(page, term)
    _require(
        descriptor.get("query") == page.data.get("query"),
        "Registry publication source receipt mismatch",
    )
    return page

def _candidate(searches, accession):
    matches = [
        row for row in merge_search_results(searches)
        if _accession(row.get("accession")) == accession
    ]
    _require(
        len(matches) == 1,
        "Dataset has no unique explicit registry source",
    )
    return matches[0]


def _check_papers(papers, candidate, root, as_of, expected_root):
    references = extract_publication_links(candidate)["references"]
    _require(
        isinstance(papers, list) and papers,
        "Published metadata requires papers",
    )
    identities = set()
    for paper in papers:
        _require(isinstance(paper, dict), "Invalid paper metadata")
        identity = paper.get("doi"), paper.get("pmid")
        _require(
            all(value is None or isinstance(value, str) for value in identity),
            "Invalid paper identifier",
        )
        _require(
            any(identity) and identity not in identities,
            "Missing or duplicated paper identity",
        )
        identities.add(identity)
        evidence = paper.get("verification_evidence")
        _require(
            isinstance(evidence, list) and evidence,
            "Paper lacks captured verification evidence",
        )
        confirmed_dois = set()
        confirmed_pmids = set()
        allowed_fields = set()
        qualified = []
        for proof in evidence:
            _require(
                isinstance(proof, dict)
                and proof.get("root") == expected_root,
                "Invalid publication evidence root",
            )
            matching = [
                reference for reference in references
                if (reference.get("doi"), reference.get("pmid"))
                == (
                    proof.get("requested_doi"),
                    proof.get("requested_pmid"),
                )
            ]
            _require(
                matching,
                "Paper identifier is not explicitly linked by the source dataset",
            )
            for reference in matching:
                allowed_fields.update(reference["source_fields"])
            result = replay_publication_lookup(
                root / expected_root,
                proof.get("receipt_path"),
                receipt_sha256=proof.get("receipt_sha256"),
                doi=proof.get("requested_doi"),
                pmid=proof.get("requested_pmid"),
                as_of=as_of,
            )
            _require(
                result.get("status") == "published"
                and isinstance(result.get("paper"), dict),
                "Paper is not confirmed published by its captured evidence",
            )
            _require(
                proof.get("provider") == result.get("provider"),
                "Paper provider disagrees with captured evidence",
            )
            verified = result["paper"]
            qualified.append(verified)
            if verified.get("doi"):
                confirmed_dois.add(verified["doi"])
            if verified.get("pmid"):
                confirmed_pmids.add(verified["pmid"])
        _require(
            confirmed_dois
            == ({paper["doi"]} if paper.get("doi") else set())
            and confirmed_pmids
            == ({paper["pmid"]} if paper.get("pmid") else set()),
            "Paper identifiers disagree with verification evidence",
        )
        fields = paper.get("association_source_fields")
        _require(
            isinstance(fields, list)
            and fields
            and all(isinstance(field, str) for field in fields)
            and set(fields) <= allowed_fields,
            "Paper association source fields are inconsistent",
        )
        metadata_fields = (
            "title", "authors", "venue", "year", "publication_date",
            "publication_dates", "work_type", "url",
        )
        _require(
            any(
                all(
                    paper.get(key) == item.get(key)
                    for key in metadata_fields
                )
                for item in qualified
            ),
            "Paper metadata differs from its captured response",
        )


def _finish(record, bodies):
    _require(
        sum(len(raw) for raw in bodies.values()) <= MAX_PACKAGE_BYTES,
        "Dataset publication evidence package exceeds 16 MiB",
    )
    record["evidence_sha256"] = {
        name: _sha(raw) for name, raw in bodies.items()
    }
    return {"record": record, "source_bodies": bodies}


def _from_discovery(root, path, accession, manifest, manifest_raw):
    _require(
        path.name == "candidates.json",
        "Use discovery candidates.json for dataset metadata",
    )
    as_of = _discovery_policy(manifest, accession)
    candidates_raw = _verified(root, path.name, manifest)
    candidates = _json(candidates_raw)
    _require(
        isinstance(candidates, list),
        "Invalid discovery candidate table",
    )
    matches = [
        row for row in candidates
        if isinstance(row, dict)
        and _accession(row.get("accession")) == accession
    ]
    _require(
        len(matches) == 1
        and matches[0].get("publication_status") == "published",
        "Dataset does not have one confirmed published candidate",
    )
    selected = matches[0]
    bodies = {SOURCE_MANIFEST: manifest_raw}
    pages = []
    searches = {}
    _require(
        isinstance(manifest.get("searches"), list)
        and manifest["searches"],
        "Discovery has no registry source receipts",
    )
    for search in manifest["searches"]:
        term = search.get("search_term")
        _require(
            isinstance(term, str) and term,
            "Invalid source search term",
        )
        for descriptor in search.get("pages", []):
            name = descriptor.get("path")
            raw = _verified(root, name, manifest)
            _require(
                _sha(raw) == descriptor.get("sha256"),
                "Registry page checksum mismatch",
            )
            parsed = _page(raw, descriptor, term)
            columns = parsed.data["result_set"]["datasets_title_list"]
            index = columns.index("dataset identifier")
            if not any(
                _accession(row[index]) == accession
                for row in parsed.data["datasets"]
            ):
                continue
            target = REGISTRY_ROOT + name
            bodies[target] = raw
            pages.append({
                **deepcopy(descriptor),
                "path": target,
                "source_path": name,
                "search_term": term,
            })
            searches.setdefault(term, []).append(parsed)

    source_candidate = _candidate(searches, accession)
    _check_papers(
        selected.get("papers"), source_candidate,
        root, as_of, "publication_evidence",
    )
    papers = deepcopy(selected["papers"])
    for paper in papers:
        for proof in paper["verification_evidence"]:
            receipt_name = proof["receipt_path"]
            raw = _verified(
                root, "publication_evidence/" + receipt_name, manifest
            )
            _require(
                _sha(raw) == proof["receipt_sha256"],
                "Publication receipt checksum mismatch",
            )
            bodies[HTTP_ROOT + "/" + receipt_name] = raw
            receipt = _json(raw)
            for attempt in receipt["attempts"]:
                name = attempt["body_path"]
                body = _verified(
                    root, "publication_evidence/" + name, manifest
                )
                _require(
                    _sha(body) == attempt["body_sha256"],
                    "Publication body checksum mismatch",
                )
                bodies[HTTP_ROOT + "/" + name] = body
            proof["root"] = HTTP_ROOT

    record = {
        "schema_version": 1,
        "kind": "xlms_dataset_publications",
        "dataset_accession": accession,
        "publication_status": "published",
        "papers": papers,
        "source_candidate": source_candidate,
        "source_discovery": {
            "manifest_path": SOURCE_MANIFEST,
            "manifest_sha256": _sha(manifest_raw),
            "candidates_sha256": _sha(candidates_raw),
            "run_id": manifest.get("run_id"),
            "namespace_id": manifest.get("namespace_id"),
            "publication_as_of": as_of,
            "registry_pages": pages,
        },
    }
    return _finish(record, bodies)


def _from_snapshot(root, path, accession, manifest):
    raw = _verified(root, path.name, manifest)
    record = _json(raw)
    entry = manifest.get("dataset_metadata")
    _require(
        isinstance(entry, dict)
        and entry.get("path") == path.name
        and entry.get("sha256") == _sha(raw)
        and entry.get("dataset_accession") == accession,
        "Export metadata descriptor mismatch",
    )
    _require(
        record.get("schema_version") == 1
        and record.get("kind") == "xlms_dataset_publications"
        and record.get("dataset_accession") == accession
        and record.get("publication_status") == "published",
        "Dataset metadata identity or schema mismatch",
    )
    pinned = record.get("evidence_sha256")
    _require(
        isinstance(pinned, dict) and SOURCE_MANIFEST in pinned,
        "Dataset metadata lacks source evidence checksums",
    )
    bodies = {}
    for name, digest in pinned.items():
        body = _verified(root, name, manifest)
        _require(
            _sha(body) == digest,
            "Packaged metadata evidence checksum mismatch",
        )
        bodies[name] = body

    source = record.get("source_discovery")
    _require(
        isinstance(source, dict)
        and source.get("manifest_path") == SOURCE_MANIFEST
        and source.get("manifest_sha256")
        == _sha(bodies[SOURCE_MANIFEST]),
        "Discovery metadata source manifest mismatch",
    )
    original = _json(bodies[SOURCE_MANIFEST])
    as_of = _discovery_policy(original, accession)
    _require(
        source.get("publication_as_of") == as_of
        and source.get("run_id") == original.get("run_id")
        and source.get("namespace_id") == original.get("namespace_id"),
        "Discovery publication context does not match its source receipt",
    )
    for name, body in bodies.items():
        if name == SOURCE_MANIFEST:
            continue
        if name.startswith(HTTP_ROOT + "/"):
            original_name = (
                "publication_evidence/" + name[len(HTTP_ROOT) + 1:]
            )
        elif name.startswith(REGISTRY_ROOT):
            original_name = name[len(REGISTRY_ROOT):]
        else:
            raise ValueError("Unsupported publication evidence path")
        _require(
            original["output_sha256"].get(original_name) == _sha(body),
            "Packaged evidence differs from the original discovery manifest",
        )

    searches = {}
    pages = source.get("registry_pages")
    _require(
        isinstance(pages, list) and pages,
        "Missing registry association evidence",
    )
    for descriptor in pages:
        name = descriptor.get("path")
        _require(
            name == REGISTRY_ROOT + descriptor.get("source_path", "")
            and name in bodies
            and _sha(bodies[name]) == descriptor.get("sha256"),
            "Packaged registry page descriptor mismatch",
        )
        term = descriptor.get("search_term")
        _require(
            isinstance(term, str) and term,
            "Invalid source search term",
        )
        searches.setdefault(term, []).append(
            _page(bodies[name], descriptor, term)
        )
    candidate = _candidate(searches, accession)
    _require(
        candidate == record.get("source_candidate"),
        "Dataset association source was changed",
    )
    _check_papers(
        record.get("papers"), candidate, root, as_of, HTTP_ROOT
    )
    return _finish(deepcopy(record), bodies)


def load_dataset_metadata(path, *, dataset_accession):
    """Verify a discovery candidate or packaged metadata offline."""
    accession = _accession(dataset_accession)
    path = Path(path).expanduser().resolve()
    root = path.parent
    manifest, manifest_raw = _manifest(root)
    document = _json(_read(root, path.name))
    if isinstance(document, list):
        return _from_discovery(
            root, path, accession, manifest, manifest_raw
        )
    _require(
        isinstance(document, dict)
        and document.get("kind") == "xlms_dataset_publications",
        "Unsupported dataset publication metadata input",
    )
    return _from_snapshot(root, path, accession, manifest)


def bind_dataset_metadata(bundle, metadata_path, *, dataset_accession):
    """Bind metadata without changing source records or confidence."""
    if metadata_path is None:
        return bundle

    accession = _accession(dataset_accession)
    provenance = bundle.get("provenance", {})
    for owner in (bundle, provenance):
        for key in ("dataset_accession", "project_accession"):
            if owner.get(key) is not None:
                _require(
                    _accession(owner[key]) == accession,
                    "Bundle dataset accession conflicts with metadata",
                )
    for record in bundle.get("records", []):
        identifier = record.get("source_record_id")
        if isinstance(identifier, str) and ":" in identifier:
            prefix = identifier.split(":", 1)[0]
            if re.fullmatch(r"PXD[0-9]+", prefix, re.IGNORECASE):
                _require(
                    _accession(prefix) == accession,
                    "Source record dataset conflicts with metadata",
                )

    path = Path(metadata_path).expanduser().resolve()
    load_dataset_metadata(path, dataset_accession=accession)
    reference = {
        "path": str(path),
        "dataset_accession": accession,
        "file_sha256": _sha(path.read_bytes()),
        "manifest_sha256": _sha(
            (path.parent / "manifest.json").read_bytes()
        ),
    }
    copied = deepcopy(bundle)
    _require(
        "_dataset_publications" not in copied,
        "Dataset metadata is already bound",
    )
    copied["_dataset_publications"] = reference
    return copied


def prepare_publication_export(bundle):
    """Prepare snapshots and CSV fields before output creation."""
    reference = bundle.get("_dataset_publications")
    if reference is None:
        return None
    _require(
        isinstance(reference, dict),
        "Invalid bound dataset metadata reference",
    )
    path = Path(reference["path"]).resolve()
    _require(
        _sha(_read(path.parent, path.name))
        == reference.get("file_sha256")
        and _sha(
            _read(path.parent, "manifest.json", 8 * 1024 * 1024)
        ) == reference.get("manifest_sha256"),
        "Bound metadata changed after validation",
    )
    packet = load_dataset_metadata(
        path, dataset_accession=reference["dataset_accession"]
    )
    record = packet["record"]
    papers = record["papers"]
    fields = {
        "xlms_dataset_accession": record["dataset_accession"],
        "xlms_publication_status": record["publication_status"],
        "xlms_paper_dois_json": json.dumps(
            list(dict.fromkeys(
                paper["doi"] for paper in papers if paper.get("doi")
            )),
            ensure_ascii=False,
        ),
        "xlms_paper_pmids_json": json.dumps(
            list(dict.fromkeys(
                paper["pmid"] for paper in papers if paper.get("pmid")
            )),
            ensure_ascii=False,
        ),
        "xlms_papers_json": json.dumps(
            papers, ensure_ascii=False, allow_nan=False
        ),
    }
    metadata_raw = _encode(record)
    bodies = dict(packet["source_bodies"])
    bodies["dataset_metadata.json"] = metadata_raw
    _require(
        sum(len(raw) for raw in bodies.values()) <= MAX_PACKAGE_BYTES,
        "Dataset publication evidence package exceeds 16 MiB",
    )
    return {
        "record": record,
        "row_fields": fields,
        "source_bodies": bodies,
        "manifest_entry": {
            "path": "dataset_metadata.json",
            "sha256": _sha(metadata_raw),
            "dataset_accession": record["dataset_accession"],
        },
    }
