"""Retrieve and save public repository metadata and file inventories."""

import csv
import hashlib
import io
import json
import math
import re
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from email.message import Message
from pathlib import Path
from typing import List, Mapping, Union
from urllib.parse import urlencode, urlsplit, unquote, unquote_to_bytes, parse_qs, urljoin, quote

import requests

from . import _write_json


PRIDE_API_URL = "https://www.ebi.ac.uk/pride/ws/archive/v3"
MASSIVE_API_URL = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1"
MASSIVE_FILES_URL = "https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp"
REGISTRY_DATASETS_URL = "https://proteomecentral.proteomexchange.org/api/proxi/v0.1/datasets"
MASSIVE_DOWNLOAD_URL = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile"
JPOST_API_URL = "https://repository.jpostdb.org/proxi/datasets"
IPROX_API_URL = "https://www.iprox.cn/proxi/datasets"
PANORAMA_ORIGIN = "https://panoramaweb.org"
PANORAMA_QUERY_URL = PANORAMA_ORIGIN + "/Panorama%20Public/query-selectRows.api"


@dataclass(frozen=True)
class RepositoryResponse:
    request_url: str
    retrieved_at_utc: str
    raw_body: bytes
    headers: Mapping[str, str]
    data: Union[dict, list]


# Existing PRIDE callers use the same response fields.
PrideResponse = RepositoryResponse


def _json_response_data(raw_body, headers):
    """Decode JSON strictly, using a declared legacy charset only if necessary."""
    try:
        return json.loads(raw_body)
    except UnicodeDecodeError:
        message = Message()
        message["Content-Type"] = requests.structures.CaseInsensitiveDict(headers).get(
            "Content-Type", ""
        )
        charset = message.get_content_charset()
        if not charset:
            raise
        try:
            text = raw_body.decode(charset, errors="strict")
        except LookupError as error:
            raise ValueError(f"Unsupported response charset: {charset}") from error
        return json.loads(text)


def _panorama_url(base, value, directory=False):
    if not isinstance(value, str) or not value:
        raise ValueError("Panorama resource has no valid URL")
    uri = urlsplit(urljoin(base, value))
    if (
        uri.scheme != "https" or uri.hostname != "panoramaweb.org"
        or uri.port not in (None, 443)
        or uri.username is not None or uri.password is not None
        or uri.query or uri.fragment
    ):
        raise ValueError("Panorama resource URL must belong to the public server")
    path = unquote(uri.path, errors="strict")
    if any(part in ("", ".", "..") for part in path.strip("/").split("/")):
        raise ValueError("Panorama resource has an invalid path")
    if directory:
        path = path.rstrip("/") + "/"
    elif path.endswith("/"):
        raise ValueError("Panorama file URL points to a directory")
    return PANORAMA_ORIGIN + quote(path, safe="/@")


def _panorama_project(data, accession):
    if not isinstance(data, dict) or not isinstance(data.get("rows"), list):
        raise ValueError("Panorama project query has no rows array")
    rows = data["rows"]
    if len(rows) != 1 or type(data.get("rowCount")) is not int or data["rowCount"] != 1:
        raise ValueError(f"Panorama query must return one current project for {accession}")
    row = rows[0]
    if not isinstance(row, dict) or row.get("pxid") != accession:
        raise ValueError("Panorama project accession mismatch")
    native_id = row.get("Id")
    if not isinstance(native_id, int) or isinstance(native_id, bool) or native_id < 1:
        raise ValueError("Panorama project has no valid native identifier")
    project_url = _panorama_url(PANORAMA_ORIGIN, row.get("_labkeyurl_Title"))
    path = unquote(urlsplit(project_url).path, errors="strict")
    suffix = "/project-begin.view"
    if not path.startswith("/Panorama Public/") or not path.endswith(suffix):
        raise ValueError("Panorama project URL does not identify a public project folder")
    root = _panorama_url(
        PANORAMA_ORIGIN, "/_webdav" + path[:-len(suffix)] + "/", directory=True
    )
    return {"native_id": str(native_id), "project_url": project_url, "root_url": root}


def _panorama_get(url, params):
    """Retry only when the server supplies a short numeric Retry-After."""
    for attempt in range(3):
        response = requests.get(url, params=params, timeout=(15, 120))
        if response.status_code == 429 and attempt < 2:
            delay = response.headers.get("Retry-After", "")
            if delay.isdigit() and int(delay) <= 60:
                time.sleep(int(delay))
                continue
        response.raise_for_status()
        return response


def fetch_panorama_project(accession: str) -> RepositoryResponse:
    """Resolve one current Panorama Public project by its PXD accession."""
    if not isinstance(accession, str) or not re.fullmatch(r"PXD\d{6,}", accession):
        raise ValueError("A valid PXD accession is required")
    response = _panorama_get(
        PANORAMA_QUERY_URL,
        params={
            "schemaName": "panoramapublic", "queryName": "ExperimentAnnotations",
            "containerFilter": "CurrentAndSubfolders",
            "query.pxid~eq": accession, "query.Version~contains": "Current",
            "maxRows": 2,
        },
    )
    retrieved_at = datetime.now(timezone.utc).isoformat()
    data = json.loads(response.content)
    _panorama_project(data, accession)
    return RepositoryResponse(
        response.url, retrieved_at, response.content, response.headers.copy(), data
    )


def _panorama_children(directory, data, context, accession):
    if (
        not isinstance(data, dict) or not isinstance(data.get("files"), list)
        or data.get("success") is False or data.get("exception")
    ):
        raise ValueError("Panorama directory response has no files array")
    in_file_store = "/@files/" in unquote(urlsplit(directory).path, errors="strict")
    folders, files = [], []
    for entry in data["files"]:
        if not isinstance(entry, dict) or not isinstance(entry.get("collection"), bool):
            raise ValueError("Panorama directory contains an invalid resource")
        collection = entry["collection"]
        url = _panorama_url(directory, entry.get("href"), directory=collection)
        if not url.startswith(context["root_url"]) or not url.startswith(directory):
            raise ValueError("Panorama resource lies outside its directory or project")
        if url == directory and collection:
            continue
        if not in_file_store:
            # Container roots expose virtual resources; submitted files live
            # in @files. Ordinary child containers may have their own stores.
            name = unquote(url.rstrip("/").rsplit("/", 1)[-1], errors="strict")
            if not collection or (name.startswith("@") and name != "@files"):
                continue
        if entry.get("canRead") is False:
            raise PermissionError(f"Panorama resource is not publicly readable: {url}")
        if collection:
            folders.append(url)
            continue
        size = entry.get("contentlength", entry.get("size"))
        if size is not None and (
            not isinstance(size, int) or isinstance(size, bool) or size < 0
        ):
            raise ValueError("Panorama file has an invalid byte size")
        relative_path = unquote(url[len(context["root_url"]):], errors="strict")
        files.append({
            "accession": url,
            "projectAccessions": [accession],
            "repositoryAccessions": [context["native_id"]],
            "relativePath": relative_path,
            "fileName": relative_path.rsplit("/", 1)[-1],
            "fileCategory": {"value": ""},
            "fileSizeBytes": size,
            "checksum": None,
            "publicFileLocations": [{"name": "HTTPS", "value": url}],
            "repositoryMetadata": entry,
        })
    return folders, files


def fetch_panorama_file_pages(accession, project):
    """Read public file stores and child containers, retaining each listing."""
    if json.loads(project.raw_body) != project.data:
        raise ValueError("Parsed project metadata differs from its raw response")
    context = _panorama_project(project.data, accession)
    pending = [context["root_url"]]
    seen, pages = set(), []
    while pending:
        directory = pending.pop()
        if directory in seen:
            continue
        response = _panorama_get(directory, {"method": "json"})
        retrieved_at = datetime.now(timezone.utc).isoformat()
        data = json.loads(response.content)
        folders, _ = _panorama_children(directory, data, context, accession)
        pages.append(RepositoryResponse(
            response.url, retrieved_at, response.content, response.headers.copy(), data
        ))
        seen.add(directory)
        pending.extend(folders)
    _panorama_files(accession, project, pages)
    return pages


def _panorama_files(accession, project, pages):
    """Require a saved response for every directory reached from the project."""
    context = _panorama_project(project.data, accession)
    by_directory = {}
    for page in pages:
        if json.loads(page.raw_body) != page.data:
            raise ValueError("Parsed file metadata differs from its raw response")
        directory = _panorama_url(
            PANORAMA_ORIGIN, page.request_url.split("?", 1)[0], directory=True
        )
        if not directory.startswith(context["root_url"]) or directory in by_directory:
            raise ValueError("Panorama directory response is repeated or outside the project")
        by_directory[directory] = page

    pending = [context["root_url"]]
    seen, file_ids, files = set(), set(), []
    has_file_store = False
    while pending:
        directory = pending.pop()
        if directory in seen:
            continue
        if directory not in by_directory:
            raise ValueError(f"Incomplete Panorama inventory: missing {directory}")
        seen.add(directory)
        has_file_store |= "/@files/" in unquote(urlsplit(directory).path, errors="strict")
        folders, records = _panorama_children(
            directory, by_directory[directory].data, context, accession
        )
        pending.extend(folders)
        for record in records:
            if record["accession"] in file_ids:
                raise ValueError("Panorama inventory contains a repeated file URL")
            file_ids.add(record["accession"])
            files.append(record)
    if not has_file_store or seen != set(by_directory):
        raise ValueError("Panorama inventory has no file store or has unrelated responses")
    return files, None


def _proxi_identifiers(data, repository, native_name, native_cv, native_pattern):
    if not isinstance(data, dict):
        raise ValueError(f"{repository} project response must be a JSON object")
    entries = data.get("accession")
    if not isinstance(entries, list):
        raise ValueError(f"{repository} project has no identifier list")

    identifiers = {}
    for name, cv_id, pattern in (
        ("pxd", "MS:1001919", r"PXD\d{6,}"),
        (native_name, native_cv, native_pattern),
    ):
        values = [
            entry.get("value") for entry in entries
            if isinstance(entry, dict) and entry.get("accession") == cv_id
        ]
        if not values or any(
            not isinstance(value, str) or not re.fullmatch(pattern, value)
            for value in values
        ):
            raise ValueError(f"{repository} project has no valid {name.upper()} identifier")
        if len(set(values)) != 1:
            raise ValueError(f"{repository} project has conflicting {name.upper()} identifiers")
        identifiers[name] = values[0]
    return identifiers


def _massive_identifiers(data):
    if isinstance(data, dict) and "accession" not in data:
        summary = data.get("datasetSummary")
        if not isinstance(summary, dict) or summary.get("hostingRepository") != "MassIVE":
            raise ValueError("Registry project is not hosted by MassIVE")
        data = {"accession": data.get("identifiers")}
    return _proxi_identifiers(data, "MassIVE", "msv", "MS:1002487", r"MSV\d{9}")


def _massive_task_id(data):
    """Return a task only for the registry's explicitly linked MassIVE page."""
    if "accession" in data:
        return None
    links = data.get("fullDatasetLinks")
    if not isinstance(links, list):
        raise ValueError("Registry project has no fullDatasetLinks array")
    tasks = set()
    for link in links:
        if not isinstance(link, dict) or link.get("accession") != "MS:1002488":
            continue
        value = link.get("value")
        if not isinstance(value, str):
            raise ValueError("Registry MassIVE project link is invalid")
        uri = urlsplit(value)
        values = parse_qs(uri.query).get("task", [])
        if (
            uri.scheme not in ("http", "https") or uri.hostname != "massive.ucsd.edu"
            or uri.username is not None or uri.password is not None
            or uri.port not in (None, 80, 443) or uri.fragment
            or uri.path != "/ProteoSAFe/dataset.jsp"
            or len(values) != 1 or not re.fullmatch(r"[0-9a-f]{32}", values[0])
        ):
            raise ValueError("Registry MassIVE project link has no valid task")
        tasks.add(values[0])
    if len(tasks) != 1:
        raise ValueError("Registry project must identify one unique MassIVE task")
    return tasks.pop()


def _jpost_identifiers(data):
    return _proxi_identifiers(data, "jPOST", "jpst", "MS:1002632", r"JPST\d{6,}")


def _iprox_identifiers(data):
    if not isinstance(data, dict):
        raise ValueError("iProX project response must be a JSON object")
    entries = data.get("accession")
    if isinstance(entries, dict):
        entries = [entries]
    if not isinstance(entries, list):
        raise ValueError("iProX project has no accession information")
    values = [
        entry.get("value") for entry in entries
        if isinstance(entry, dict) and entry.get("accession") == "MS:1001919"
    ]
    if not values or any(
        not isinstance(value, str) or not re.fullmatch(r"PXD\d{6,}", value)
        for value in values
    ) or len(set(values)) != 1:
        raise ValueError("iProX project has no unique valid PXD identifier")

    links = data.get("datasetLink")
    if isinstance(links, dict):
        links = [links]
    if not isinstance(links, list):
        raise ValueError("iProX project has no datasetLink information")
    native_ids = []
    for link in links:
        if not isinstance(link, dict) or link.get("accession") != "MS:1002837":
            continue
        if not isinstance(link.get("value"), str):
            raise ValueError("iProX dataset URI is invalid")
        uri = urlsplit(link["value"])
        identifiers = parse_qs(uri.query).get("id", [])
        if (
            uri.scheme not in ("http", "https")
            or uri.hostname not in ("iprox.org", "www.iprox.org", "iprox.cn", "www.iprox.cn")
            or uri.username is not None or uri.password is not None
            or not identifiers
            or any(not re.fullmatch(r"IPX\d{7,10}", value) for value in identifiers)
        ):
            raise ValueError("iProX dataset URI has no valid IPX identifier")
        native_ids.extend(identifiers)
    if len(set(native_ids)) != 1:
        raise ValueError("iProX project has no unique IPX identifier")
    return {"pxd": values[0], "ipx": native_ids[0]}


def _select_proxi_project(data, accession, repository):
    """Find one matching project while keeping the original response intact."""
    if isinstance(data, dict):
        projects = [data]
    elif isinstance(data, list):
        projects = data
    else:
        raise ValueError(f"{repository} response must be a project object or an array")
    matches = []
    for project in projects:
        if not isinstance(project, dict):
            raise ValueError(f"{repository} response contains an invalid project object")
        identifiers = project.get("accession")
        if isinstance(identifiers, dict):
            identifiers = [identifiers]
        if not isinstance(identifiers, list):
            raise ValueError(f"{repository} project has no identifier list")
        if any(
            isinstance(entry, dict)
            and entry.get("accession") == "MS:1001919"
            and entry.get("value") == accession
            for entry in identifiers
        ):
            matches.append(project)
    if len(matches) != 1:
        raise ValueError(
            f"{repository} response must contain one unique project matching {accession}; "
            f"found {len(matches)}"
        )
    return matches[0]


def _select_jpost_project(data, accession):
    return _select_proxi_project(data, accession, "jPOST")


def _jpost_files(accession, data):
    """Normalize reported file links without inventing sizes or a total."""
    data = _select_jpost_project(data, accession)
    identifiers = _jpost_identifiers(data)
    if identifiers["pxd"] != accession:
        raise ValueError("jPOST project accession mismatch")
    entries = data.get("datasetFiles")
    if not isinstance(entries, list):
        raise ValueError("jPOST response has no datasetFiles array")
    return _reported_links(
        "jPOST", accession, identifiers["jpst"], entries,
        "ftp.jpostdb.org", ("ftp", "ftps"),
    )


def _iprox_files(accession, data):
    data = _select_proxi_project(data, accession, "iProX")
    identifiers = _iprox_identifiers(data)
    if identifiers["pxd"] != accession:
        raise ValueError("iProX project accession mismatch")
    entries = data.get("dataFiles")
    if isinstance(entries, dict):
        entries = [entries]
    if not isinstance(entries, list):
        raise ValueError("iProX response has no dataFiles array or object")
    return _reported_links(
        "iProX", accession, identifiers["ipx"], entries,
        "download.iprox.org", ("http", "https"),
    )


def _reported_links(repository, accession, native_id, entries, host, schemes):
    """Keep reported URLs and paths without inferring missing file metadata."""
    prefix = f"/{native_id}/"
    files = []
    seen = {}
    for entry_index, entry in enumerate(entries):
        if not isinstance(entry, dict) or not isinstance(entry.get("value"), str):
            raise ValueError(f"{repository} file link has no valid URI")
        uri = entry["value"]
        parsed = urlsplit(uri)
        path_bytes = unquote_to_bytes(parsed.path)
        unresolved_path = False
        try:
            path = path_bytes.decode("utf-8", errors="strict")
        except UnicodeDecodeError:
            # This is a catalog entry, not a filename to repair or download.
            path = parsed.path
            unresolved_path = True
        if (
            parsed.scheme not in schemes
            or parsed.hostname != host
            or parsed.username is not None or parsed.password is not None
            or parsed.query or parsed.fragment
            or not path.startswith(prefix)
            or not path_bytes.startswith(prefix.encode("ascii"))
        ):
            raise ValueError(f"{repository} file link does not belong to {native_id}: {uri}")
        relative_path = path[len(prefix):]
        relative_bytes = path_bytes[len(prefix.encode("ascii")):]
        if any(part in (b"", b".", b"..") for part in relative_bytes.split(b"/")):
            raise ValueError(f"{repository} file link has an invalid relative path")
        if uri in seen:
            file_index, first_entry_index = seen[uri]
            previous = files[file_index]
            if repository != "jPOST" or previous["repositoryMetadata"] != entry:
                raise ValueError(f"Repeated {repository} file link has conflicting metadata: {uri}")
            indices = previous.setdefault("reportedLinkIndices", [first_entry_index])
            indices.append(entry_index)
            previous["reportedLinkCount"] = len(indices)
            continue
        seen[uri] = (len(files), entry_index)
        files.append({
            "accession": uri,
            "projectAccessions": [accession],
            "repositoryAccessions": [native_id],
            "relativePath": relative_path,
            "fileName": relative_path.rsplit("/", 1)[-1],
            "fileCategory": {"value": entry.get("name", "")},
            "fileSizeBytes": None,
            "checksum": None,
            "publicFileLocations": [{
                "name": "FTP" if repository == "jPOST" else parsed.scheme.upper(),
                "value": uri,
            }],
            "repositoryMetadata": entry,
        })
        if unresolved_path:
            files[-1]["pathDecodingStatus"] = "unresolved"
    return files, None


def _parse_massive_file_page(raw_body):
    """Read the embedded JSON without evaluating the page's JavaScript."""
    text = raw_body.decode("utf-8")
    marker = re.search(r"\bvar\s+dataset_files\s*=\s*", text)
    if marker is None:
        raise ValueError("MassIVE page has no embedded file listing")
    data, _ = json.JSONDecoder().raw_decode(text, marker.end())
    if not isinstance(data, dict):
        raise ValueError("MassIVE file listing must be a JSON object")
    return data


def _massive_file_records(msv, data):
    records = data.get("row_data")
    total = data.get("total_rows")
    if not isinstance(records, list):
        raise ValueError("MassIVE file listing has no row_data array")
    if not isinstance(total, int) or isinstance(total, bool) or total < 0:
        raise ValueError("MassIVE file listing has no valid total_rows count")
    if len(records) != total:
        raise ValueError(
            f"MassIVE file count mismatch: received {len(records)}, expected {total}"
        )

    seen = set()
    prefix = f"f.{msv}/"
    for record in records:
        if not isinstance(record, dict) or record.get("type") != "File":
            raise ValueError("MassIVE listing contains an invalid file record")
        descriptor = record.get("file_descriptor")
        if not isinstance(descriptor, str) or not descriptor.startswith(prefix):
            raise ValueError(f"MassIVE file does not belong to {msv}")
        relative_path = descriptor[len(prefix):]
        if any(part in ("", ".", "..") for part in relative_path.split("/")):
            raise ValueError("MassIVE file has an invalid relative path")
        if descriptor in seen:
            raise ValueError(f"Repeated MassIVE file: {descriptor}")
        seen.add(descriptor)
        size = record.get("size")
        if not isinstance(size, int) or isinstance(size, bool) or size < 0:
            raise ValueError(f"MassIVE file has no valid byte size: {descriptor}")
    return records, total


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
        "source_kind": "library_export" if response.request_url.startswith("pridepy://") else "http_response",
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


def _repository_files(repository, accession, project, file_pages):
    """Validate raw sources and produce the common file records."""
    if _json_response_data(project.raw_body, project.headers) != project.data:
        raise ValueError("Parsed project metadata differs from its raw response")
    if file_pages and isinstance(file_pages[0].data, dict) and file_pages[0].data.get("kind") == EXPORT_KIND:
        return _pridepy_files(repository, accession, project, file_pages)
    if repository == "PRIDE":
        if not isinstance(project.data, dict) or project.data.get("accession") != accession:
            raise ValueError("PRIDE project accession mismatch")
        return _validated_files(accession, file_pages)
    if repository == "jPOST":
        if file_pages:
            raise ValueError("jPOST file links must come from its project response")
        return _jpost_files(accession, project.data)
    if repository == "iProX":
        if file_pages:
            raise ValueError("iProX file links must come from its project response")
        return _iprox_files(accession, project.data)
    if repository == "PanoramaPublic":
        return _panorama_files(accession, project, file_pages)
    if repository != "MassIVE":
        raise ValueError(f"Unsupported repository: {repository}")

    identifiers = _massive_identifiers(project.data)
    if identifiers["pxd"] != accession:
        raise ValueError("MassIVE project accession mismatch")
    task = _massive_task_id(project.data)
    if task is not None:
        source = urlsplit(project.request_url)
        if (
            source.scheme != "https"
            or source.hostname != "proteomecentral.proteomexchange.org"
            or source.port not in (None, 443)
            or source.username is not None or source.password is not None
            or source.path != f"/api/proxi/v0.1/datasets/{accession}"
            or parse_qs(source.query).get("resultType") != ["extended"]
            or source.fragment
        ):
            raise ValueError("Fallback project metadata has an invalid registry source")
    if len(file_pages) != 1:
        raise ValueError("A complete MassIVE main-site file response is required")
    page = file_pages[0]
    if task is not None:
        source = urlsplit(page.request_url)
        if (
            source.scheme != "https" or source.hostname != "massive.ucsd.edu"
            or source.port not in (None, 443)
            or source.username is not None or source.password is not None
            or source.path != "/ProteoSAFe/dataset_files.jsp"
            or parse_qs(source.query).get("task") != [task] or source.fragment
        ):
            raise ValueError("MassIVE file listing does not match the registry task")
    if _parse_massive_file_page(page.raw_body) != page.data:
        raise ValueError("Parsed file metadata differs from its raw response")
    msv = identifiers["msv"]
    records, total = _massive_file_records(msv, page.data)
    prefix = f"f.{msv}/"
    files = []
    for record in records:
        descriptor = record["file_descriptor"]
        relative_path = descriptor[len(prefix):]
        files.append({
            **record,
            "accession": descriptor,
            "projectAccessions": [accession],
            "repositoryAccessions": [msv],
            "relativePath": relative_path,
            "fileName": relative_path.rsplit("/", 1)[-1],
            "fileCategory": {"value": record.get("collection", "")},
            "fileSizeBytes": record["size"],
            "publicFileLocations": [{
                "name": "HTTPS",
                "value": MASSIVE_DOWNLOAD_URL + "?" + urlencode({
                    "forceDownload": "true", "file": descriptor,
                }),
            }],
        })
    return files, total


def write_pride_outputs(accession, project, file_pages, output_dir):
    """Write a PRIDE snapshot using the shared repository output format."""
    return write_repository_outputs("PRIDE", accession, project, file_pages, output_dir)


def write_repository_outputs(repository, accession, project, file_pages, output_dir):
    """Write a complete repository project snapshot into a new directory.

    Parameters
    ----------
    repository : str
        Supported repository name: PRIDE, MassIVE, jPOST, iProX or PanoramaPublic.
    accession : str
        ProteomeXchange project accession.
    project : RepositoryResponse
        Retrieved project metadata.
    file_pages : list of RepositoryResponse
        Complete retrieved file-list responses.
    output_dir : str or pathlib.Path
        New directory for the snapshot.

    Returns
    -------
    output_dir : pathlib.Path
    """
    files, reported_total = _repository_files(repository, accession, project, file_pages)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "raw").mkdir()

    project_source = _save_response(
        output_dir, Path("raw") / "project.json", project
    )
    file_sources = [
        _save_response(
            output_dir,
            Path("raw") / (
                f"files_page_{number:04d}.html" if repository == "MassIVE" and not page.request_url.startswith("pridepy://")
                else f"files_page_{number:04d}.json"
            ),
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
        "source": repository,
        "project_accession": accession,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "screening_status": "pending",
        "file_count": len(files),
        "reported_file_count": reported_total,
        "project_response": project_source,
        "file_responses": file_sources,
    }
    is_library = bool(file_pages and isinstance(file_pages[0].data, dict) and file_pages[0].data.get("kind") == EXPORT_KIND)
    if is_library:
        export = file_pages[0].data
        manifest.update({"backend": "pridepy", "backend_version": export["backend_version"],
                         "repository_accession": export["repository_accession"],
                         "file_inventory_scope": "provider_returned_files",
                         "provider_record_count": len(export["records"]),
                         "duplicate_record_count": len(export["records"]) - len(files)})
    elif repository == "MassIVE":
        manifest["repository_accession"] = _massive_identifiers(project.data)["msv"]
        manifest["project_metadata_source"] = (
            "ProteomeCentral" if _massive_task_id(project.data) is not None
            else "MassIVE PROXI"
        )
    elif repository == "jPOST":
        manifest["repository_accession"] = _jpost_identifiers(
            _select_jpost_project(project.data, accession)
        )["jpst"]
        manifest["file_inventory_scope"] = "repository_reported_links"
        reported_links = _select_jpost_project(project.data, accession)["datasetFiles"]
        manifest["reported_link_count"] = len(reported_links)
        manifest["duplicate_link_count"] = len(reported_links) - len(files)
    elif repository == "iProX":
        manifest["repository_accession"] = _iprox_identifiers(
            _select_proxi_project(project.data, accession, repository)
        )["ipx"]
        manifest["file_inventory_scope"] = "repository_reported_links"
    elif repository == "PanoramaPublic":
        context = _panorama_project(project.data, accession)
        manifest["repository_accession"] = context["native_id"]
        manifest["repository_project_url"] = context["project_url"]
        manifest["file_inventory_scope"] = "recursive_webdav_listing"
    if repository in ("jPOST", "iProX") or is_library:
        manifest["unresolved_path_count"] = sum(
            record.get("pathDecodingStatus") == "unresolved" for record in files
        )
    _write_json(output_dir / "manifest.json", manifest)
    return output_dir


FILE_SCREENING_RULE_VERSION = "metadata-file-roles-v2"


def classify_inventory_file(record):
    """Suggest a file role from metadata, without inspecting its contents.

    The returned role and download priority are screening suggestions. They do
    not confirm crosslinks, search database use, FDR, or scientific eligibility.
    First downloads prioritize XL-hinted CSV/TSV/TXT tables and FASTA candidates.
    Other result candidates remain available for review and later format-specific
    processing. Reported URLs and filenames are never repaired by this function.
    """
    name = record.get("fileName")
    if not isinstance(name, str) or not name.strip():
        relative = record.get("relativePath", "")
        name = relative.rsplit("/", 1)[-1] if isinstance(relative, str) else ""
    name = name.casefold()
    archive = next((suffix for suffix in (
        ".tar.gz", ".tar.bz2", ".tar.xz", ".tgz", ".tbz2", ".txz",
        ".zip", ".tar", ".7z", ".rar",
    ) if name.endswith(suffix)), "")
    compression = archive or next((suffix for suffix in (".gz", ".bz2", ".xz")
                                   if name.endswith(suffix)), "")
    logical_name = name[:-len(compression)] if compression else name
    extension = "." + logical_name.rsplit(".", 1)[-1] if "." in logical_name else ""
    xl_hint = bool(re.search(
        r"(?:^|[^a-z0-9])(?:cross[-_ ]?links?(?:ed|ing)?|p[-_]?link[0-9]*|"
        r"xisearch|xifdr[0-9.]*|merox|xlinkx|kojak)(?:[^a-z0-9]|$)", logical_name
    ))
    evidence = []
    roles = set()
    spectra_kind = ""
    if extension in (".fasta", ".fa", ".faa", ".fas"):
        roles.add("search_fasta_candidate")
        evidence.append(f"FASTA filename extension: {extension}")
    elif extension in (".raw", ".wiff") or logical_name.endswith(".wiff.scan"):
        roles.add("mass_spectra")
        spectra_kind = "raw"
        evidence.append("Raw mass-spectrometry filename extension")
    elif extension in (".mgf", ".ms2", ".mzml", ".mzxml"):
        roles.add("mass_spectra")
        spectra_kind = "peak_list" if extension in (".mgf", ".ms2") else "converted"
        evidence.append(f"Spectrum filename extension: {extension}")
    elif extension in (".mzid", ".mzidentml", ".mztab"):
        roles.add("identification_result_candidate")
        evidence.append(f"Identification result filename extension: {extension}")
    elif xl_hint and extension in (".csv", ".tsv", ".txt", ".xls", ".xlsx"):
        roles.add("identification_result_candidate")
        evidence.append("Table filename contains a crosslink or XL-MS tool hint")
    elif extension in (".pdf", ".doc", ".docx", ".odt", ".rtf", ".ppt", ".pptx"):
        roles.add("supplementary_candidate")
        evidence.append(f"Supporting document filename extension: {extension}")
    elif (any(word in logical_name for word in ("supplementary", "supplemental", "suppl_", "supporting_"))
          or logical_name.startswith("readme")):
        roles.add("supplementary_candidate")
        evidence.append("Filename suggests supplementary material or a readme")

    category = record.get("fileCategory", "")
    if isinstance(category, dict):
        category = category.get("value", "")
    if not isinstance(category, str):
        category = ""
    category_key = re.sub(r"[^a-z0-9]+", " ", category.casefold()).strip()
    category_roles = {
        "raw": "mass_spectra", "peak": "mass_spectra", "ccms peak": "mass_spectra",
        "associated raw file uri": "mass_spectra", "peak list file uri": "mass_spectra",
        "result": "identification_result_candidate", "search": "identification_result_candidate",
        "search engine output file uri": "identification_result_candidate",
        "fasta": "search_fasta_candidate", "sequence": "search_fasta_candidate",
        "sequence database file uri": "search_fasta_candidate", "search database": "search_fasta_candidate",
        "supplementary": "supplementary_candidate", "supplemental": "supplementary_candidate",
    }
    if category:
        evidence.append(f"Repository category: {category}")
    category_role = category_roles.get(category_key)
    if category_role:
        roles.add(category_role)
    if not evidence:
        evidence.append("No specific role identified from the available metadata")

    if archive:
        role, status = "archive", "archive_contents_unknown"
        evidence.append(f"Archive suffix {archive}; members have not been inspected")
    elif len(roles) > 1:
        role, status = "uncertain", "conflicting_metadata"
        evidence.append("Filename and repository category suggest different roles")
    elif roles:
        role, status = next(iter(roles)), "candidate"
    else:
        role, status = "uncertain", "insufficient_metadata"
    if role == "mass_spectra" and not spectra_kind:
        spectra_kind = "repository_reported"

    urls = []
    locations = record.get("publicFileLocations") or []
    if isinstance(locations, dict):
        locations = [locations]
    for location in locations if isinstance(locations, list) else []:
        uri = location.get("value") if isinstance(location, dict) else None
        if not isinstance(uri, str):
            continue
        try:
            parsed = urlsplit(uri)
            valid = parsed.scheme in ("http", "https", "ftp", "ftps") and parsed.hostname
        except ValueError:
            valid = False
        if valid and uri not in urls:
            urls.append(uri)
    blockers = []
    if not urls:
        blockers.append("no_public_download_url")
    if record.get("pathDecodingStatus") == "unresolved":
        blockers.append("unresolved_path_encoding")
    size = record.get("fileSizeBytes")
    if isinstance(size, (int, float)) and not isinstance(size, bool) and size == 0:
        blockers.append("zero_byte_file")
    if blockers:
        priority = "blocked"
    elif role == "search_fasta_candidate":
        priority = "first"
    elif (role == "identification_result_candidate" and xl_hint
          and extension in (".csv", ".tsv", ".txt")):
        priority = "first"
    elif role == "identification_result_candidate":
        priority = "review"
        evidence.append("Result candidate retained; the initial download pass prioritizes XL-hinted text tables")
    elif role == "mass_spectra":
        priority = "defer"
    else:
        priority = "review"
    return {
        "file_role": role,
        "classification_status": status,
        "possible_roles": sorted(roles),
        "xl_specific_hint": xl_hint,
        "spectra_kind": spectra_kind,
        "download_priority": priority,
        "download_blockers": blockers,
        "evidence": evidence,
        "download_urls": urls,
    }


def probe_pilot_content(body, role):
    """Describe bounded file contents without accepting scientific identifications."""
    try:
        encoding = "utf-16" if body.startswith((b"\xff\xfe", b"\xfe\xff")) else "utf-8-sig"
        text = body.decode(encoding, errors="strict")
    except UnicodeDecodeError as error:
        return {"status": "text_encoding_unresolved", "error": str(error)}
    if re.match(r"\s*(?:<!doctype\s+html|<html\b)", text, flags=re.IGNORECASE):
        return {"status": "html_response"}
    if role == "fasta_candidate":
        if not text.lstrip().startswith(">"):
            return {"status": "not_fasta"}
        ids, lengths, seen = [], [], set()
        count = duplicates = empty = invalid = length = 0
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if count:
                    empty += length == 0
                count += 1
                identifier = line[1:].split()[0] if line[1:].split() else ""
                invalid += not bool(identifier)
                duplicates += identifier in seen
                seen.add(identifier)
                if len(ids) < 5:
                    ids.append(identifier[:200])
                if count > 1 and len(lengths) < 5:
                    lengths.append(length)
                length = 0
            else:
                sequence = "".join(line.split()).upper()
                invalid += bool(re.search(r"[^A-Z*.-]", sequence))
                length += len(sequence)
        empty += length == 0
        if len(lengths) < 5:
            lengths.append(length)
        return {
            "status": "fasta_candidate" if count and not empty and not invalid else "fasta_needs_review",
            "sequence_count": count, "empty_sequence_count": empty,
            "invalid_line_or_header_count": invalid, "duplicate_identifier_count": duplicates,
            "example_ids": ids, "example_sequence_lengths": lengths,
            "search_database_association": "unconfirmed",
        }

    best = None
    for delimiter in (",", "\t", ";"):
        try:
            reader = csv.reader(io.StringIO(text), delimiter=delimiter)
            for number, row in enumerate(reader):
                if number >= 30:
                    break
                labels = [re.sub(r"[^a-z0-9]", "", cell.casefold()) for cell in row]
                hints = sum(bool(re.search(r"protein|accession|peptide|residue|position|link|site|score|fdr|decoy|scan|spectrum", label))
                            for label in labels)
                if len(row) > 1 and hints >= 2:
                    score = (hints, len(row), -number)
                    if best is None or score > best[0]:
                        best = (score, delimiter, number, row)
        except csv.Error:
            continue
    if best is None:
        return {
            "status": "text_control" if role.endswith("_control") else "unstructured_text",
            "nonempty_line_count": sum(bool(line.strip()) for line in text.splitlines()),
            "preview": text[:500],
        }
    _, delimiter, header_index, header = best
    rows_after = ragged = 0
    examples = []
    try:
        for number, row in enumerate(csv.reader(io.StringIO(text), delimiter=delimiter)):
            if number <= header_index or not any(cell.strip() for cell in row):
                continue
            rows_after += 1
            ragged += len(row) != len(header)
            if len(examples) < 3:
                examples.append([cell[:300] for cell in row[:20]])
    except csv.Error as error:
        return {"status": "table_parse_error", "error": str(error)}
    status = "table_candidate" if rows_after else "header_only_table"
    if role.endswith("_control"):
        status = "table_control"
    return {
        "status": status, "delimiter": delimiter,
        "header_record_index": header_index, "columns": [cell[:300] for cell in header[:200]],
        "rows_after_header": rows_after, "ragged_row_count": ragged,
        "example_rows": examples, "column_mapping": "unconfirmed",
    }


def _pilot_transfer_url(item, reported_url):
    """Use the verified iProX HTTPS host while retaining the reported URI."""
    uri = urlsplit(reported_url)
    if (item["repository"] == "iProX" and uri.scheme == "https"
            and uri.hostname == "download.iprox.org" and uri.port in (None, 443)
            and uri.username is None and uri.password is None):
        return uri._replace(netloc="download.iprox.cn").geturl(), "Verified iProX HTTPS download-host alias"
    return reported_url, ""


PRIDEPY_VERSION = "0.0.16"
PRIDEPY_REPOSITORIES = ("PRIDE", "MassIVE", "jPOST", "iProX")
EXPORT_KIND = "pridepy_file_inventory"


def _pridepy_version():
    from importlib.metadata import PackageNotFoundError, version

    try:
        installed = version("pridepy")
    except PackageNotFoundError as error:
        raise RuntimeError('Install the dataset dependency: python -m pip install "pridepy==0.0.16"') from error
    if installed != PRIDEPY_VERSION:
        raise RuntimeError(f"This adapter requires pridepy=={PRIDEPY_VERSION}; found {installed}")
    return installed


def _pridepy_project(accession):
    from pridepy.project.project import Project

    return Project().get_by_accession(accession)


def _pridepy_provider(repository):
    from pridepy.download.pride import PrideProvider
    from pridepy.download.massive import MassiveProvider
    from pridepy.download.jpost import JpostProvider
    from pridepy.download.iprox import IproxProvider

    return {"PRIDE": PrideProvider, "MassIVE": MassiveProvider,
            "jPOST": JpostProvider, "iProX": IproxProvider}[repository]()


def _pridepy_download(url, output_dir):
    from pridepy.download.client import Client

    Client.download_files_by_url(
        urls=[url], output_folder=str(output_dir),
        skip_if_downloaded_already=False, parallel_files=1, checksum_check=False,
    )


def _library_response(data, repository, accession, operation):
    version = _pridepy_version()
    body = (json.dumps(data, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode("utf-8")
    return RepositoryResponse(
        f"pridepy://{version}/{repository}/{accession}/{operation}",
        datetime.now(timezone.utc).isoformat(), body,
        {"Content-Type": "application/json", "X-CLAUDIO-Source-Kind": "library_export"},
        json.loads(body),
    )


def _registry_native_id(data, repository, accession):
    if not isinstance(data, dict) or data.get("datasetSummary", {}).get("hostingRepository") != repository:
        raise ValueError("Registry hosting repository differs from the candidate")
    code, pattern = {"MassIVE": ("MS:1002487", r"MSV\d{9}"),
                     "jPOST": ("MS:1002632", r"JPST\d{6}"),
                     "iProX": ("MS:1002836", r"IPX\d{7,10}")}[repository]
    identifiers = _proxi_identifiers(
        {"accession": data.get("identifiers")}, repository, "native", code, pattern,
    )
    if identifiers["pxd"] != accession:
        raise ValueError("Registry project accession mismatch")
    return identifiers["native"]


def _fetch_registry_project(accession, repository):
    # The registry binds a PXD to the native accession required by pridepy.
    response = requests.get(
        f"{REGISTRY_DATASETS_URL}/{accession}",
        params={"resultType": "extended"}, timeout=(15, 120),
    )
    response.raise_for_status()
    data = _json_response_data(response.content, response.headers)
    _registry_native_id(data, repository, accession)
    return RepositoryResponse(response.url, datetime.now(timezone.utc).isoformat(),
                              response.content, response.headers.copy(), data)


def fetch_pride_project(accession):
    _pridepy_version()
    data = _pridepy_project(accession)
    if not isinstance(data, dict) or data.get("accession") != accession:
        raise ValueError("PRIDE project accession mismatch")
    return _library_response(data, "PRIDE", accession, "get_by_accession")


def _fetch_library_files(repository, accession, native):
    _pridepy_version()
    records = _pridepy_provider(repository).list_files(native)
    if not isinstance(records, list):
        raise ValueError(f"pridepy returned no valid file list for {native}")
    return [_library_response({
        "kind": EXPORT_KIND, "backend": "pridepy", "backend_version": PRIDEPY_VERSION,
        "repository": repository, "project_accession": accession,
        "repository_accession": native, "records": records,
    }, repository, native, "list_files")]


def fetch_pride_file_pages(accession, *, page_size=100):
    if type(page_size) is not int or page_size < 1:
        raise ValueError("page_size must be a positive integer")
    # pridepy owns listing/pagination; page_size remains accepted by old callers.
    return _fetch_library_files("PRIDE", accession, accession)


def fetch_massive_project(accession):
    return _fetch_registry_project(accession, "MassIVE")


def fetch_jpost_project(accession):
    return _fetch_registry_project(accession, "jPOST")


def fetch_iprox_project(accession):
    return _fetch_registry_project(accession, "iProX")


def fetch_repository_inventory(candidate, *, on_stage=None):
    repository, accession = candidate["repository"], candidate["accession"]
    if on_stage:
        on_stage("project_metadata")
    if repository == "PanoramaPublic":
        project = fetch_panorama_project(accession)
        if on_stage:
            on_stage("file_list")
        return project, fetch_panorama_file_pages(accession, project)
    if repository not in PRIDEPY_REPOSITORIES:
        raise ValueError(f"Unsupported repository: {repository}")
    _pridepy_version()
    if repository == "PRIDE":
        project, native = fetch_pride_project(accession), accession
    else:
        project = _fetch_registry_project(accession, repository)
        native = _registry_native_id(project.data, repository, accession)
    if on_stage:
        on_stage("file_list")
    return project, _fetch_library_files(repository, accession, native)


def _pridepy_files(repository, accession, project, pages):
    """Normalize library records; the saved export retains their original values."""
    if len(pages) != 1:
        raise ValueError("A single complete pridepy export is required")
    page = pages[0]
    data = page.data
    if (_json_response_data(page.raw_body, page.headers) != data
            or data.get("kind") != EXPORT_KIND or data.get("backend") != "pridepy"
            or data.get("backend_version") != PRIDEPY_VERSION
            or data.get("repository") != repository or data.get("project_accession") != accession
            or not isinstance(data.get("records"), list)):
        raise ValueError("Invalid or mismatched pridepy export")
    if repository == "PRIDE":
        if not isinstance(project.data, dict) or project.data.get("accession") != accession:
            raise ValueError("PRIDE project accession mismatch")
        native = accession
    else:
        native = _registry_native_id(project.data, repository, accession)
    if data.get("repository_accession") != native:
        raise ValueError("pridepy native project identity mismatch")
    files, seen = [], {}
    for index, record in enumerate(data["records"]):
        if not isinstance(record, dict) or not isinstance(record.get("fileName"), str) or not record["fileName"]:
            raise ValueError("pridepy file record has no valid filename")
        locations = record.get("publicFileLocations") or []
        if not isinstance(locations, list):
            raise ValueError("pridepy file record has invalid public locations")
        urls = []
        unresolved = False
        for location in locations:
            url = location.get("value") if isinstance(location, dict) else None
            if not isinstance(url, str) or not url:
                raise ValueError("pridepy file location has no valid URL")
            uri = urlsplit(url)
            if uri.scheme in ("http", "https", "ftp", "ftps"):
                if not uri.hostname or uri.username is not None or uri.password is not None:
                    raise ValueError("pridepy file location is invalid")
                if repository != "PRIDE":
                    identity_path = unquote(uri.path + "?" + uri.query)
                    if not re.search(r"(?<![A-Za-z0-9])" + re.escape(native) + r"(?![A-Za-z0-9])", identity_path):
                        raise ValueError("pridepy file URL belongs to another dataset")
                try:
                    path = unquote_to_bytes(uri.path).decode("utf-8", errors="strict")
                except UnicodeDecodeError:
                    unresolved = True
                else:
                    if any(part in (".", "..") for part in path.split("/")):
                        raise ValueError("pridepy file location has an invalid path")
            urls.append(url)
        if repository == "PRIDE":
            if accession not in (record.get("projectAccessions") or []):
                raise ValueError("PRIDE file belongs to another project")
            file_id = record.get("accession")
            if not isinstance(file_id, str) or not file_id:
                raise ValueError("PRIDE file has no file accession")
        else:
            if record.get("accession") != native:
                raise ValueError("pridepy file has a different native dataset accession")
            # Several providers put the dataset accession in every file record.
            identity = json.dumps([record.get("relativePath"), sorted(urls)], ensure_ascii=False, separators=(",", ":"))
            file_id = native + ":" + hashlib.sha256(identity.encode("utf-8")).hexdigest()
        if file_id in seen:
            old_index, original_record = seen[file_id]
            if original_record != record:
                raise ValueError("Repeated pridepy file has conflicting metadata")
            files[old_index]["providerRecordIndices"].append(index)
            continue
        seen[file_id] = (len(files), record)
        normalized = json.loads(json.dumps(record, allow_nan=False))
        size = normalized.get("fileSizeBytes")
        if size is not None and (type(size) is not int or size < 0):
            raise ValueError("pridepy file size must be a non-negative integer or unknown")
        normalized.update({"accession": file_id, "projectAccessions": [accession],
                           "repositoryAccessions": [native], "fileSizeBytes": size,
                           "providerRecordIndices": [index]})
        if unresolved:
            normalized["pathDecodingStatus"] = "unresolved"
        files.append(normalized)
    return files, None


def download_pilot_file(item, destination, file_limit, reservation):
    """Let pridepy transfer a known-size planned file into an isolated staging directory.

    Limits are preflight plan limits, not limits on network traffic or retries.
    No raw HTTP response metadata is exposed by pridepy's public download API.
    """
    destination = Path(destination)
    staging = destination.with_name(destination.name + ".part")
    expected = item.get("reported_size_bytes")
    reported_url = item["download_urls"][0]
    transfer_url, rewrite_reason = _pilot_transfer_url(item, reported_url)
    result = {
        "pilot_file_id": item["pilot_file_id"], "project_accession": item["project_accession"],
        "repository": item["repository"], "reported_file_name": item["reported_file_name"],
        "pilot_role": item["pilot_role"], "requested_url": reported_url, "transfer_url": transfer_url,
        "url_rewrite_reason": rewrite_reason, "response_url": "", "response_headers": {},
        "response_metadata_available": False, "transport_attempts": [],
        "download_backend": "pridepy", "backend_version": PRIDEPY_VERSION,
        "accounting_mode": "planned_file_sizes", "download_status": "failed", "local_path": "",
        "confirmed_payload_bytes": 0, "charged_payload_bytes": 0, "file_sha256": None,
        "observed_file_size_bytes": None, "error_type": "", "error_message": "",
        "content_probe": {"status": "not_run"},
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    if expected is None:
        result.update(download_status="skipped_unknown_size", error_type="UnknownFileSize",
                      error_message="File size is unknown; no transfer was started")
    elif type(expected) is not int or expected <= 0:
        result.update(error_type="InvalidFileSize", error_message="Reported file size must be positive")
    elif expected > file_limit:
        result.update(download_status="skipped_file_limit", error_type="PlannedFileLimit",
                      error_message="Reported file size exceeds the per-file plan limit")
    elif expected > reservation:
        result.update(download_status="budget_exhausted", error_type="PlannedBudgetLimit",
                      error_message="Reported file size exceeds the remaining plan budget")
    else:
        created_staging = False
        try:
            _pridepy_version()
            if destination.exists() or staging.exists():
                raise ValueError("An existing payload or staging directory must not be overwritten")
            destination.parent.mkdir(parents=True, exist_ok=True)
            staging.mkdir()
            created_staging = True
            result["charged_payload_bytes"] = expected
            attempt = {"scope": "library_call", "internal_attempt_count": None,
                       "transfer_url": transfer_url, "backend": "pridepy"}
            result["transport_attempts"].append(attempt)
            _pridepy_download(transfer_url, staging)
            produced = list(staging.iterdir())
            if len(produced) != 1 or produced[0].is_symlink() or not produced[0].is_file():
                raise ValueError("pridepy did not produce exactly one regular payload file")
            payload = produced[0]
            size = payload.stat().st_size
            result["observed_file_size_bytes"] = size
            if size != expected:
                raise ValueError("Downloaded file size differs from the reported size")
            digest = hashlib.sha256()
            with payload.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            checksum = digest.hexdigest()
            if item.get("expected_download_file_sha256") not in (None, checksum):
                raise ValueError("Downloaded file SHA256 differs from the expected checksum")
            payload.replace(destination)
            result.update(download_status="downloaded", local_path=str(destination.resolve()),
                          confirmed_payload_bytes=size, file_sha256=checksum)
            attempt["outcome"] = "downloaded"
            try:
                result["content_probe"] = probe_pilot_content(destination.read_bytes(), item["pilot_role"])
            except (OSError, ValueError) as error:
                result["content_probe"] = {"status": "probe_error", "error": str(error)}
        except Exception as error:
            result["error_type"], result["error_message"] = type(error).__name__, str(error)
            if result["transport_attempts"]:
                result["transport_attempts"][-1]["outcome"] = "failed"
        finally:
            if created_staging:
                # The explicit-URL API writes only flat files into this private directory.
                for path in staging.iterdir():
                    if path.is_symlink() or path.is_file():
                        path.unlink()
                staging.rmdir()
    result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    return result
