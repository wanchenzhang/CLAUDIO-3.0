
"""Bounded official publication requests with checksummed offline replay.

This layer verifies HTTP evidence, not a dataset-to-paper association.
An injected session must make exactly one HTTP attempt per get() call.
Byte budgets bound captured decoded response bodies, not network traffic.
"""

import hashlib
import json
import math
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
import time
from urllib.parse import quote, urlencode
from uuid import uuid4

import requests

from .publications import normalize_doi, normalize_pmid
from .publication_verification import _as_of_date, qualify_crossref_publication
from .pubmed_verification import qualify_pubmed_publication

MAX_BODY = 2 * 1024 * 1024
RETRY_STATUSES = {429, 500, 502, 503, 504}


def _request_identity(doi, pmid):
    normalized_doi = normalize_doi(doi) if doi is not None else None
    normalized_pmid = normalize_pmid(pmid) if pmid is not None else None
    if ((doi is not None and normalized_doi is None)
            or (pmid is not None and normalized_pmid is None)
            or (normalized_doi is None and normalized_pmid is None)):
        raise ValueError("A valid publication DOI or PMID is required")
    if normalized_pmid is not None:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode({
            "db": "pubmed", "id": normalized_pmid, "retmode": "xml",
            "rettype": "abstract", "tool": "claudio",
        })
        return "pubmed", normalized_pmid, url, normalized_doi, normalized_pmid
    url = "https://api.crossref.org/works/" + quote(normalized_doi, safe="")
    return "crossref", normalized_doi, url, normalized_doi, None


def _stamp():
    return datetime.now(timezone.utc).isoformat()


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _dump_new(path, value):
    import os
    raw = (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())
    return _digest(raw)


def _unresolved(provider, doi, pmid, reason, evidence=None):
    return {
        "provider": provider, "requested_doi": doi, "requested_pmid": pmid,
        "verification_scope": "captured_http_metadata", "status": "unresolved",
        "reason": reason, "paper": None, "dataset_association_verified": False,
        "evidence": evidence,
    }


def _bound_file(root, name, limit):
    if not isinstance(name, str) or not name or "\\" in name:
        raise ValueError("Invalid evidence path")
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("Evidence path must stay inside its directory")
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root) or not resolved.is_file():
        raise ValueError("Evidence path is missing or outside its directory")
    if resolved.stat().st_size > limit:
        raise ValueError("Evidence file exceeds the replay size limit")
    return resolved.read_bytes()


def replay_publication_lookup(evidence_dir, receipt_path, *, receipt_sha256,
                              doi=None, pmid=None, as_of):
    """Requalify pinned HTTP evidence without making any network requests."""
    provider, identifier, url, doi, pmid = _request_identity(doi, pmid)
    _as_of_date(as_of)
    root = Path(evidence_dir).resolve()
    raw = _bound_file(root, receipt_path, 128 * 1024)
    if _digest(raw) != receipt_sha256:
        raise ValueError("Publication receipt checksum mismatch")
    try:
        receipt = json.loads(raw)
    except (ValueError, UnicodeError) as error:
        raise ValueError("Invalid publication receipt JSON") from error
    if not isinstance(receipt, dict) or receipt.get("schema_version") != 1 or receipt.get("kind") != "publication_http_receipt":
        raise ValueError("Unsupported publication receipt")
    if (receipt.get("provider"), receipt.get("identifier"), receipt.get("request_url")) != (provider, identifier, url):
        raise ValueError("Publication receipt does not match the requested identifier")
    attempts = receipt.get("attempts")
    if not isinstance(attempts, list) or not 1 <= len(attempts) <= 3:
        raise ValueError("Invalid publication receipt attempts")
    body = None
    for attempt in attempts:
        if not isinstance(attempt, dict) or attempt.get("request_url") != url:
            raise ValueError("Invalid publication request in receipt")
        body = _bound_file(root, attempt.get("body_path"), MAX_BODY)
        if _digest(body) != attempt.get("body_sha256") or len(body) != attempt.get("body_bytes"):
            raise ValueError("Publication response body checksum or length mismatch")
    evidence = {"receipt_path": Path(receipt_path).as_posix(), "receipt_sha256": receipt_sha256}
    if receipt.get("failure_reason"):
        return _unresolved(provider, doi, pmid, receipt["failure_reason"], evidence)
    final = attempts[-1]
    if (final.get("status_code") != 200 or final.get("response_url") != url
            or final.get("complete") is not True or final.get("error") is not None):
        raise ValueError("Publication receipt does not contain a complete official response")
    if provider == "crossref":
        try:
            payload = json.loads(body)
        except (ValueError, UnicodeError):
            return _unresolved(provider, doi, pmid, "invalid_crossref_json", evidence)
        result = qualify_crossref_publication(payload, requested_doi=doi, as_of=as_of)
    else:
        result = qualify_pubmed_publication(body, requested_pmid=pmid, requested_doi=doi, as_of=as_of)
    return {**result, "verification_scope": "captured_http_metadata",
            "dataset_association_verified": False, "evidence": evidence}


class PublicationLookup:
    """Own one run's bounded request budget, immutable evidence, and raw cache."""

    def __init__(self, evidence_dir, *, session=None, max_requests=100,
                 max_response_bytes=MAX_BODY, max_total_bytes=20 * MAX_BODY,
                 max_attempts=3, max_elapsed_seconds=1800,
                 sleeper=time.sleep, clock=time.monotonic):
        for name, value, upper in (("max_requests", max_requests, 10000),
                                   ("max_response_bytes", max_response_bytes, MAX_BODY),
                                   ("max_total_bytes", max_total_bytes, 100 * MAX_BODY),
                                   ("max_attempts", max_attempts, 3)):
            if type(value) is not int or not 1 <= value <= upper:
                raise ValueError("Invalid " + name)
        if (isinstance(max_elapsed_seconds, bool) or not isinstance(max_elapsed_seconds, (int, float))
                or not math.isfinite(max_elapsed_seconds) or not 0 < max_elapsed_seconds <= 3600):
            raise ValueError("Invalid max_elapsed_seconds")
        self.root = Path(evidence_dir).resolve()
        self.root.mkdir(parents=True, exist_ok=False)
        (self.root / "bodies").mkdir()
        (self.root / "receipts").mkdir()
        self.session = session if session is not None else requests.Session()
        if session is None:
            self.session.mount("https://", requests.adapters.HTTPAdapter(max_retries=0))
        self.max_requests = max_requests
        self.max_response_bytes = max_response_bytes
        self.max_total_bytes = max_total_bytes
        self.max_attempts = max_attempts
        self.max_elapsed_seconds = max_elapsed_seconds
        self.request_count = 0
        self.captured_body_bytes = 0
        self._sleep = sleeper
        self._clock = clock
        self._started = clock()
        self._last_request = None
        self._cache = {}
        self._closed = False

    def close(self):
        if not self._closed:
            self.session.close()
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.close()

    def _budget_reason(self):
        if self.request_count >= self.max_requests:
            return "request_budget_exhausted"
        if self.captured_body_bytes >= self.max_total_bytes:
            return "total_byte_limit"
        if self._clock() - self._started >= self.max_elapsed_seconds:
            return "time_budget_exhausted"
        return None

    def _capture(self, url):
        import os
        body_path = "bodies/" + uuid4().hex + ".bin"
        attempt = {
            "request_url": url, "started_at_utc": _stamp(), "response_url": None,
            "status_code": None, "headers": {}, "body_path": body_path,
            "body_bytes": 0, "complete": False, "error": None, "error_type": None,
        }
        response = None
        digest = hashlib.sha256()
        self.request_count += 1
        self._last_request = self._clock()
        with (self.root / body_path).open("xb") as handle:
            try:
                response = self.session.get(
                    url, timeout=(10, 30), allow_redirects=False, stream=True,
                    headers={"User-Agent": "CLAUDIO-publication-lookup/1.0",
                             "Accept": "application/json, application/xml, text/xml"},
                )
                attempt["response_url"] = response.url
                attempt["status_code"] = response.status_code
                attempt["headers"] = {
                    str(key).lower(): str(value) for key, value in response.headers.items()
                    if str(key).lower() in {"content-type", "content-length", "retry-after", "date", "etag", "last-modified"}
                }
                cap = min(self.max_response_bytes, self.max_total_bytes - self.captured_body_bytes)
                for chunk in response.iter_content(chunk_size=min(16384, cap + 1)):
                    if self._clock() - self._started >= self.max_elapsed_seconds:
                        attempt["error"] = "time_budget_exhausted"
                        break
                    if not chunk:
                        continue
                    remaining = cap - attempt["body_bytes"]
                    saved = chunk[:remaining]
                    handle.write(saved)
                    digest.update(saved)
                    attempt["body_bytes"] += len(saved)
                    self.captured_body_bytes += len(saved)
                    if len(chunk) > remaining:
                        attempt["error"] = ("total_byte_limit" if self.captured_body_bytes >= self.max_total_bytes else "response_byte_limit")
                        break
                else:
                    attempt["complete"] = True
                if response.url != url:
                    attempt["error"] = "unexpected_response_url"
            except requests.exceptions.RequestException as error:
                attempt["error"] = "request_failed"
                attempt["error_type"] = type(error).__name__
            finally:
                if response is not None:
                    response.close()
                handle.flush()
                os.fsync(handle.fileno())
        attempt["body_sha256"] = digest.hexdigest()
        attempt["completed_at_utc"] = _stamp()
        return attempt

    @staticmethod
    def _retry_delay(attempt, number):
        raw = attempt["headers"].get("retry-after")
        if raw is not None:
            try:
                seconds = float(raw)
                if not math.isfinite(seconds) or seconds < 0:
                    raise ValueError("Invalid Retry-After")
                return seconds
            except ValueError:
                try:
                    date = parsedate_to_datetime(raw)
                    if date.tzinfo is None:
                        date = date.replace(tzinfo=timezone.utc)
                    return max(0.0, (date - datetime.now(timezone.utc)).total_seconds())
                except (TypeError, ValueError, OverflowError):
                    pass
        return float(2 ** number)

    def lookup(self, *, doi=None, pmid=None, as_of):
        """Capture an official response, then qualify its pinned stored bytes."""
        if self._closed:
            raise ValueError("Publication lookup is closed")
        provider, identifier, url, doi, pmid = _request_identity(doi, pmid)
        _as_of_date(as_of)
        key = (provider, identifier)
        if key in self._cache:
            path, digest = self._cache[key]
            return replay_publication_lookup(self.root, path, receipt_sha256=digest, doi=doi, pmid=pmid, as_of=as_of)
        attempts = []
        reason = None
        for number in range(self.max_attempts):
            reason = self._budget_reason()
            if reason:
                break
            if self._last_request is not None:
                delay = 0.4 - (self._clock() - self._last_request)
                if delay > 0:
                    self._sleep(delay)
            reason = self._budget_reason()
            if reason:
                break
            attempt = self._capture(url)
            attempts.append(attempt)
            reason = attempt["error"]
            status = attempt["status_code"]
            if reason is None and status == 200:
                break
            if reason is None:
                reason = "redirect_not_followed" if isinstance(status, int) and 300 <= status < 400 else "http_status_error"
            retry = reason == "request_failed" or (reason == "http_status_error" and status in RETRY_STATUSES)
            if not retry or number + 1 >= self.max_attempts:
                break
            budget = self._budget_reason()
            if budget:
                reason = budget
                break
            delay = self._retry_delay(attempt, number)
            if delay > 30:
                reason = "retry_wait_exceeds_limit"
                break
            if self._clock() - self._started + delay >= self.max_elapsed_seconds:
                reason = "time_budget_exhausted"
                break
            self._sleep(delay)
        if not attempts:
            return _unresolved(provider, doi, pmid, reason)
        receipt = {
            "schema_version": 1, "kind": "publication_http_receipt",
            "provider": provider, "identifier": identifier, "request_url": url,
            "attempts": attempts, "failure_reason": reason, "completed_at_utc": _stamp(),
        }
        path = "receipts/" + uuid4().hex + ".json"
        digest = _dump_new(self.root / path, receipt)
        if reason is None:
            self._cache[key] = (path, digest)
        return replay_publication_lookup(self.root, path, receipt_sha256=digest, doi=doi, pmid=pmid, as_of=as_of)
