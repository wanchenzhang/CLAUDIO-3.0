
"""Persistent discovery emission history with immutable receipts and recovery.

Publication qualification and dataset-to-paper association are caller duties.
Only caller-qualified published candidates with paper metadata may be emitted.
Network requests must finish before publish() acquires the history transaction.
"""

from contextlib import contextmanager
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
from uuid import UUID, uuid4


SCHEMA_VERSION = 1
MANIFEST_KIND = "published_incremental_discovery"
STATUSES = {"published", "unpublished", "unverified", "pending", "unresolved", "ineligible"}


def default_history_directory():
    """Return the state directory without creating or initializing it."""
    value = os.environ.get("XDG_STATE_HOME")
    base = Path(value).expanduser() if value else Path.home() / ".local/state"
    if not base.is_absolute():
        raise ValueError("XDG_STATE_HOME must be an absolute directory")
    return base / "claudio/discovery"


def _stamp():
    return datetime.now(timezone.utc).isoformat()


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sync_directory(path):
    if os.name == "posix":
        descriptor = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


def _json_new(path, value):
    """Expose an immutable JSON receipt only after all bytes are synced."""
    path = Path(path)
    raw = (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    # Keep the temporary entry outside the directory whose files are checksummed.
    temporary = path.parent.parent / (".claudio-json-" + uuid4().hex + ".tmp")
    created = False
    try:
        with temporary.open("xb") as handle:
            created = True
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        # Linking publishes the complete inode atomically and refuses overwrite.
        os.link(temporary, path)
        _sync_directory(path.parent)
    finally:
        if created:
            temporary.unlink(missing_ok=True)


def _json_read(path, limit):
    path = Path(path)
    if not path.is_file() or path.stat().st_size > limit:
        raise ValueError("Missing or oversized JSON evidence: " + str(path))
    try:
        return json.loads(path.read_bytes())
    except (ValueError, UnicodeError) as error:
        raise ValueError("Invalid JSON evidence: " + str(path)) from error


def _accession(value):
    if not isinstance(value, str):
        raise ValueError("A dataset accession must be a string")
    value = value.strip().upper()
    if re.fullmatch(r"[A-Z][A-Z0-9._-]*", value) is None:
        raise ValueError("Invalid dataset accession")
    return value


def _eligible(candidates):
    if not isinstance(candidates, list):
        raise ValueError("Discovery candidates must be a list")
    accessions = set()
    eligible = []
    for candidate in candidates:
        if not isinstance(candidate, dict):
            raise ValueError("A candidate must be an object")
        accession = _accession(candidate.get("accession"))
        if accession in accessions:
            raise ValueError("Duplicate dataset accession")
        accessions.add(accession)
        status = candidate.get("publication_status")
        if not isinstance(status, str) or status not in STATUSES:
            raise ValueError("Candidate publication status is missing or unsupported")
        if status == "published":
            papers = candidate.get("papers")
            if not isinstance(papers, list) or not papers or any(not isinstance(paper, dict) or not paper for paper in papers):
                raise ValueError("Published candidates must carry paper metadata")
            eligible.append(deepcopy(candidate))
    json.dumps(eligible, allow_nan=False)
    return eligible


def initialize_discovery_history(directory):
    """Explicitly create a new namespace; never reset an existing directory."""
    directory = Path(directory).expanduser().absolute()
    directory.mkdir(parents=True, exist_ok=False)
    directory = directory.resolve()
    namespace = str(uuid4())
    connection = sqlite3.connect(directory / "history.sqlite3")
    try:
        connection.execute("PRAGMA journal_mode=DELETE")
        connection.execute("PRAGMA synchronous=FULL")
        connection.execute("PRAGMA foreign_keys=ON")
        connection.executescript("""
            CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
            CREATE TABLE runs (
                run_id TEXT PRIMARY KEY,
                output_directory TEXT NOT NULL UNIQUE,
                status TEXT NOT NULL CHECK(status IN ('prepared', 'complete')),
                manifest_sha256 TEXT,
                created_at_utc TEXT NOT NULL
            );
            CREATE TABLE seen_accessions (
                accession TEXT PRIMARY KEY,
                run_id TEXT NOT NULL REFERENCES runs(run_id)
            );
            PRAGMA user_version=1;
        """)
        connection.executemany("INSERT INTO metadata(key, value) VALUES (?, ?)",
                               [("namespace_id", namespace), ("schema_version", str(SCHEMA_VERSION))])
        connection.commit()
    finally:
        connection.close()
    marker = {"kind": "claudio_discovery_history_namespace", "schema_version": SCHEMA_VERSION,
              "namespace_id": namespace, "created_at_utc": _stamp()}
    _json_new(directory / "namespace.json", marker)
    _sync_directory(directory)
    _sync_directory(directory.parent)
    return {"history_directory": str(directory), "namespace_id": namespace}


class DiscoveryHistory:
    """Serialize final selection and record only completed, verified outputs."""

    def __init__(self, directory):
        self.directory = Path(directory).expanduser().resolve()
        database = self.directory / "history.sqlite3"
        marker_path = self.directory / "namespace.json"
        if not database.is_file() or not marker_path.is_file():
            raise ValueError("Discovery history is missing; explicit initialization is required")
        self._connection = None
        try:
            marker = _json_read(marker_path, 65536)
            if (not isinstance(marker, dict) or marker.get("kind") != "claudio_discovery_history_namespace"
                    or marker.get("schema_version") != SCHEMA_VERSION):
                raise ValueError("Unsupported history namespace marker")
            namespace = marker.get("namespace_id")
            if not isinstance(namespace, str) or str(UUID(namespace)) != namespace:
                raise ValueError("Invalid history namespace identifier")
            self._connection = sqlite3.connect(database.as_uri() + "?mode=rw", uri=True,
                                               timeout=5, isolation_level=None)
            self._connection.row_factory = sqlite3.Row
            self._connection.execute("PRAGMA busy_timeout=5000")
            self._connection.execute("PRAGMA foreign_keys=ON")
            check = self._connection.execute("PRAGMA quick_check").fetchall()
            if len(check) != 1 or check[0][0] != "ok":
                raise ValueError("Discovery history database integrity check failed")
            metadata = dict(self._connection.execute("SELECT key, value FROM metadata").fetchall())
            version = self._connection.execute("PRAGMA user_version").fetchone()[0]
            if version != SCHEMA_VERSION or metadata.get("schema_version") != str(SCHEMA_VERSION):
                raise ValueError("Unsupported discovery history schema")
            if metadata.get("namespace_id") != namespace:
                raise ValueError("History namespace marker and database disagree")
            expected = {"runs": {"run_id", "output_directory", "status", "manifest_sha256", "created_at_utc"},
                        "seen_accessions": {"accession", "run_id"}}
            for table, fields in expected.items():
                actual = {row[1] for row in self._connection.execute("PRAGMA table_info(" + table + ")")}
                if actual != fields:
                    raise ValueError("Unexpected discovery history table schema")
            self._connection.execute("PRAGMA journal_mode=DELETE")
            self._connection.execute("PRAGMA synchronous=FULL")
            self.namespace_id = namespace
        except (sqlite3.Error, ValueError, OSError, TypeError, AttributeError) as error:
            self.close()
            raise ValueError("Cannot open discovery history: " + str(error)) from error

    def close(self):
        if self._connection is not None:
            self._connection.close()
            self._connection = None

    def __enter__(self):
        return self

    def __exit__(self, *exception):
        self.close()

    @contextmanager
    def _transaction(self):
        if self._connection is None:
            raise ValueError("Discovery history is closed")
        self._connection.execute("BEGIN IMMEDIATE")
        try:
            yield
            self._connection.commit()
        except BaseException:
            self._connection.rollback()
            raise

    @staticmethod
    def _files(directory, *, sync=False):
        files = {}
        directories = [directory]
        if directory.is_symlink() or not directory.is_dir() or directory.resolve() != directory:
            raise ValueError("Output directory changed or contains a symlink")
        for path in sorted(directory.rglob("*")):
            if path.is_symlink():
                raise ValueError("Output evidence must not contain symlinks")
            if path.is_dir():
                directories.append(path)
            elif path.is_file():
                name = path.relative_to(directory).as_posix()
                if name == "manifest.json":
                    continue
                files[name] = _sha256(path)
                if sync:
                    with path.open("ab") as handle:
                        os.fsync(handle.fileno())
            else:
                raise ValueError("Unsupported output evidence file")
        if sync:
            for path in reversed(directories):
                _sync_directory(path)
        return files

    def _read_complete(self, row):
        directory = Path(row["output_directory"])
        manifest_path = directory / "manifest.json"
        manifest = _json_read(manifest_path, 2 * 1024 * 1024)
        if not isinstance(manifest, dict):
            raise ValueError("Invalid complete discovery receipt")
        if (manifest.get("schema_version") != 2 or manifest.get("kind") != MANIFEST_KIND
                or manifest.get("status") != "complete"
                or manifest.get("namespace_id") != self.namespace_id
                or manifest.get("run_id") != row["run_id"]
                or manifest.get("publication_policy") != "confirmed_published_only"):
            raise ValueError("Discovery receipt has a different run, namespace, or policy")
        recorded = manifest.get("output_sha256")
        if not isinstance(recorded, dict) or "candidates.json" not in recorded:
            raise ValueError("Discovery receipt is missing candidate output checksums")
        for name, digest in recorded.items():
            relative = Path(name) if isinstance(name, str) else None
            if (relative is None or not name or "\\" in name or relative.is_absolute()
                    or ".." in relative.parts or name == "manifest.json"
                    or not isinstance(digest, str) or re.fullmatch("[0-9a-f]{64}", digest) is None):
                raise ValueError("Invalid output checksum path in discovery receipt")
        if self._files(directory) != recorded:
            raise ValueError("Completed discovery output checksum mismatch")
        candidates = _json_read(directory / "candidates.json", 64 * 1024 * 1024)
        selected = _eligible(candidates)
        if len(selected) != len(candidates):
            raise ValueError("Complete output contains unqualified candidates")
        accessions = [_accession(candidate["accession"]) for candidate in selected]
        if (manifest.get("emitted_accessions") != accessions
                or type(manifest.get("candidate_count")) is not int
                or manifest["candidate_count"] != len(accessions)):
            raise ValueError("Discovery receipt candidate count or accession mismatch")
        digest = _sha256(manifest_path)
        if row["status"] == "complete" and row["manifest_sha256"] != digest:
            raise ValueError("Completed discovery manifest checksum mismatch")
        return manifest, digest

    def _record_completed(self, row, manifest, digest):
        for accession in manifest["emitted_accessions"]:
            previous = self._connection.execute(
                "SELECT run_id FROM seen_accessions WHERE accession=?", (accession,)).fetchone()
            if previous is not None:
                raise ValueError("Completed outputs conflict for dataset accession " + accession)
            self._connection.execute("INSERT INTO seen_accessions(accession, run_id) VALUES (?, ?)",
                                     (accession, row["run_id"]))
        self._connection.execute("UPDATE runs SET status='complete', manifest_sha256=? WHERE run_id=?",
                                 (digest, row["run_id"]))

    def _recover_locked(self):
        recovered = []
        verified = 0
        rows = self._connection.execute("SELECT * FROM runs ORDER BY rowid").fetchall()
        for row in rows:
            path = Path(row["output_directory"]) / "manifest.json"
            if row["status"] == "prepared" and not path.exists():
                continue
            manifest, digest = self._read_complete(row)
            if row["status"] == "prepared":
                self._record_completed(row, manifest, digest)
                recovered.append(row["run_id"])
            else:
                saved = {item[0] for item in self._connection.execute(
                    "SELECT accession FROM seen_accessions WHERE run_id=?", (row["run_id"],))}
                if saved != set(manifest["emitted_accessions"]):
                    raise ValueError("History accession entries differ from the completed receipt")
                verified += 1
        invalid = self._connection.execute("""
            SELECT COUNT(*) FROM seen_accessions AS seen
            LEFT JOIN runs AS run ON seen.run_id=run.run_id
            WHERE run.run_id IS NULL OR run.status!='complete'
        """).fetchone()[0]
        if invalid:
            raise ValueError("History has accession entries without completed outputs")
        return {"recovered_run_ids": recovered, "verified_completed_runs": verified}

    def recover(self):
        """Reconcile completed filesystem outputs before the next selection."""
        with self._transaction():
            return self._recover_locked()

    def seen_accessions(self):
        with self._transaction():
            self._recover_locked()
            return {row[0] for row in self._connection.execute("SELECT accession FROM seen_accessions")}

    def prepare_run(self, output_directory):
        """Durably reserve a new output path without consuming any candidates."""
        output = Path(output_directory).expanduser().resolve()
        if output.exists():
            raise FileExistsError("Discovery output directory already exists")
        if output == self.directory or self.directory in output.parents or output in self.directory.parents:
            raise ValueError("Discovery outputs and history must use separate directories")
        run_id = str(uuid4())
        with self._transaction():
            self._recover_locked()
            if self._connection.execute("SELECT 1 FROM runs WHERE output_directory=?", (str(output),)).fetchone():
                raise ValueError("Discovery output path is already reserved")
            self._connection.execute("""
                INSERT INTO runs(run_id, output_directory, status, manifest_sha256, created_at_utc)
                VALUES (?, ?, 'prepared', NULL, ?)
            """, (run_id, str(output), _stamp()))
        return {"run_id": run_id, "namespace_id": self.namespace_id, "output_directory": str(output)}

    def _summary(self, row, manifest, digest):
        return {"run_id": row["run_id"], "namespace_id": self.namespace_id,
                "output_directory": row["output_directory"], "candidate_count": manifest["candidate_count"],
                "emitted_accessions": deepcopy(manifest["emitted_accessions"]), "manifest_sha256": digest}

    def publish(self, run_id, candidates, write_outputs):
        """Select unseen qualified candidates, write files, then register them.

        The callback receives a new directory and a copy of the selected rows.
        It writes candidates.json and other evidence, and returns extra manifest
        fields. It must not write manifest.json or make network requests.
        Failed partial directories are retained and require a new run path.
        """
        eligible = _eligible(candidates)
        with self._transaction():
            self._recover_locked()
            row = self._connection.execute("SELECT * FROM runs WHERE run_id=?", (run_id,)).fetchone()
            if row is None:
                raise ValueError("Unknown prepared discovery run")
            if row["status"] == "complete":
                manifest, digest = self._read_complete(row)
                return self._summary(row, manifest, digest)
            seen = {item[0] for item in self._connection.execute("SELECT accession FROM seen_accessions")}
            selected = [candidate for candidate in eligible if _accession(candidate["accession"]) not in seen]
            directory = Path(row["output_directory"])
            directory.mkdir(parents=True, exist_ok=False)
            _sync_directory(directory.parent)
            extra = write_outputs(directory, deepcopy(selected))
            if not isinstance(extra, dict):
                raise ValueError("The discovery writer must return extra manifest fields")
            manifest_path = directory / "manifest.json"
            if manifest_path.exists():
                raise ValueError("The discovery writer must not write the complete manifest")
            written = _json_read(directory / "candidates.json", 64 * 1024 * 1024)
            if written != selected:
                raise ValueError("Written candidates do not match the selected records")
            manifest = {
                "schema_version": 2, "kind": MANIFEST_KIND, "status": "complete",
                "namespace_id": self.namespace_id, "run_id": run_id, "created_at_utc": _stamp(),
                "publication_policy": "confirmed_published_only",
                "candidate_count": len(selected),
                "emitted_accessions": [_accession(candidate["accession"]) for candidate in selected],
                "output_sha256": self._files(directory, sync=True),
            }
            if set(extra) & set(manifest):
                raise ValueError("The discovery writer tried to override reserved manifest fields")
            manifest.update(deepcopy(extra))
            _json_new(manifest_path, manifest)
            manifest, digest = self._read_complete(row)
            self._record_completed(row, manifest, digest)
            return self._summary(row, manifest, digest)
