
"""Offline migration contracts for original raw discovery snapshots."""

import copy
import hashlib
import importlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from click.testing import CliRunner
import requests

from claudio.dataset_catalog.discovery_history import DiscoveryHistory, initialize_discovery_history
from claudio.dataset_catalog.proteomecentral import CatalogPage, write_discovery_outputs
from claudio.dataset_catalog.publication_lookup import PublicationLookup as RealLookup
from test_dataset_catalog_publication_lookup import Clock, crossref_body
from test_dataset_catalog_published_discovery import ProviderSession, record, page


DOI = "10.1234/example"
LEGACY_DATE = "2026-10-31T12:00:00+00:00"


class TestDiscoverySeeding(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.history = self.root / "history"
        initialize_discovery_history(self.history)
        self.calls = []
        guard = patch(
            "requests.sessions.Session.request",
            side_effect=AssertionError(
                "Offline tests prohibit real network requests"
            ),
        )
        guard.start()
        self.addCleanup(guard.stop)

    def api(self):
        try:
            return importlib.import_module(
                "claudio.dataset_catalog.discovery_seed"
            )
        except ModuleNotFoundError as error:
            if error.name != "claudio.dataset_catalog.discovery_seed":
                raise
            self.fail(
                "Legacy discovery history seeding is not implemented"
            )

    def legacy(self, records=None, name="legacy", created=LEGACY_DATE):
        records = [record()] if records is None else records
        source = page("XL-MS", records)
        source = CatalogPage(
            source.request_url, created, source.raw_body, source.data
        )
        destination = self.root / name
        write_discovery_outputs({"XL-MS": [source]}, destination)
        manifest = json.loads(
            (destination / "manifest.json").read_bytes()
        )
        manifest["created_at_utc"] = created
        (destination / "manifest.json").write_text(
            json.dumps(manifest) + "\n", encoding="utf-8"
        )
        return destination

    def factory(self, responses=None):
        def create(directory, **options):
            clock = Clock()
            return RealLookup(
                directory,
                session=ProviderSession(self.calls, responses),
                clock=clock,
                sleeper=clock.sleep,
                **options,
            )
        return create

    def seed(self, source, name="seed", responses=None, **options):
        api = self.api()
        with patch.object(
            api, "PublicationLookup",
            side_effect=self.factory(responses),
        ):
            return api.seed_legacy_discovery(
                source,
                self.root / name,
                history_dir=self.history,
                **options,
            )

    def read(self, directory, name):
        return json.loads((Path(directory) / name).read_bytes())

    def seen(self):
        with DiscoveryHistory(self.history) as history:
            return history.seen_accessions()

    def test_seeds_only_original_candidates_published_by_the_legacy_date(self):
        output = self.seed(
            self.legacy([record(), record("PXD000002", "Pending")])
        )
        self.assertEqual(self.seen(), {"PXD000001"})
        report = self.read(output, "seed_report.json")
        self.assertEqual(
            report["seeded_accessions"], ["PXD000001"]
        )
        self.assertEqual(
            report["not_seeded_accessions"], ["PXD000002"]
        )
        self.assertEqual(
            report["publication_as_of"], "2026-10-31"
        )
        self.assertEqual(report["source_candidates"], 2)
        self.assertEqual(len(self.calls), 1)

    def test_repeated_seeding_does_not_reemit_or_requery_seen_candidates(self):
        source = self.legacy()
        first = self.seed(source, "first")
        second = self.seed(source, "second")
        self.assertEqual(
            self.read(first, "seed_report.json")["seeded_count"], 1
        )
        self.assertEqual(
            self.read(second, "seed_report.json")["seeded_count"], 0
        )
        self.assertEqual(self.seen(), {"PXD000001"})
        self.assertEqual(len(self.calls), 1)

    def test_paper_published_after_the_legacy_run_is_not_backdated(self):
        payload = json.loads(crossref_body())
        payload["message"]["published-online"] = {
            "date-parts": [[2026, 11, 5]]
        }
        output = self.seed(
            self.legacy(), responses={DOI: payload}
        )
        self.assertEqual(
            self.read(output, "seed_report.json")["seeded_count"], 0
        )
        self.assertEqual(self.seen(), set())
        self.assertEqual(
            self.read(output, "seed_audit.json")[0][
                "publication_status"
            ],
            "ineligible",
        )

    def test_previously_pending_dataset_can_enter_a_later_fresh_discovery(self):
        source = self.legacy(
            [record(publications="Pending")]
        )
        self.seed(source)
        self.assertEqual(self.seen(), set())
        from claudio.dataset_catalog import published_discovery
        with patch.object(
            published_discovery,
            "fetch_dataset_pages",
            return_value=[page("XL-MS", [record()])],
        ), patch.object(
            published_discovery,
            "PublicationLookup",
            side_effect=self.factory(),
        ):
            output = published_discovery.run_published_discovery(
                ["XL-MS"],
                self.root / "new",
                history_dir=self.history,
                as_of="2026-11-09",
            )
        self.assertEqual(
            len(self.read(output, "candidates.json")), 1
        )

    def test_seeding_does_not_add_papers_missing_from_the_original_sources(self):
        source = self.legacy([record(publications="")])
        output = self.seed(source)
        self.assertEqual(
            self.read(output, "seed_report.json")["seeded_count"], 0
        )
        self.assertEqual(self.calls, [])
        self.assertEqual(self.seen(), set())

    def test_original_snapshot_bytes_are_unchanged(self):
        source = self.legacy()
        before = {
            p.relative_to(source).as_posix(): p.read_bytes()
            for p in source.rglob("*") if p.is_file()
        }
        self.seed(source)
        after = {
            p.relative_to(source).as_posix(): p.read_bytes()
            for p in source.rglob("*") if p.is_file()
        }
        self.assertEqual(after, before)

    def test_changed_raw_response_is_rejected_before_lookup(self):
        source = self.legacy()
        raw = next((source / "raw").iterdir())
        raw.write_bytes(raw.read_bytes() + b"\n")
        with self.assertRaisesRegex(
            ValueError, "(?i)checksum"
        ):
            self.seed(source)
        self.assertEqual(self.calls, [])
        self.assertEqual(self.seen(), set())

    def test_candidate_json_must_equal_the_recomputed_raw_search_merge(self):
        source = self.legacy()
        rows = self.read(source, "candidates.json")
        rows[0]["publications"] = (
            '<a href="https://doi.org/10.1234/added">'
            'Added later</a>'
        )
        (source / "candidates.json").write_text(
            json.dumps(rows), encoding="utf-8"
        )
        with self.assertRaisesRegex(
            ValueError, "(?i)candidate|merge"
        ):
            self.seed(source)
        self.assertEqual(self.calls, [])

    def test_manifest_candidate_count_must_match(self):
        source = self.legacy()
        manifest = self.read(source, "manifest.json")
        manifest["candidate_count"] = 9
        (source / "manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        with self.assertRaisesRegex(
            ValueError, "(?i)count"
        ):
            self.seed(source)

    def test_inventory_copy_is_not_accepted_as_an_original_discovery(self):
        source = self.legacy()
        (source / "manifest.json").write_text(
            json.dumps({
                "schema_version": 1,
                "source": "PRIDE",
                "status": "complete",
            }),
            encoding="utf-8",
        )
        with self.assertRaisesRegex(
            ValueError, "(?i)original|discovery"
        ):
            self.seed(source)
        self.assertEqual(self.calls, [])

    def test_incomplete_original_manifest_is_rejected(self):
        source = self.legacy()
        manifest = self.read(source, "manifest.json")
        manifest["status"] = "failed"
        (source / "manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        with self.assertRaises(ValueError):
            self.seed(source)
        self.assertEqual(self.seen(), set())

    def test_source_paths_cannot_escape_the_legacy_snapshot(self):
        source = self.legacy()
        manifest = self.read(source, "manifest.json")
        manifest["searches"][0]["pages"][0]["path"] = (
            "../outside.json"
        )
        (source / "manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        with self.assertRaisesRegex(
            ValueError, "(?i)path|snapshot"
        ):
            self.seed(source)

    def test_timezone_is_required_for_the_legacy_run_date(self):
        source = self.legacy(
            created="2026-10-31T12:00:00"
        )
        with self.assertRaisesRegex(
            ValueError, "(?i)date|timezone"
        ):
            self.seed(source)
        self.assertEqual(self.calls, [])

    def test_future_dated_source_snapshot_is_rejected(self):
        source = self.legacy(
            created="2099-01-01T12:00:00+00:00"
        )
        with self.assertRaisesRegex(
            ValueError, "(?i)future|date"
        ):
            self.seed(source)

    def test_unresolved_http_lookups_do_not_consume_history(self):
        source = self.legacy()
        output = self.seed(
            source,
            responses={
                DOI: requests.exceptions.Timeout("Offline")
            },
            publication_request_limit=1,
        )
        self.assertEqual(
            self.read(output, "seed_report.json")["seeded_count"], 0
        )
        self.assertEqual(
            self.read(output, "seed_audit.json")[0][
                "publication_status"
            ],
            "unresolved",
        )
        self.assertEqual(self.seen(), set())

    def test_preprints_are_not_seeded(self):
        payload = json.loads(crossref_body())
        payload["message"]["type"] = "posted-content"
        output = self.seed(
            self.legacy(), responses={DOI: payload}
        )
        self.assertEqual(
            self.read(output, "seed_report.json")["seeded_count"], 0
        )
        self.assertEqual(self.seen(), set())

    def test_paper_doi_is_not_the_history_key(self):
        output = self.seed(
            self.legacy([record(), record("PXD000002")])
        )
        self.assertEqual(
            self.read(output, "seed_report.json")["seeded_count"], 2
        )
        self.assertEqual(
            self.seen(), {"PXD000001", "PXD000002"}
        )
        self.assertEqual(len(self.calls), 1)

    def test_missing_history_is_not_auto_initialized(self):
        api = self.api()
        with self.assertRaisesRegex(
            ValueError, "(?i)initializ|missing"
        ):
            api.seed_legacy_discovery(
                self.legacy(),
                self.root / "seed",
                history_dir=self.root / "missing",
            )
        self.assertFalse((self.root / "missing").exists())

    def test_existing_output_is_preserved(self):
        source = self.legacy()
        output = self.root / "seed"
        output.mkdir()
        (output / "keep").write_bytes(b"keep")
        with self.assertRaises(FileExistsError):
            self.seed(source)
        self.assertEqual(
            (output / "keep").read_bytes(), b"keep"
        )
        self.assertEqual(self.calls, [])

    def test_report_and_manifest_pin_source_and_publication_evidence(self):
        source = self.legacy()
        output = self.seed(source)
        manifest = self.read(output, "manifest.json")
        report = self.read(output, "seed_report.json")
        self.assertEqual(
            manifest["operation"], "legacy_history_seed"
        )
        self.assertEqual(
            report["source_manifest_sha256"],
            hashlib.sha256(
                (source / "manifest.json").read_bytes()
            ).hexdigest(),
        )
        self.assertEqual(
            report["source_candidates_sha256"],
            hashlib.sha256(
                (source / "candidates.json").read_bytes()
            ).hexdigest(),
        )
        self.assertIn(
            "legacy_source/manifest.json",
            manifest["output_sha256"],
        )
        self.assertIn(
            "seed_audit.json", manifest["output_sha256"]
        )
        self.assertTrue(any(
            name.startswith("publication_evidence/bodies/")
            for name in manifest["output_sha256"]
        ))
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(
                hashlib.sha256(
                    (output / name).read_bytes()
                ).hexdigest(),
                digest,
            )

    def test_output_write_failure_does_not_consume_history(self):
        api = self.api()
        real_write = api._write_json

        def fail_report(path, value):
            if Path(path).name == "seed_report.json":
                raise OSError("Disk failure")
            return real_write(path, value)

        with patch.object(
            api, "_write_json", side_effect=fail_report
        ):
            with self.assertRaises(OSError):
                self.seed(self.legacy())
        self.assertEqual(self.seen(), set())
        self.assertFalse(
            (self.root / "seed/manifest.json").exists()
        )

    def test_cli_can_seed_validated_legacy_results(self):
        api = self.api()
        from claudio.dataset_catalog.cli import main
        source = self.legacy()
        with patch.object(
            api,
            "PublicationLookup",
            side_effect=self.factory(),
        ):
            result = CliRunner().invoke(main, [
                "seed-discovery-history",
                "--source-run", str(source),
                "--history", str(self.history),
                "--output", str(self.root / "seed"),
            ])
        self.assertEqual(
            result.exit_code, 0,
            f"{result.output}\n{result.exception}",
        )
        self.assertEqual(self.seen(), {"PXD000001"})
        self.assertIn(
            "Seeded historical published datasets: 1",
            result.output,
        )

    def test_cli_rejects_invalid_publication_budget(self):
        from claudio.dataset_catalog.cli import main
        result = CliRunner().invoke(main, [
            "seed-discovery-history",
            "--source-run", str(self.legacy()),
            "--history", str(self.history),
            "--output", str(self.root / "seed"),
            "--publication-request-limit", "0",
        ])
        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn(
            "--publication-request-limit", result.output
        )


if __name__ == "__main__":
    unittest.main()
