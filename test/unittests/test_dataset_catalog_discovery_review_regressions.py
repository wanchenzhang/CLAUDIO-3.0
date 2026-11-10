"""Regression tests for final-review discovery boundaries."""

from datetime import datetime, timedelta, timezone
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from claudio.dataset_catalog import discovery_history, published_discovery
import test_dataset_catalog_discovery_history as history_tests
import test_dataset_catalog_published_discovery as discovery_tests
from test_dataset_catalog_publication_lookup import crossref_body


class InterruptedWrite:
    def __init__(self, handle):
        self.handle = handle

    def __enter__(self):
        self.handle.__enter__()
        return self

    def __exit__(self, *exception):
        return self.handle.__exit__(*exception)

    def __getattr__(self, name):
        return getattr(self.handle, name)

    def write(self, raw):
        self.handle.write(raw[:max(1, len(raw) // 3)])
        self.handle.flush()
        raise OSError("simulated_manifest_write_failure")


class TestDiscoveryReviewRegressions(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.history = self.root / "history"
        discovery_history.initialize_discovery_history(self.history)
        guard = patch(
            "requests.sessions.Session.request",
            side_effect=AssertionError("Live HTTP is forbidden"),
        )
        guard.start()
        self.addCleanup(guard.stop)

    def interrupted_publish(self, failure):
        original_json_new = discovery_history._json_new
        original_open = Path.open

        def interrupted_open(path, *arguments, **options):
            handle = original_open(path, *arguments, **options)
            mode = (
                arguments[0] if arguments
                else options.get("mode", "r")
            )
            return (
                InterruptedWrite(handle)
                if mode == "xb" else handle
            )

        def fail_manifest(path, value):
            if Path(path).name != "manifest.json":
                return original_json_new(path, value)
            if failure == "write":
                with patch.object(
                    Path, "open", new=interrupted_open
                ):
                    return original_json_new(path, value)
            with patch.object(
                discovery_history.os,
                "fsync",
                side_effect=OSError(
                    "simulated_manifest_fsync_failure"
                ),
            ):
                return original_json_new(path, value)

        output = self.root / "interrupted"
        with discovery_history.DiscoveryHistory(
            self.history
        ) as ledger:
            run = ledger.prepare_run(output)
            with patch.object(
                discovery_history,
                "_json_new",
                side_effect=fail_manifest,
            ):
                with self.assertRaisesRegex(
                    OSError, "simulated_manifest_"
                ):
                    ledger.publish(
                        run["run_id"],
                        [history_tests.candidate()],
                        history_tests.write_outputs,
                    )
            self.assertFalse(
                (output / "manifest.json").exists(),
                "Interrupted manifest write exposed the final receipt",
            )
            self.assertTrue(
                (output / "candidates.json").is_file()
            )

        with discovery_history.DiscoveryHistory(
            self.history
        ) as reopened:
            self.assertEqual(
                reopened.seen_accessions(), set()
            )
            retry = reopened.prepare_run(
                self.root / "next_run"
            )
            result = reopened.publish(
                retry["run_id"],
                [history_tests.candidate()],
                history_tests.write_outputs,
            )
            self.assertEqual(
                result["candidate_count"], 1
            )
            self.assertEqual(
                reopened.seen_accessions(),
                {"PXD000001"},
            )

    def test_interrupted_manifest_write_does_not_poison_history(self):
        self.interrupted_publish("write")

    def test_failed_manifest_fsync_does_not_publish_a_receipt(self):
        self.interrupted_publish("fsync")

    def test_json_publication_never_overwrites_existing_receipt(self):
        directory = self.root / "immutable"
        directory.mkdir()
        path = directory / "manifest.json"
        original = b'{"immutable":"original"}\n'
        path.write_bytes(original)
        with self.assertRaises(FileExistsError):
            discovery_history._json_new(
                path, {"replacement": True}
            )
        self.assertEqual(path.read_bytes(), original)

    def test_future_discovery_date_is_rejected_before_network_or_history(self):
        fixture = discovery_tests.TestPublishedDiscoveryIntegration(
            "test_pipeline_excludes_seen_published_candidates_on_second_run"
        )
        self.addCleanup(fixture.doCleanups)
        fixture.setUp()
        future = (
            datetime.now(timezone.utc).date()
            + timedelta(days=366)
        )
        payload = json.loads(crossref_body())
        payload["message"]["published-online"] = {
            "date-parts": [[
                future.year, future.month, future.day
            ]],
        }
        destination = fixture.root / "future_discovery"

        with patch.object(
            published_discovery,
            "fetch_dataset_pages",
            side_effect=fixture.fetcher([
                discovery_tests.record()
            ]),
        ), patch.object(
            published_discovery,
            "PublicationLookup",
            side_effect=fixture.factory({
                discovery_tests.DOI: payload
            }),
        ):
            with self.assertRaises(
                ValueError,
                msg="Future publication qualification date was accepted",
            ):
                published_discovery.run_published_discovery(
                    ["XL-MS"],
                    destination,
                    history_dir=fixture.history,
                    as_of=future.isoformat(),
                )

        self.assertFalse(destination.exists())
        self.assertEqual(fixture.http_calls, [])
        self.assertEqual(fixture.search_calls, [])
        self.assertEqual(fixture.seen(), set())


if __name__ == "__main__":
    unittest.main()
