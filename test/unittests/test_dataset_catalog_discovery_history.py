
"""Offline history contracts: complete outputs, incremental selection, and recovery."""

import copy
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib
import json
import os
from pathlib import Path
import sqlite3
import tempfile
from threading import Barrier
import unittest
import contextlib
from unittest.mock import patch


def candidate(accession="PXD000001", status="published", doi="10.1234/shared"):
    return {
        "accession": accession, "title": "An XL-MS dataset", "search_terms": ["XL-MS"],
        "publication_status": status,
        "papers": ([{"doi": doi, "pmid": None, "title": "An article", "venue": "A journal",
                    "year": 2025, "authors": []}] if status == "published" else []),
    }


def write_outputs(destination, candidates):
    (destination / "candidates.json").write_text(
        json.dumps(candidates, allow_nan=False) + "\n", encoding="utf-8")
    (destination / "audit.json").write_text('{"source":"offline fixture"}\n', encoding="utf-8")
    return {"source": "offline_fixture"}


class TestDiscoveryHistory(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.directory = self.root / "history"
        self.run_count = 0

    def api(self):
        try:
            return importlib.import_module("claudio.dataset_catalog.discovery_history")
        except ModuleNotFoundError as error:
            if error.name != "claudio.dataset_catalog.discovery_history":
                raise
            self.fail("Persistent discovery history is not implemented")

    def history(self):
        api = self.api()
        if not self.directory.exists():
            api.initialize_discovery_history(self.directory)
        history = api.DiscoveryHistory(self.directory)
        self.addCleanup(history.close)
        return history

    def publish(self, history, candidates, writer=write_outputs):
        self.run_count += 1
        run = history.prepare_run(self.root / ("run_" + str(self.run_count)))
        return history.publish(run["run_id"], candidates, writer)

    def test_successful_published_output_is_excluded_on_the_next_run(self):
        history = self.history()
        first = self.publish(history, [candidate()])
        second = self.publish(history, [candidate()])
        self.assertEqual(first["emitted_accessions"], ["PXD000001"])
        self.assertEqual(second["emitted_accessions"], [])
        self.assertEqual(second["candidate_count"], 0)
        self.assertEqual(history.seen_accessions(), {"PXD000001"})
        saved = json.loads((Path(second["output_directory"]) / "candidates.json").read_bytes())
        self.assertEqual(saved, [])

    def test_missing_history_requires_explicit_initialization(self):
        with self.assertRaisesRegex(ValueError, "(?i)initializ|missing"):
            self.api().DiscoveryHistory(self.directory)
        self.assertFalse(self.directory.exists())

    def test_initialization_cannot_reset_an_existing_history(self):
        history = self.history()
        self.publish(history, [candidate()])
        history.close()
        before = {p.name: p.read_bytes() for p in self.directory.iterdir() if p.is_file()}
        with self.assertRaises(FileExistsError):
            self.api().initialize_discovery_history(self.directory)
        after = {p.name: p.read_bytes() for p in self.directory.iterdir() if p.is_file()}
        self.assertEqual(before, after)

    def test_deleted_database_is_not_silently_recreated(self):
        history = self.history()
        history.close()
        database = self.directory / "history.sqlite3"
        database.unlink()
        with self.assertRaisesRegex(ValueError, "(?i)missing|initializ"):
            self.api().DiscoveryHistory(self.directory)
        self.assertFalse(database.exists())

    def test_corrupt_database_is_rejected(self):
        history = self.history()
        history.close()
        (self.directory / "history.sqlite3").write_bytes(b"not a sqlite database")
        with self.assertRaises(ValueError):
            self.api().DiscoveryHistory(self.directory)

    def test_namespace_marker_must_match_the_database(self):
        history = self.history()
        history.close()
        marker = self.directory / "namespace.json"
        data = json.loads(marker.read_bytes())
        data["namespace_id"] = "00000000-0000-4000-8000-000000000000"
        marker.write_text(json.dumps(data), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "(?i)namespace"):
            self.api().DiscoveryHistory(self.directory)

    def test_unpublished_candidate_can_enter_after_later_publication(self):
        history = self.history()
        first = self.publish(history, [candidate(status="pending")])
        self.assertEqual(first["candidate_count"], 0)
        self.assertEqual(history.seen_accessions(), set())
        second = self.publish(history, [candidate()])
        self.assertEqual(second["emitted_accessions"], ["PXD000001"])

    def test_lookup_failures_and_ineligible_records_do_not_consume_history(self):
        history = self.history()
        records = [candidate("PXD000001", "unresolved"), candidate("PXD000002", "ineligible"),
                   candidate("PXD000003", "unverified"), candidate("PXD000004", "unpublished")]
        result = self.publish(history, records)
        self.assertEqual(result["candidate_count"], 0)
        self.assertEqual(history.seen_accessions(), set())

    def test_accession_not_paper_doi_is_the_history_key(self):
        history = self.history()
        self.publish(history, [candidate("PXD000001")])
        second = self.publish(history, [candidate("PXD000002")])
        self.assertEqual(second["emitted_accessions"], ["PXD000002"])
        self.assertEqual(history.seen_accessions(), {"PXD000001", "PXD000002"})

    def test_changed_search_terms_and_new_papers_do_not_reemit_seen_accession(self):
        history = self.history()
        self.publish(history, [candidate()])
        changed = candidate(doi="10.1234/new")
        changed["search_terms"] = ["crosslinking"]
        changed["title"] = "Updated title"
        self.assertEqual(self.publish(history, [changed])["candidate_count"], 0)

    def test_accession_keys_are_normalized_across_runs(self):
        history = self.history()
        first = self.publish(history, [candidate(" pxd000001 ")])
        second = self.publish(history, [candidate("PXD000001")])
        self.assertEqual(first["emitted_accessions"], ["PXD000001"])
        self.assertEqual(second["candidate_count"], 0)

    def test_prepare_does_not_mark_any_candidates_as_seen(self):
        history = self.history()
        run = history.prepare_run(self.root / "prepared")
        self.assertTrue(run["run_id"])
        self.assertFalse((self.root / "prepared").exists())
        self.assertEqual(history.seen_accessions(), set())

    def test_failure_before_complete_manifest_does_not_consume_candidates(self):
        history = self.history()

        def interrupted(destination, selected):
            (destination / "partial.txt").write_bytes(b"partial")
            raise OSError("Interrupted output")

        with self.assertRaises(OSError):
            self.publish(history, [candidate()], interrupted)
        self.assertFalse((self.root / "run_1/manifest.json").exists())
        self.assertEqual(history.seen_accessions(), set())
        self.assertEqual(self.publish(history, [candidate()])["candidate_count"], 1)
        self.assertEqual((self.root / "run_1/partial.txt").read_bytes(), b"partial")

    def test_complete_output_is_recovered_after_database_registration_failure(self):
        history = self.history()
        run = history.prepare_run(self.root / "crashed")
        with patch.object(history, "_record_completed", side_effect=RuntimeError("Simulated crash")):
            with self.assertRaises(RuntimeError):
                history.publish(run["run_id"], [candidate()], write_outputs)
        self.assertTrue((self.root / "crashed/manifest.json").is_file())
        history.close()
        reopened = self.api().DiscoveryHistory(self.directory)
        self.addCleanup(reopened.close)
        recovered = reopened.recover()
        self.assertEqual(recovered["recovered_run_ids"], [run["run_id"]])
        self.assertEqual(reopened.seen_accessions(), {"PXD000001"})
        self.assertEqual(self.publish(reopened, [candidate()])["candidate_count"], 0)

    def test_completed_run_replay_does_not_invoke_the_writer_again(self):
        history = self.history()
        run = history.prepare_run(self.root / "completed")
        first = history.publish(run["run_id"], [candidate()], write_outputs)

        def must_not_write(*args):
            raise AssertionError("Completed output must be immutable")

        second = history.publish(run["run_id"], [candidate()], must_not_write)
        self.assertEqual(second, first)

    def test_partial_output_is_preserved_and_cannot_be_overwritten(self):
        history = self.history()
        run = history.prepare_run(self.root / "partial")
        (self.root / "partial").mkdir()
        marker = self.root / "partial/keep.txt"
        marker.write_bytes(b"keep")
        with self.assertRaises(FileExistsError):
            history.publish(run["run_id"], [candidate()], write_outputs)
        self.assertEqual(marker.read_bytes(), b"keep")
        self.assertEqual(history.seen_accessions(), set())

    def test_existing_output_directory_is_rejected_before_preparation(self):
        history = self.history()
        target = self.root / "existing"
        target.mkdir()
        (target / "keep").write_bytes(b"keep")
        with self.assertRaises(FileExistsError):
            history.prepare_run(target)
        self.assertEqual((target / "keep").read_bytes(), b"keep")

    def test_invalid_candidate_metadata_does_not_create_a_complete_output(self):
        history = self.history()
        invalid = candidate()
        invalid["papers"] = []
        with self.assertRaises(ValueError):
            self.publish(history, [invalid])
        self.assertEqual(history.seen_accessions(), set())
        self.assertFalse((self.root / "run_1").exists())

    def test_duplicate_accessions_are_rejected_even_with_different_case(self):
        history = self.history()
        with self.assertRaisesRegex(ValueError, "(?i)duplicate"):
            self.publish(history, [candidate("PXD000001"), candidate("pxd000001")])
        self.assertEqual(history.seen_accessions(), set())

    def test_written_candidates_must_match_the_selected_records(self):
        history = self.history()

        def wrong_writer(destination, selected):
            write_outputs(destination, [candidate("PXD999999")])
            return {}

        with self.assertRaisesRegex(ValueError, "(?i)candidates|selected"):
            self.publish(history, [candidate()], wrong_writer)
        self.assertFalse((self.root / "run_1/manifest.json").exists())
        self.assertEqual(history.seen_accessions(), set())

    def test_writer_cannot_override_history_manifest_fields(self):
        history = self.history()

        def wrong_writer(destination, selected):
            write_outputs(destination, selected)
            return {"namespace_id": "foreign"}

        with self.assertRaisesRegex(ValueError, "(?i)reserved|manifest"):
            self.publish(history, [candidate()], wrong_writer)
        self.assertEqual(history.seen_accessions(), set())

    def test_complete_receipt_pins_all_output_files(self):
        history = self.history()
        result = self.publish(history, [candidate()])
        directory = Path(result["output_directory"])
        raw = (directory / "manifest.json").read_bytes()
        self.assertEqual(hashlib.sha256(raw).hexdigest(), result["manifest_sha256"])
        manifest = json.loads(raw)
        self.assertEqual(set(manifest["output_sha256"]), {"candidates.json", "audit.json"})
        self.assertEqual(manifest["emitted_accessions"], ["PXD000001"])
        for name, digest in manifest["output_sha256"].items():
            self.assertEqual(hashlib.sha256((directory / name).read_bytes()).hexdigest(), digest)

    def test_changed_completed_output_causes_a_visible_error(self):
        history = self.history()
        result = self.publish(history, [candidate()])
        (Path(result["output_directory"]) / "audit.json").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            history.seen_accessions()

    def test_changed_completed_manifest_causes_a_visible_error(self):
        history = self.history()
        result = self.publish(history, [candidate()])
        manifest = Path(result["output_directory"]) / "manifest.json"
        manifest.write_bytes(manifest.read_bytes() + b"\n")
        with self.assertRaisesRegex(ValueError, "(?i)checksum"):
            history.recover()

    def test_deleted_seen_entry_is_detected_instead_of_reemitting(self):
        history = self.history()
        self.publish(history, [candidate()])
        with contextlib.closing(sqlite3.connect(self.directory / "history.sqlite3")) as connection:
            connection.execute("DELETE FROM seen_accessions")
            connection.commit()
        with self.assertRaisesRegex(ValueError, "(?i)history|accession"):
            history.seen_accessions()

    def test_foreign_namespace_complete_manifest_is_not_recovered(self):
        history = self.history()
        run = history.prepare_run(self.root / "foreign")
        with patch.object(history, "_record_completed", side_effect=RuntimeError("Simulated crash")):
            with self.assertRaises(RuntimeError):
                history.publish(run["run_id"], [candidate()], write_outputs)
        path = self.root / "foreign/manifest.json"
        manifest = json.loads(path.read_bytes())
        manifest["namespace_id"] = "another namespace"
        path.write_text(json.dumps(manifest), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "(?i)namespace|receipt"):
            history.recover()

    def test_corrupt_prepared_complete_manifest_is_not_treated_as_partial(self):
        history = self.history()
        history.prepare_run(self.root / "invalid")
        (self.root / "invalid").mkdir()
        (self.root / "invalid/manifest.json").write_bytes(b"{invalid")
        with self.assertRaises(ValueError):
            history.recover()

    def test_two_concurrent_runs_emit_the_same_accession_only_once(self):
        api = self.api()
        api.initialize_discovery_history(self.directory)
        barrier = Barrier(2)

        def run(number):
            with api.DiscoveryHistory(self.directory) as history:
                prepared = history.prepare_run(self.root / ("concurrent_" + str(number)))
                barrier.wait(timeout=10)
                return history.publish(prepared["run_id"], [candidate()], write_outputs)

        with ThreadPoolExecutor(max_workers=2) as executor:
            results = list(executor.map(run, (1, 2)))
        self.assertEqual(sorted(result["candidate_count"] for result in results), [0, 1])
        self.assertEqual(self.history().seen_accessions(), {"PXD000001"})

    def test_candidates_are_not_mutated_and_new_record_order_is_preserved(self):
        history = self.history()
        original = [candidate("PXD000003"), candidate("PXD000001"), candidate("PXD000002")]
        before = copy.deepcopy(original)
        result = self.publish(history, original)
        self.assertEqual(original, before)
        self.assertEqual(result["emitted_accessions"], ["PXD000003", "PXD000001", "PXD000002"])

    def test_default_history_location_respects_xdg_state_home(self):
        with patch.dict(os.environ, {"XDG_STATE_HOME": str(self.root / "state")}):
            self.assertEqual(self.api().default_history_directory(), self.root / "state/claudio/discovery")
        self.assertFalse((self.root / "state").exists())


if __name__ == "__main__":
    unittest.main()
