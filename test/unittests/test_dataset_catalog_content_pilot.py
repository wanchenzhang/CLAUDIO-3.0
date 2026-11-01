import csv
import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from click.testing import CliRunner

from claudio.dataset_catalog.cli import main
from claudio.dataset_catalog.inventory import inspect_content_plan, prepare_content_retry
from claudio.dataset_catalog.repositories import download_pilot_file, probe_pilot_content


TABLE = b"Protein1,Protein2,Peptide1,Peptide2,Score\nP1,P2,AK,BK,10\n"


def file_item(size=None):
    return {
        "pilot_file_id": "pilot_001", "project_accession": "PXD000001", "repository": "PRIDE",
        "reported_file_name": "crosslinks.csv", "reported_size_bytes": size,
        "pilot_role": "result_table_candidate", "expected_download_file_sha256": None,
        "download_urls": ["https://example.org/crosslinks.csv"],
        "local_relative_path": "files/PXD000001/pilot_001.csv",
    }


def save_payload(url, folder):
    (Path(folder) / "payload.csv").write_bytes(TABLE)


class TestContentPlanInspection(unittest.TestCase):
    def setUp(self):
        patcher = patch("claudio.dataset_catalog.repositories._pridepy_version", return_value="0.0.16")
        patcher.start()
        self.addCleanup(patcher.stop)

    def make_plan(self, root, count=2, total=10000, unknown_size=False):
        inventory = root / "inventory"
        inventory.mkdir()
        records = [{
            "accession": f"file-{number}", "projectAccessions": ["PXD000001"],
            "fileName": f"crosslinks_{number}.csv", "fileCategory": {"value": "RESULT"},
            "fileSizeBytes": None if unknown_size else len(TABLE),
            "publicFileLocations": [{"value": f"https://example.org/crosslinks_{number}.csv"}],
        } for number in range(1, count + 1)]
        body = json.dumps(records).encode()
        (inventory / "files.json").write_bytes(body)
        digest = hashlib.sha256(body).hexdigest()
        rows, items = [], []
        for number, record in enumerate(records, 1):
            row = {
                "project_accession": "PXD000001", "repository": "PRIDE",
                "file_accession": record["accession"], "file_name": record["fileName"],
                "file_size_bytes": "" if unknown_size else str(len(TABLE)), "download_priority": "first",
                "download_blockers": "[]", "download_urls": json.dumps([record["publicFileLocations"][0]["value"]]),
                "source_inventory_path": str(inventory), "source_file_record_index": str(number - 1),
                "source_files_sha256": digest,
            }
            rows.append(row)
            item = file_item(None if unknown_size else len(TABLE))
            item.update({
                "pilot_file_id": f"pilot_{number:03d}", "file_accession": record["accession"],
                "reported_file_name": record["fileName"],
                "download_urls": [record["publicFileLocations"][0]["value"]],
                "local_relative_path": f"files/PXD000001/pilot_{number:03d}.csv",
                "source_candidate_record_number": number, "source_candidate_record": row,
                "source_candidate_record_sha256": hashlib.sha256(json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest(),
                "source_inventory_path": str(inventory), "source_file_record_index": number - 1,
                "inventory_files_json_sha256": digest,
            })
            items.append(item)
        catalog = root / "candidates.csv"
        with catalog.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        plan = {
            "schema_version": 1, "kind": "xlms_file_content_pilot", "counts": {"files": count},
            "source": {"server_candidate_csv": str(catalog), "candidate_record_count": count},
            "proposed_download_limits": {"per_file_bytes": 1024, "total_payload_bytes": total},
            "files": items,
        }
        path = root / "plan.json"
        path.write_text(json.dumps(plan), encoding="utf-8")
        return path

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_saves_each_outcome_and_reuses_verified_completed_downloads(self, mock_get):
        mock_get.side_effect = save_payload
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root)
            output = root / "run"
            results = inspect_content_plan(plan, output)
            self.assertEqual([row["download_status"] for row in results], ["downloaded", "downloaded"])
            self.assertTrue((output / "manifest.json").is_file())
            self.assertTrue((output / "results.csv").is_file())
            mock_get.reset_mock()
            self.assertEqual(inspect_content_plan(plan, output, resume=True), results)
            mock_get.assert_not_called()
            Path(results[0]["local_path"]).write_bytes(b"changed")
            before = (output / "content_checkpoint.json").read_bytes()
            with self.assertRaises(ValueError):
                inspect_content_plan(plan, output, resume=True)
            self.assertEqual((output / "content_checkpoint.json").read_bytes(), before)

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_source_mismatch_is_rejected_before_creating_output_or_requesting_files(self, mock_get):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan_path = self.make_plan(root)
            plan = json.loads(plan_path.read_bytes())
            plan["files"][0]["download_urls"] = ["https://example.org/another.csv"]
            plan_path.write_text(json.dumps(plan), encoding="utf-8")
            with self.assertRaises(ValueError):
                inspect_content_plan(plan_path, root / "run")
            self.assertFalse((root / "run").exists())
            mock_get.assert_not_called()

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_resume_after_a_saved_result_preserves_budget_and_does_not_refetch_it(self, mock_get):
        mock_get.side_effect = save_payload
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root)
            output = root / "run"
            def interrupt(completed, total, result):
                if completed == 1:
                    raise KeyboardInterrupt()
            with self.assertRaises(KeyboardInterrupt):
                inspect_content_plan(plan, output, on_progress=interrupt)
            mock_get.reset_mock()
            results = inspect_content_plan(plan, output, resume=True)
            self.assertEqual(len(results), 2)
            self.assertEqual(mock_get.call_count, 1)
            checkpoint = json.loads((output / "content_checkpoint.json").read_bytes())
            self.assertEqual(checkpoint["charged_payload_bytes"], 2 * len(TABLE))

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_interrupted_active_transfer_keeps_its_reservation_when_retried(self, mock_get):
        mock_get.side_effect = KeyboardInterrupt()
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root, count=1)
            output = root / "run"
            with self.assertRaises(KeyboardInterrupt):
                inspect_content_plan(plan, output)
            checkpoint = json.loads((output / "content_checkpoint.json").read_bytes())
            self.assertEqual(checkpoint["charged_payload_bytes"], len(TABLE))
            self.assertIsNotNone(checkpoint["active"])
            mock_get.reset_mock()
            mock_get.side_effect = save_payload
            result = inspect_content_plan(plan, output, resume=True)
            self.assertEqual(result[0]["download_status"], "downloaded")
            checkpoint = json.loads((output / "content_checkpoint.json").read_bytes())
            self.assertEqual(len(checkpoint["interrupted_attempts"]), 1)
            self.assertEqual(checkpoint["charged_payload_bytes"], 2 * len(TABLE))
            self.assertEqual(mock_get.call_count, 1)

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_budget_is_enforced_and_partial_files_are_not_accepted(self, mock_get):
        mock_get.side_effect = save_payload
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root, total=len(TABLE) + 1)
            output = root / "run"
            results = inspect_content_plan(plan, output)
            self.assertEqual(results[0]["download_status"], "downloaded")
            self.assertEqual(results[1]["download_status"], "budget_exhausted")
            self.assertLessEqual(sum(row["charged_payload_bytes"] for row in results), len(TABLE) + 1)
            self.assertEqual(len(list((output / "files").rglob("*.csv"))), 1)

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_retry_plan_keeps_only_failures_and_carries_remaining_budget(self, mock_get):
        outcomes = [RuntimeError("server unavailable"), TABLE]
        def transfer(url, folder):
            body = outcomes.pop(0)
            if isinstance(body, Exception):
                raise body
            (Path(folder) / "payload.csv").write_bytes(body)
        mock_get.side_effect = transfer
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan_path = self.make_plan(root)
            parent = root / "run"
            first = inspect_content_plan(plan_path, parent)
            self.assertEqual([row["download_status"] for row in first], ["failed", "downloaded"])
            before = {str(path.relative_to(parent)): path.read_bytes() for path in parent.rglob("*") if path.is_file()}
            mock_get.reset_mock()
            retry_path = root / "retry_plan.json"
            prepare_content_retry(parent, retry_path)
            mock_get.assert_not_called()
            retry = json.loads(retry_path.read_bytes())
            original = json.loads(plan_path.read_bytes())
            self.assertEqual(retry["counts"]["files"], 1)
            self.assertEqual(retry["files"], original["files"][:1])
            self.assertEqual(retry["proposed_download_limits"]["total_payload_bytes"], 10000 - 2 * len(TABLE))
            self.assertEqual(retry["retry_context"]["parent_charged_payload_bytes"], 2 * len(TABLE))
            self.assertEqual(before, {str(path.relative_to(parent)): path.read_bytes()
                                      for path in parent.rglob("*") if path.is_file()})
            mock_get.side_effect = save_payload
            results = inspect_content_plan(retry_path, root / "retry_run")
            self.assertEqual(results[0]["download_status"], "downloaded", results)

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_retry_plan_is_not_written_when_the_original_budget_is_exhausted(self, mock_get):
        mock_get.side_effect = RuntimeError("server unavailable")
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root, count=1, total=len(TABLE))
            parent = root / "run"
            inspect_content_plan(plan, parent)
            with self.assertRaisesRegex(ValueError, "budget"):
                prepare_content_retry(parent, root / "retry_plan.json")
            self.assertFalse((root / "retry_plan.json").exists())

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_inspect_accepts_a_plan_and_rejects_combined_input_modes(self, mock_get):
        mock_get.side_effect = save_payload
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root, count=1)
            result = CliRunner().invoke(main, ["inspect", "--plan", str(plan), "--output", str(root / "run"), "--delay", "0"])
            self.assertEqual(result.exit_code, 0, result.output)
            self.assertIn("Downloaded: 1", result.output)
            result = CliRunner().invoke(main, ["inspect", "--plan", str(plan), "--catalog", str(plan), "--output", str(root / "other")])
            self.assertEqual(result.exit_code, 2, result.output)
            self.assertFalse((root / "other").exists())


    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_unknown_sizes_are_saved_as_skipped_and_resume_without_requests(self, download):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root, count=1, unknown_size=True)
            results = inspect_content_plan(plan, root / "run")
            self.assertEqual(results[0]["download_status"], "skipped_unknown_size")
            self.assertEqual(results[0]["charged_payload_bytes"], 0)
            self.assertEqual(inspect_content_plan(plan, root / "run", resume=True), results)
            manifest = json.loads((root / "run/manifest.json").read_bytes())
            self.assertEqual(manifest["accounting_mode"], "planned_file_sizes")
        download.assert_not_called()

    @patch("claudio.dataset_catalog.repositories._pridepy_download")
    def test_completed_legacy_checkpoint_is_reused_without_changing_reports(self, download):
        download.side_effect = save_payload
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            plan = self.make_plan(root, count=1)
            output = root / "run"
            results = inspect_content_plan(plan, output)
            checkpoint_path = output / "content_checkpoint.json"
            checkpoint = json.loads(checkpoint_path.read_bytes())
            checkpoint["schema_version"] = 1
            checkpoint.pop("accounting_mode")
            checkpoint_path.write_text(json.dumps(checkpoint))
            manifest_path = output / "manifest.json"
            manifest = json.loads(manifest_path.read_bytes())
            manifest["output_sha256"]["content_checkpoint.json"] = hashlib.sha256(checkpoint_path.read_bytes()).hexdigest()
            manifest_path.write_text(json.dumps(manifest))
            before = {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()}
            download.reset_mock()
            self.assertEqual(inspect_content_plan(plan, output, resume=True), results)
            self.assertEqual(before, {p.relative_to(output): p.read_bytes() for p in output.rglob("*") if p.is_file()})
            download.assert_not_called()


class TestPilotContentProbe(unittest.TestCase):
    def test_checks_fasta_structure_without_claiming_search_database_identity(self):
        result = probe_pilot_content(b">P1 description\nMKKAA\n>P2\nAKK\n", "fasta_candidate")
        self.assertEqual(result["status"], "fasta_candidate")
        self.assertEqual(result["sequence_count"], 2)
        self.assertEqual(result["example_ids"], ["P1", "P2"])
        self.assertEqual(result["search_database_association"], "unconfirmed")

    def test_handles_a_title_before_the_table_header_and_header_only_files(self):
        data = b"Filtered cross-linked sites\n" + TABLE
        result = probe_pilot_content(data, "result_table_candidate")
        self.assertEqual(result["status"], "table_candidate")
        self.assertEqual(result["header_record_index"], 1)
        self.assertEqual(result["rows_after_header"], 1)
        empty = probe_pilot_content(TABLE.splitlines()[0] + b"\n", "result_table_candidate")
        self.assertEqual(empty["status"], "header_only_table")

    def test_keeps_controls_and_html_separate_from_result_candidates(self):
        control = probe_pilot_content(b"python convert.py input.txt output.tsv\n", "possible_instructions_control")
        self.assertEqual(control["status"], "text_control")
        library = probe_pilot_content(TABLE, "possible_library_control")
        self.assertEqual(library["status"], "table_control")
        html = probe_pilot_content(b"<!doctype html><html>Login</html>", "result_table_candidate")
        self.assertEqual(html["status"], "html_response")
        for result in (control, library, html):
            self.assertNotIn("fdr_passed", result)
            self.assertNotIn("eligible", result)


if __name__ == "__main__":
    unittest.main()
