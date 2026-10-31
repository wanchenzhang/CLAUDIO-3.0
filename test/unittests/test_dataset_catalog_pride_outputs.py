import csv
import hashlib
import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from claudio.dataset_catalog.pride import PrideResponse
from claudio.dataset_catalog.pride_outputs import write_pride_outputs


def make_snapshot(url, data, headers=None):
    return PrideResponse(
        request_url=url,
        retrieved_at_utc="2026-10-31T16:00:00+00:00",
        raw_body=(
            json.dumps(data, separators=(",", ":")) + "\n\n"
        ).encode("utf-8"),
        headers={"Content-Type": "application/json", **(headers or {})},
        data=data,
    )


class TestDatasetCatalogPrideOutputs(unittest.TestCase):
    def test_writes_project_files_and_verifiable_response_sources(self):
        base_url = (
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD035844"
        )
        project = make_snapshot(base_url, {
            "accession": "PXD035844",
            "title": "Test cross-linking dataset",
            "references": [{"doi": "10.1073/pnas.2219418120"}],
        })
        records = [
            {
                "accession": "file-1",
                "projectAccessions": ["PXD035844"],
                "fileName": "sample.raw",
                "fileCategory": {"value": "RAW"},
                "fileSizeBytes": 100,
                "checksum": "a" * 40,
                "publicFileLocations": [
                    {"value": "ftp://example.org/sample.raw"}
                ],
            },
            {
                "accession": "file-2",
                "projectAccessions": ["PXD035844", "PXD013928"],
                "fileName": "results.zip",
                "fileCategory": {"value": "SEARCH"},
                "fileSizeBytes": 200,
                "checksum": "b" * 40,
                "publicFileLocations": [
                    {"value": "ftp://example.org/results.zip"}
                ],
            },
        ]
        file_pages = [
            make_snapshot(
                f"{base_url}/files?pageSize=1&page=0",
                [records[0]], {"total_records": "2"},
            ),
            make_snapshot(
                f"{base_url}/files?pageSize=1&page=1",
                [records[1]], {"total_records": "2"},
            ),
        ]

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "PXD035844"
            write_pride_outputs(
                "PXD035844", project, file_pages, output_dir
            )

            self.assertTrue((output_dir / "manifest.json").is_file())
            self.assertEqual(
                json.loads(
                    (output_dir / "project.json").read_text(encoding="utf-8")
                ),
                project.data,
            )
            self.assertEqual(
                json.loads(
                    (output_dir / "files.json").read_text(encoding="utf-8")
                ),
                records,
            )

            with (output_dir / "files.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(
                [row["project_accession"] for row in rows],
                ["PXD035844", "PXD035844"],
            )
            self.assertEqual(
                [row["file_accession"] for row in rows],
                ["file-1", "file-2"],
            )
            self.assertEqual(rows[0]["file_name"], "sample.raw")
            self.assertEqual(rows[0]["file_size_bytes"], "100")
            self.assertEqual(rows[0]["file_category"], "RAW")
            self.assertEqual(rows[0]["checksum"], "a" * 40)
            self.assertEqual(
                json.loads(rows[0]["public_file_locations"]),
                records[0]["publicFileLocations"],
            )
            self.assertEqual(
                json.loads(rows[1]["project_accessions"]),
                ["PXD035844", "PXD013928"],
            )

            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["project_accession"], "PXD035844")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_count"], 2)
            self.assertEqual(manifest["reported_file_count"], 2)

            saved_sources = [
                manifest["project_response"], *manifest["file_responses"]
            ]
            self.assertEqual(len(saved_sources), 3)
            for saved, snapshot in zip(
                saved_sources, [project, *file_pages]
            ):
                relative_path = Path(saved["path"])
                self.assertFalse(relative_path.is_absolute())
                self.assertEqual(
                    (output_dir / relative_path).read_bytes(), snapshot.raw_body
                )
                self.assertEqual(
                    saved["sha256"],
                    hashlib.sha256(snapshot.raw_body).hexdigest(),
                )
                self.assertEqual(saved["request_url"], snapshot.request_url)
                self.assertEqual(
                    saved["retrieved_at_utc"], snapshot.retrieved_at_utc
                )
                self.assertEqual(saved["headers"], dict(snapshot.headers))
    def _inputs(self):
        base_url = (
            "https://www.ebi.ac.uk/pride/ws/archive/v3/projects/PXD035844"
        )
        project = make_snapshot(base_url, {
            "accession": "PXD035844", "title": "Test dataset",
        })
        page = make_snapshot(
            f"{base_url}/files?pageSize=100&page=0",
            [{
                "accession": "file-1",
                "projectAccessions": ["PXD035844"],
                "fileName": "sample.raw",
            }],
            {"total_records": "1"},
        )
        return project, [page]

    def test_writes_empty_inventory_with_a_pending_screening_status(self):
        project, pages = self._inputs()
        empty_page = make_snapshot(
            pages[0].request_url, [], {"total_records": "0"}
        )

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            write_pride_outputs(
                "PXD035844", project, [empty_page], output_dir
            )

            self.assertEqual(
                json.loads((output_dir / "files.json").read_text(encoding="utf-8")),
                [],
            )
            with (output_dir / "files.csv").open(
                encoding="utf-8", newline=""
            ) as handle:
                reader = csv.DictReader(handle)
                self.assertIn("file_accession", reader.fieldnames)
                self.assertEqual(list(reader), [])

            manifest = json.loads(
                (output_dir / "manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_count"], 0)
            self.assertEqual(len(manifest["file_responses"]), 1)

    def test_missing_file_responses_are_not_treated_as_zero_files(self):
        project, _ = self._inputs()
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"

            with self.assertRaises(ValueError):
                write_pride_outputs("PXD035844", project, [], output_dir)

            self.assertFalse(output_dir.exists())

    def test_rejects_project_or_file_identity_mismatches_before_writing(self):
        for source in ("project", "file"):
            with self.subTest(source=source):
                project, pages = self._inputs()
                if source == "project":
                    project = make_snapshot(
                        project.request_url, {"accession": "PXD013928"}
                    )
                else:
                    record = dict(pages[0].data[0])
                    record["projectAccessions"] = ["PXD013928"]
                    pages = [make_snapshot(
                        pages[0].request_url, [record], {"total_records": "1"}
                    )]

                with TemporaryDirectory() as temp_dir:
                    output_dir = Path(temp_dir) / "run"
                    with self.assertRaises(ValueError):
                        write_pride_outputs(
                            "PXD035844", project, pages, output_dir
                        )
                    self.assertFalse(output_dir.exists())

    def test_rejects_incomplete_file_inventory_before_writing(self):
        project, pages = self._inputs()
        incomplete = make_snapshot(
            pages[0].request_url, pages[0].data, {"total_records": "2"}
        )

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            with self.assertRaisesRegex(ValueError, "file count mismatch"):
                write_pride_outputs(
                    "PXD035844", project, [incomplete], output_dir
                )
            self.assertFalse(output_dir.exists())

    def test_rejects_duplicate_file_records_before_writing(self):
        project, pages = self._inputs()
        repeated = make_snapshot(
            pages[0].request_url, pages[0].data * 2, {"total_records": "2"}
        )

        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            with self.assertRaisesRegex(ValueError, "Repeated PRIDE file"):
                write_pride_outputs(
                    "PXD035844", project, [repeated], output_dir
                )
            self.assertFalse(output_dir.exists())

    def test_rejects_metadata_that_differs_from_the_saved_raw_response(self):
        for source in ("project", "file"):
            with self.subTest(source=source):
                project, pages = self._inputs()
                if source == "project":
                    project.data["title"] = "Edited metadata"
                else:
                    pages[0].data[0]["fileName"] = "edited.raw"

                with TemporaryDirectory() as temp_dir:
                    output_dir = Path(temp_dir) / "run"
                    with self.assertRaisesRegex(ValueError, "raw response"):
                        write_pride_outputs(
                            "PXD035844", project, pages, output_dir
                        )
                    self.assertFalse(output_dir.exists())

    def test_refuses_to_overwrite_an_existing_inventory(self):
        project, pages = self._inputs()
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            old_csv = output_dir / "files.csv"
            old_csv.write_text("previous inventory\n", encoding="utf-8")

            with self.assertRaises(FileExistsError):
                write_pride_outputs(
                    "PXD035844", project, pages, output_dir
                )

            self.assertEqual(
                old_csv.read_text(encoding="utf-8"), "previous inventory\n"
            )
            self.assertFalse((output_dir / "manifest.json").exists())

    def test_write_failure_does_not_publish_a_complete_manifest(self):
        project, pages = self._inputs()
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "run"
            with patch(
                "claudio.dataset_catalog.pride_outputs.csv.DictWriter.writeheader",
                side_effect=OSError("No space left on device"),
            ):
                with self.assertRaises(OSError):
                    write_pride_outputs(
                        "PXD035844", project, pages, output_dir
                    )

            self.assertTrue((output_dir / "files.json").is_file())
            self.assertFalse((output_dir / "manifest.json").exists())

if __name__ == "__main__":
    unittest.main()