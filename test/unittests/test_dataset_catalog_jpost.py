import csv
import hashlib
import json
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import requests

from claudio.dataset_catalog.repositories import RepositoryResponse, write_repository_outputs


def metadata():
    return {
        "title": "Cross-linking dataset",
        "accession": [
            {"accession": "MS:1001919", "value": "PXD000001"},
            {"accession": "MS:1002632", "value": "JPST000001"},
        ],
        "datasetFiles": [
            {
                "accession": "PRIDE:0000404", "name": "Associated raw file URI",
                "value": "ftp://ftp.jpostdb.org/JPST000001/raw/sample.raw",
            },
            {
                "accession": "PRIDE:0000408", "name": "Search engine output file URI",
                "value": "ftp://ftp.jpostdb.org/JPST000001/results/sample%20%26%20%231.tsv",
            },
        ],
    }


def response(data):
    result = requests.Response()
    result.status_code = 200
    result.url = (
        "https://repository.jpostdb.org/proxi/datasets"
        "?accession=PXD000001&resultType=full"
    )
    result.headers["Content-Type"] = "application/json"
    result._content = (json.dumps(data) + "\n\n").encode("utf-8")
    return result


def snapshot(data):
    result = response(data)
    return RepositoryResponse(
        result.url, "2026-10-31T15:00:00+00:00", result.content,
        result.headers.copy(), data,
    )


class TestJpostInventory(unittest.TestCase):
    def candidate(self):
        return {
            "accession": "PXD000001", "repository": "jPOST", "title": "Test dataset",
            "review_status": "pending", "search_terms": ["cross-linking"],
        }


    def test_writes_reported_file_links_and_keeps_unknown_counts_and_sizes_empty(self):
        data = metadata()
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            files = json.loads((output_dir / "files.json").read_bytes())

            self.assertEqual(manifest["source"], "jPOST")
            self.assertEqual(manifest["repository_accession"], "JPST000001")
            self.assertEqual(manifest["file_count"], 2)
            self.assertIsNone(manifest["reported_file_count"])
            self.assertEqual(manifest["file_inventory_scope"], "repository_reported_links")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_responses"], [])
            source = manifest["project_response"]
            self.assertEqual((output_dir / source["path"]).read_bytes(), snapshot(data).raw_body)
            self.assertEqual(source["sha256"], hashlib.sha256(snapshot(data).raw_body).hexdigest())
            self.assertEqual(files[1]["fileName"], "sample & #1.tsv")
            self.assertEqual(files[1]["repositoryMetadata"], data["datasetFiles"][1])
            self.assertEqual(files[1]["publicFileLocations"][0]["value"], data["datasetFiles"][1]["value"])
            self.assertIsNone(files[1]["fileSizeBytes"])
            self.assertNotEqual(files[0]["accession"], files[1]["accession"])
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["file_size_bytes"], "")
            self.assertEqual(rows[0]["checksum"], "")

    def test_accepts_an_explicitly_empty_link_list_but_rejects_a_missing_one(self):
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "empty"
            data = metadata()
            data["datasetFiles"] = []
            write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_count"], 0)
            self.assertIsNone(manifest["reported_file_count"])
            data.pop("datasetFiles")
            missing = Path(temp_dir) / "missing"
            with self.assertRaisesRegex(ValueError, "datasetFiles"):
                write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], missing)
            self.assertFalse(missing.exists())

    def test_rejects_unusable_file_links_before_creating_outputs(self):
        invalid = [
            "ftp://ftp.jpostdb.org/JPST000002/sample.raw",
            "ftp://example.org/JPST000001/sample.raw",
            "ftp://ftp.jpostdb.org/JPST000001/%2e%2e/sample.raw",
            "ftp://ftp.jpostdb.org/JPST000001/",
            None,
        ]
        for url in invalid:
            with self.subTest(url=url), TemporaryDirectory() as temp_dir:
                data = metadata()
                data["datasetFiles"][0]["value"] = url
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], output_dir)
                self.assertFalse(output_dir.exists())

    def test_rejects_conflicting_duplicate_metadata_and_edited_sources(self):
        for problem in ("conflicting_duplicate", "edited"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                data = metadata()
                if problem == "conflicting_duplicate":
                    duplicate = deepcopy(data["datasetFiles"][0])
                    duplicate["name"] = "Different file category"
                    data["datasetFiles"].append(duplicate)
                project = snapshot(data)
                if problem == "edited":
                    project.data["datasetFiles"][0]["value"] += ".changed"
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("jPOST", "PXD000001", project, [], output_dir)
                self.assertFalse(output_dir.exists())


class TestJpostRepeatedLinks(unittest.TestCase):
    def test_groups_identical_links_without_discarding_the_reported_occurrences(self):
        data = metadata()
        data["datasetFiles"].extend([
            deepcopy(data["datasetFiles"][0]), deepcopy(data["datasetFiles"][0]),
        ])
        with TemporaryDirectory() as temporary:
            output = Path(temporary) / "project"
            write_repository_outputs("jPOST", "PXD000001", snapshot(data), [], output)
            files = json.loads((output / "files.json").read_bytes())
            self.assertEqual(len(files), 2)
            self.assertEqual(files[0]["reportedLinkCount"], 3)
            self.assertEqual(files[0]["reportedLinkIndices"], [0, 2, 3])
            self.assertEqual(files[0]["repositoryMetadata"], data["datasetFiles"][0])
            self.assertEqual(json.loads((output / "raw/project.json").read_bytes()), data)
            manifest = json.loads((output / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_count"], 2)
            self.assertEqual(manifest["reported_link_count"], 4)
            self.assertEqual(manifest["duplicate_link_count"], 2)
            self.assertIsNone(manifest["reported_file_count"])


if __name__ == "__main__":
    unittest.main()
