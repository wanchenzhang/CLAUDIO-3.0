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
        "accession": {"accession": "MS:1001919", "value": "PXD000001"},
        "title": "iProX cross-linking dataset",
        "datasetLink": {
            "accession": "MS:1002837", "name": "iProX dataset URI",
            "value": "https://www.iprox.org/page/project.html?id=IPX0000001000",
        },
        "dataFiles": [
            {
                "accession": "MS:1002846", "name": "Associated raw file URI",
                "value": "https://download.iprox.org/IPX0000001000/IPX0000001001/sample.raw",
            },
            {
                "accession": "MS:1002849", "name": "Search engine output file URI",
                "value": "http://download.iprox.org/IPX0000001000/IPX0000001002/sample%20results.tsv",
            },
        ],
    }


def response(data):
    result = requests.Response()
    result.status_code = 200
    result.url = "https://www.iprox.cn/proxi/datasets/PXD000001"
    result.headers["Content-Type"] = "application/json"
    result._content = (json.dumps(data) + "\n\n").encode("utf-8")
    return result


def snapshot(data):
    raw = response(data)
    return RepositoryResponse(
        raw.url, "2026-10-31T16:00:00+00:00", raw.content,
        raw.headers.copy(), data,
    )


class TestIproxInventory(unittest.TestCase):
    def candidate(self):
        return {
            "accession": "PXD000001", "repository": "iProX", "title": "Test dataset",
            "search_terms": ["crosslinking"], "review_status": "pending",
        }


    def test_writes_reported_links_with_subproject_paths_and_unknown_sizes(self):
        data = metadata()
        project = snapshot(data)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("iProX", "PXD000001", project, [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            files = json.loads((output_dir / "files.json").read_bytes())

            self.assertEqual(manifest["source"], "iProX")
            self.assertEqual(manifest["repository_accession"], "IPX0000001000")
            self.assertEqual(manifest["file_count"], 2)
            self.assertIsNone(manifest["reported_file_count"])
            self.assertEqual(manifest["file_inventory_scope"], "repository_reported_links")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_responses"], [])
            self.assertEqual(files[0]["relativePath"], "IPX0000001001/sample.raw")
            self.assertEqual(files[1]["fileName"], "sample results.tsv")
            self.assertEqual(files[1]["repositoryMetadata"], data["dataFiles"][1])
            self.assertEqual(files[1]["publicFileLocations"][0]["value"], data["dataFiles"][1]["value"])
            self.assertIsNone(files[0]["fileSizeBytes"])
            source = manifest["project_response"]
            self.assertEqual((output_dir / source["path"]).read_bytes(), project.raw_body)
            self.assertEqual(source["sha256"], hashlib.sha256(project.raw_body).hexdigest())
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["file_size_bytes"], "")
            self.assertEqual(rows[0]["checksum"], "")

    def test_keeps_explicitly_empty_file_lists_distinct_from_missing_metadata(self):
        with TemporaryDirectory() as temp_dir:
            data = metadata()
            data["dataFiles"] = []
            output_dir = Path(temp_dir) / "empty"
            write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_count"], 0)
            self.assertIsNone(manifest["reported_file_count"])
            data.pop("dataFiles")
            output_dir = Path(temp_dir) / "missing"
            with self.assertRaisesRegex(ValueError, "dataFiles"):
                write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
            self.assertFalse(output_dir.exists())

    def test_accepts_a_single_file_object_without_changing_the_original_metadata(self):
        data = metadata()
        data["dataFiles"] = data["dataFiles"][0]
        before = deepcopy(data)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            self.assertEqual(manifest["file_count"], 1)
            self.assertEqual(json.loads((output_dir / "project.json").read_bytes()), before)
            self.assertEqual(data, before)

    def test_rejects_wrong_roots_invalid_links_and_duplicates_before_writing(self):
        for problem in ("root", "host", "traversal", "duplicate"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                data = metadata()
                if problem == "root":
                    data["dataFiles"][0]["value"] = "https://download.iprox.org/IPX0000002000/sample.raw"
                elif problem == "host":
                    data["dataFiles"][0]["value"] = "https://example.org/IPX0000001000/sample.raw"
                elif problem == "traversal":
                    data["dataFiles"][0]["value"] = "https://download.iprox.org/IPX0000001000/%2e%2e/sample.raw"
                else:
                    data["dataFiles"].append(deepcopy(data["dataFiles"][0]))
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output_dir)
                self.assertFalse(output_dir.exists())


class TestIproxUndecodablePaths(unittest.TestCase):
    def test_preserves_an_undecodable_uri_path_without_renaming_the_file(self):
        data = metadata()
        uri = "https://download.iprox.org/IPX0000001000/IPX0000001001/results_1%FDR.plink"
        data["dataFiles"][0]["value"] = uri
        with TemporaryDirectory() as temporary:
            output = Path(temporary) / "project"
            write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output)
            files = json.loads((output / "files.json").read_bytes())
            self.assertEqual(files[0]["publicFileLocations"][0]["value"], uri)
            self.assertEqual(files[0]["fileName"], "results_1%FDR.plink")
            self.assertEqual(files[0]["relativePath"], "IPX0000001001/results_1%FDR.plink")
            self.assertEqual(files[0]["pathDecodingStatus"], "unresolved")
            self.assertEqual(files[0]["repositoryMetadata"], data["dataFiles"][0])
            self.assertEqual(json.loads((output / "raw/project.json").read_bytes()), data)
            manifest = json.loads((output / "manifest.json").read_bytes())
            self.assertEqual(manifest["unresolved_path_count"], 1)

    def test_encoded_parent_directory_stays_invalid_even_with_non_utf8_bytes(self):
        data = metadata()
        data["dataFiles"][0]["value"] = (
            "https://download.iprox.org/IPX0000001000/%2e%2e/results_%FD.plink"
        )
        with TemporaryDirectory() as temporary:
            output = Path(temporary) / "project"
            with self.assertRaises(ValueError):
                write_repository_outputs("iProX", "PXD000001", snapshot(data), [], output)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
