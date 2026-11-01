import json
import csv
import hashlib
import unittest
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
from urllib.parse import parse_qs, urlsplit

import requests

from claudio.dataset_catalog.repositories import RepositoryResponse, write_repository_outputs


def massive_metadata():
    return {
        "title": "Cross-linking dataset",
        "accession": [
            {"accession": "MS:1001919", "value": "PXD000001"},
            {"accession": "MS:1002487", "value": "MSV000012345"},
        ],
    }


def massive_file(path="raw/sample.raw", **fields):
    return {
        "file_descriptor": "f.MSV000012345/" + path,
        "name": path.rsplit("/", 1)[-1],
        "collection": path.split("/", 1)[0],
        "size": 123,
        "type": "File",
        **fields,
    }


def http_response(data=None, files=None, total=None):
    response = requests.Response()
    response.status_code = 200
    if files is None:
        response.url = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets/PXD000001"
        response.headers["Content-Type"] = "application/json"
        response._content = (json.dumps(data) + "\n").encode("utf-8")
    else:
        response.url = "https://massive.ucsd.edu/ProteoSAFe/dataset_files.jsp?accession=MSV000012345"
        response.headers["Content-Type"] = "text/html; charset=UTF-8"
        payload = {
            "source": "ProteoSAFe",
            "total_rows": len(files) if total is None else total,
            "row_data": files,
        }
        response._content = (
            "<html><script>\nvar dataset_files = "
            + json.dumps(payload) + ";\n</script></html>\n"
        ).encode("utf-8")
    return response


def snapshot(response, data):
    return RepositoryResponse(
        response.url, "2026-10-31T15:00:00+00:00", response.content,
        response.headers.copy(), data,
    )


class TestMassiveInventory(unittest.TestCase):
    def candidate(self):
        return {
            "accession": "PXD000001", "repository": "MassIVE",
            "title": "Cross-linking dataset", "search_terms": ["XL-MS"],
            "review_status": "pending",
        }

    def inputs(self, records=None):
        metadata = massive_metadata()
        records = [massive_file()] if records is None else records
        project = snapshot(http_response(metadata), metadata)
        data = {"source": "ProteoSAFe", "total_rows": len(records), "row_data": records}
        page = snapshot(http_response(files=records), data)
        return project, [page]

    def test_writes_common_outputs_with_native_identity_and_original_html(self):
        records = [massive_file("raw/a/sample.raw"), massive_file("result/b/sample & #1.tsv")]
        project, pages = self.inputs(records)
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("MassIVE", "PXD000001", project, pages, output_dir)

            manifest = json.loads((output_dir / "manifest.json").read_bytes())
            files = json.loads((output_dir / "files.json").read_bytes())
            self.assertEqual(manifest["source"], "MassIVE")
            self.assertEqual(manifest["repository_accession"], "MSV000012345")
            self.assertEqual(manifest["project_accession"], "PXD000001")
            self.assertEqual(manifest["screening_status"], "pending")
            self.assertEqual(manifest["file_count"], 2)
            self.assertEqual([f["accession"] for f in files], [r["file_descriptor"] for r in records])
            self.assertEqual(files[1]["projectAccessions"], ["PXD000001"])
            self.assertEqual(files[1]["relativePath"], "result/b/sample & #1.tsv")
            location = files[1]["publicFileLocations"][0]["value"]
            self.assertEqual(parse_qs(urlsplit(location).query)["file"], [records[1]["file_descriptor"]])
            saved = manifest["file_responses"][0]
            self.assertTrue(saved["path"].endswith(".html"))
            self.assertEqual((output_dir / saved["path"]).read_bytes(), pages[0].raw_body)
            self.assertEqual(saved["sha256"], hashlib.sha256(pages[0].raw_body).hexdigest())
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[1]["file_name"], "sample & #1.tsv")

    def test_writes_a_header_for_a_confirmed_empty_inventory(self):
        project, pages = self.inputs([])
        with TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / "project"
            write_repository_outputs("MassIVE", "PXD000001", project, pages, output_dir)
            self.assertEqual(json.loads((output_dir / "files.json").read_bytes()), [])
            with (output_dir / "files.csv").open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                self.assertIn("file_accession", reader.fieldnames)
                self.assertEqual(list(reader), [])

    def test_rejects_bad_sources_before_creating_outputs(self):
        for problem in ("project", "missing_files", "raw_html", "different_msv"):
            with self.subTest(problem=problem), TemporaryDirectory() as temp_dir:
                project, pages = self.inputs()
                if problem == "project":
                    project.data["accession"][0]["value"] = "PXD000002"
                elif problem == "missing_files":
                    pages = []
                elif problem == "raw_html":
                    pages[0].data["row_data"][0]["size"] = 999
                else:
                    metadata = deepcopy(project.data)
                    metadata["accession"][1]["value"] = "MSV000054321"
                    project = snapshot(http_response(metadata), metadata)
                output_dir = Path(temp_dir) / "project"
                with self.assertRaises(ValueError):
                    write_repository_outputs("MassIVE", "PXD000001", project, pages, output_dir)
                self.assertFalse(output_dir.exists())


if __name__ == "__main__":
    unittest.main()
