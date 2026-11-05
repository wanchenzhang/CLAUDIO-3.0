import hashlib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from claudio.dataset_catalog import repositories


class TestPridepyDownloads(unittest.TestCase):
    def item(self, size=8, repository="PRIDE"):
        return {"pilot_file_id": "pilot_001", "project_accession": "PXD000001",
                "repository": repository, "reported_file_name": "sites.csv", "pilot_role": "result_table_candidate",
                "download_urls": ["https://example.org/download?file=sites.csv"], "reported_size_bytes": size}

    def transfer(self, body):
        def save(url, folder):
            (Path(folder) / "download").write_bytes(body)
        return save

    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    @patch.object(repositories, "_pridepy_download")
    def test_download_is_delegated_and_atomically_published(self, download, version):
        download.side_effect = self.transfer(b"a,b\n1,2\n")
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "sites.csv"
            result = repositories.download_pilot_file(self.item(), target, 64, 8)
            self.assertEqual(result["download_status"], "downloaded", result)
            self.assertEqual(target.read_bytes(), b"a,b\n1,2\n")
            self.assertEqual(result["file_sha256"], hashlib.sha256(target.read_bytes()).hexdigest())
            self.assertEqual(result["charged_payload_bytes"], 8)
            self.assertEqual(result["download_backend"], "pridepy")
            self.assertFalse(result["response_metadata_available"])
            self.assertFalse(target.with_name("sites.csv.part").exists())
        download.assert_called_once()

    @patch.object(repositories, "_pridepy_download")
    def test_unknown_size_does_not_start_transfer(self, download):
        with TemporaryDirectory() as tmp:
            result = repositories.download_pilot_file(self.item(None), Path(tmp) / "sites.csv", 64, 64)
            self.assertEqual(result["download_status"], "skipped_unknown_size")
            self.assertEqual(result["charged_payload_bytes"], 0)
        download.assert_not_called()

    @patch.object(repositories, "_pridepy_download")
    def test_known_size_must_fit_plan_limits(self, download):
        for limit, reservation, status in [(7, 8, "skipped_file_limit"), (64, 7, "budget_exhausted")]:
            with self.subTest(status=status), TemporaryDirectory() as tmp:
                result = repositories.download_pilot_file(self.item(), Path(tmp) / "sites.csv", limit, reservation)
                self.assertEqual(result["download_status"], status)
                self.assertEqual(result["charged_payload_bytes"], 0)
        download.assert_not_called()

    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    @patch.object(repositories, "_pridepy_download")
    def test_partial_failure_never_leaves_a_final_cache_file(self, download, version):
        def fail(url, folder):
            (Path(folder) / "download").write_bytes(b"partial")
            raise RuntimeError("server unavailable")
        download.side_effect = fail
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "sites.csv"
            result = repositories.download_pilot_file(self.item(), target, 64, 8)
            self.assertEqual(result["download_status"], "failed")
            self.assertFalse(target.exists())
            self.assertFalse(target.with_name("sites.csv.part").exists())
            self.assertEqual(result["charged_payload_bytes"], 8)
            self.assertEqual(result["confirmed_payload_bytes"], 0)

    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    @patch.object(repositories, "_pridepy_download")
    def test_unexpected_size_or_checksum_is_not_published(self, download, version):
        for body, checksum in [(b"short", None), (b"a,b\n1,2\n", "0" * 64)]:
            with self.subTest(body=body), TemporaryDirectory() as tmp:
                download.side_effect = self.transfer(body)
                item = self.item()
                item["expected_download_file_sha256"] = checksum
                target = Path(tmp) / "sites.csv"
                result = repositories.download_pilot_file(item, target, 64, 8)
                self.assertEqual(result["download_status"], "failed")
                self.assertFalse(target.exists())

    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    @patch.object(repositories, "_pridepy_download")
    def test_already_present_payload_is_not_overwritten(self, download, version):
        with TemporaryDirectory() as tmp:
            target = Path(tmp) / "sites.csv"
            target.write_bytes(b"original")
            result = repositories.download_pilot_file(self.item(), target, 64, 8)
            self.assertEqual(result["download_status"], "failed")
            self.assertEqual(target.read_bytes(), b"original")
        download.assert_not_called()

    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    @patch.object(repositories, "_pridepy_download")
    def test_iprox_reported_url_is_kept_when_verified_host_is_used(self, download, version):
        download.side_effect = self.transfer(b"a,b\n1,2\n")
        item = self.item(repository="iProX")
        item["download_urls"] = ["https://download.iprox.org/IPX0005101000/sites.csv"]
        with TemporaryDirectory() as tmp:
            result = repositories.download_pilot_file(item, Path(tmp) / "sites.csv", 64, 8)
            self.assertEqual(result["download_status"], "downloaded", result)
            self.assertEqual(result["requested_url"], item["download_urls"][0])
            self.assertEqual(download.call_args.args[0], "https://download.iprox.cn/IPX0005101000/sites.csv")
    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    def test_massive_ftp_download_uses_tls(self, version):
        url = (
            "ftp://massive-ftp.ucsd.edu/"
            "v11/MSV000100391/search/sites.csv"
        )
        item = self.item(repository="MassIVE")
        item["download_urls"] = [url]

        def save(ftp_urls, output_folder, use_tls=False, **kwargs):
            if ftp_urls != [url] or not use_tls:
                raise RuntimeError("421 TLS is required")
            (Path(output_folder) / "sites.csv").write_bytes(b"a,b\n1,2\n")

        with (
            patch(
                "pridepy.download.client.Client.download_files_by_url",
                side_effect=RuntimeError("421 TLS is required"),
            ),
            patch(
                "pridepy.download.client.Client.download_ftp_urls",
                side_effect=save,
            ),
            TemporaryDirectory() as directory,
        ):
            target = Path(directory) / "sites.csv"
            result = repositories.download_pilot_file(
                item, target, 64, 8
            )

            self.assertEqual(result["download_status"], "downloaded", result)
            self.assertEqual(target.read_bytes(), b"a,b\n1,2\n")
            self.assertEqual(result["requested_url"], url)
            self.assertEqual(result["transfer_url"], url)
            self.assertEqual(result["confirmed_payload_bytes"], 8)
            self.assertFalse(target.with_name("sites.csv.part").exists())

    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    def test_pride_ftp_keeps_the_existing_download_entrypoint(self, version):
        url = "ftp://ftp.pride.ebi.ac.uk/pride/data/sites.csv"
        item = self.item()
        item["download_urls"] = [url]

        def save(urls, output_folder, **kwargs):
            if urls != [url]:
                raise RuntimeError("Unexpected download URL")
            (Path(output_folder) / "sites.csv").write_bytes(b"a,b\n1,2\n")

        with (
            patch(
                "pridepy.download.client.Client.download_files_by_url",
                side_effect=save,
            ),
            patch(
                "pridepy.download.client.Client.download_ftp_urls",
                side_effect=AssertionError("Unexpected transport switch"),
            ),
            TemporaryDirectory() as directory,
        ):
            target = Path(directory) / "sites.csv"
            result = repositories.download_pilot_file(
                item, target, 64, 8
            )

            self.assertEqual(result["download_status"], "downloaded", result)
            self.assertEqual(target.read_bytes(), b"a,b\n1,2\n")