import unittest
from unittest.mock import Mock, patch

from claudio.dataset_catalog import repositories


class TestPridepyDependencies(unittest.TestCase):
    @patch("importlib.metadata.version", return_value="0.0.16")
    def test_records_the_pinned_library_version(self, version):
        self.assertEqual(repositories._pridepy_version(), "0.0.16")

    @patch("importlib.metadata.version", return_value="0.0.15")
    def test_different_versions_need_an_explicit_adapter_update(self, version):
        with self.assertRaisesRegex(RuntimeError, "0.0.16"):
            repositories._pridepy_version()

    @patch.object(repositories, "_pridepy_provider")
    def test_invalid_page_size_is_rejected_before_listing(self, provider):
        for value in [0, -1, True, 1.5]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                repositories.fetch_pride_file_pages("PXD000001", page_size=value)
        provider.assert_not_called()

    @patch.object(repositories, "_pridepy_version", return_value="0.0.16")
    @patch.object(repositories, "_pridepy_provider")
    def test_listing_is_delegated_without_reimplementing_pagination(self, provider, version):
        provider.return_value.list_files.return_value = []
        pages = repositories.fetch_pride_file_pages("PXD000001")
        self.assertEqual(len(pages), 1)
        self.assertEqual(pages[0].data["records"], [])
        provider.return_value.list_files.assert_called_once_with("PXD000001")
