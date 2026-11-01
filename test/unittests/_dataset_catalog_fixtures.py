"""Small fake pridepy backend for the existing catalog workflow fixtures."""
from unittest.mock import Mock, patch

from claudio.dataset_catalog import repositories


def install_pridepy_fixture(case):
    def request(url, **kwargs):
        if not isinstance(repositories.requests.get, Mock):
            raise AssertionError("This fixture must never make a real HTTP request")
        response = repositories.requests.get(url, **kwargs)
        response.raise_for_status()
        return response.json()

    class Provider:
        def list_files(self, accession):
            return request(f"{repositories.PRIDE_API_URL}/projects/{accession}/files",
                           params={"pageSize": 100, "page": 0})

    def provider(repository):
        if repository != "PRIDE":
            raise AssertionError("Use explicit native-provider fixtures in repository tests")
        return Provider()

    for target, replacement in [
        ("_pridepy_version", lambda: "0.0.16"),
        ("_pridepy_project", lambda accession: request(f"{repositories.PRIDE_API_URL}/projects/{accession}")),
        ("_pridepy_provider", provider),
    ]:
        patcher = patch.object(repositories, target, replacement)
        patcher.start()
        case.addCleanup(patcher.stop)
