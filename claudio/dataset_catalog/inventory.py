"""Collect project inventories while retaining candidate provenance."""

from copy import deepcopy
from datetime import datetime, timezone

import requests

from .pride import fetch_pride_project, fetch_pride_file_pages
from .pride_outputs import write_pride_outputs


def inspect_candidate(candidate, output_dir):
    """Inspect one candidate and return its technical collection status.

    Parameters
    ----------
    candidate : dict
        Candidate record from the discovery catalog.
    output_dir : str or pathlib.Path
        New directory for this project's inventory.

    Returns
    -------
    result : dict
        Candidate metadata with inspection status and error information.
    """
    result = deepcopy(candidate)
    result.update({
        "inspection_status": "adapter_pending",
        "file_count": None,
        "inventory_path": "",
        "error_stage": "",
        "error_type": "",
        "error_message": "",
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "finished_at_utc": "",
    })

    if candidate.get("repository") != "PRIDE":
        result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
        return result

    stage = "candidate"
    try:
        accession = candidate.get("accession")
        if not isinstance(accession, str) or not accession.strip():
            raise ValueError("Candidate has no valid project accession")

        stage = "project_metadata"
        project = fetch_pride_project(accession)

        stage = "file_list"
        pages = fetch_pride_file_pages(accession)

        stage = "output"
        destination = write_pride_outputs(
            accession, project, pages, output_dir
        )
        result.update({
            "inspection_status": "collected",
            "file_count": sum(len(page.data) for page in pages),
            "inventory_path": str(destination),
        })
    except (requests.exceptions.RequestException, ValueError, OSError) as error:
        result.update({
            "inspection_status": "failed",
            "error_stage": stage,
            "error_type": type(error).__name__,
            "error_message": str(error),
        })

    result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    return result