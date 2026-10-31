"""Retrieve project metadata and file listings from PRIDE."""

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Mapping, Union

import requests


PRIDE_API_URL = "https://www.ebi.ac.uk/pride/ws/archive/v3"


@dataclass(frozen=True)
class PrideResponse:
    request_url: str
    retrieved_at_utc: str
    raw_body: bytes
    headers: Mapping[str, str]
    data: Union[dict, list]


def fetch_pride_project(accession: str) -> PrideResponse:
    """Fetch a PRIDE project and retain its original response.

    Parameters
    ----------
    accession : str
        ProteomeXchange project accession.

    Returns
    -------
    project : PrideResponse
    """
    response = requests.get(
        f"{PRIDE_API_URL}/projects/{accession}", timeout=(15, 120)
    )
    response.raise_for_status()
    retrieved_at = datetime.now(timezone.utc).isoformat()
    data = json.loads(response.content)

    if not isinstance(data, dict):
        raise ValueError("PRIDE project response must be a JSON object")
    if data.get("accession") != accession:
        raise ValueError(
            f"PRIDE project accession mismatch: requested {accession}, "
            f"received {data.get('accession')}"
        )

    return PrideResponse(
        request_url=response.url,
        retrieved_at_utc=retrieved_at,
        raw_body=response.content,
        headers=response.headers.copy(),
        data=data,
    )