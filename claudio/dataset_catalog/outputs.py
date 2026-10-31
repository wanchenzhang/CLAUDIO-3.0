"""Write dataset discovery results and their provenance files."""

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from .discovery import merge_search_results


CSV_COLUMNS = [
    "accession",
    "title",
    "repository",
    "species",
    "search_terms",
    "review_status",
    "metadata_conflict",
    "metadata_variant_count",
]


def _write_json(path, data):
    path.write_text(
        json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_discovery_outputs(searches, output_dir):
    """Write a discovery snapshot into a new directory.

    Parameters
    ----------
    searches : dict
        Search terms mapped to their retrieved CatalogPage objects.
    output_dir : str or pathlib.Path
        New directory for the snapshot.

    Returns
    -------
    output_dir : pathlib.Path
    """
    if not searches or any(not pages for pages in searches.values()):
        raise ValueError("Each search must include at least one response page")

    candidates = merge_search_results(searches)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "raw").mkdir()

    manifest = {
        "schema_version": 1,
        "source": "ProteomeCentral",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "candidate_count": len(candidates),
        "searches": [],
    }

    for search_number, (search_term, pages) in enumerate(
        searches.items(), start=1
    ):
        search = {
            "search_term": search_term,
            "reported_total": pages[0].data["result_set"]["n_available_rows"],
            "downloaded_records": sum(
                len(page.data["datasets"]) for page in pages
            ),
            "pages": [],
        }

        for page_number, page in enumerate(pages, start=1):
            relative_path = Path("raw") / (
                f"search_{search_number:03d}_page_{page_number:04d}.json"
            )
            (output_dir / relative_path).write_bytes(page.raw_body)
            search["pages"].append({
                "path": relative_path.as_posix(),
                "request_url": page.request_url,
                "retrieved_at_utc": page.retrieved_at_utc,
                "query": page.data["query"],
                "sha256": hashlib.sha256(page.raw_body).hexdigest(),
            })

        manifest["searches"].append(search)

    _write_json(output_dir / "candidates.json", candidates)

    with (output_dir / "candidates.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for candidate in candidates:
            row = {
                column: candidate.get(column, "")
                for column in CSV_COLUMNS
            }
            row["search_terms"] = json.dumps(
                candidate["search_terms"], ensure_ascii=False
            )
            row["metadata_variant_count"] = len(
                candidate["metadata_variants"]
            )
            writer.writerow(row)

    _write_json(output_dir / "manifest.json", manifest)
    return output_dir