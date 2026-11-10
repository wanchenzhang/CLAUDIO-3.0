
"""Publication columns retained in per-interaction assessment outputs."""

PUBLICATION_COLUMNS = (
    "xlms_dataset_accession",
    "xlms_publication_status",
    "xlms_paper_dois_json",
    "xlms_paper_pmids_json",
    "xlms_papers_json",
)

PUBLICATION_SOURCE_COLUMNS = (
    "dataset_accession",
    "xlms_source_record_id",
    "source_record_id",
)


def publication_output_columns(data, base_columns):
    """Keep available paper fields and their row identifiers."""
    columns = list(base_columns)
    if not any(
        name in data.columns for name in PUBLICATION_COLUMNS
    ):
        return columns
    for name in (
        PUBLICATION_COLUMNS + PUBLICATION_SOURCE_COLUMNS
    ):
        if name in data.columns and name not in columns:
            columns.append(name)
    return columns
