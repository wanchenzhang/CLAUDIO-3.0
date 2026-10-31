"""Combine dataset candidates found by different searches."""


def merge_search_results(searches):
    """Merge candidates and retain metadata variants with their search sources."""
    candidates = {}

    for search_term, pages in searches.items():
        for page in pages:
            columns = page.data["result_set"]["datasets_title_list"]

            for values in page.data["datasets"]:
                if len(values) != len(columns):
                    raise ValueError(
                        "Dataset row does not match its column headers"
                    )

                metadata = dict(zip(columns, values))
                accession = metadata.pop("dataset identifier")

                if accession not in candidates:
                    candidates[accession] = {
                        **metadata,
                        "accession": accession,
                        "search_terms": [],
                        "review_status": "pending",
                        "metadata_conflict": False,
                        "metadata_variants": [],
                    }

                candidate = candidates[accession]
                if search_term not in candidate["search_terms"]:
                    candidate["search_terms"].append(search_term)

                variants = candidate["metadata_variants"]
                for variant in variants:
                    if variant["metadata"] == metadata:
                        if search_term not in variant["search_terms"]:
                            variant["search_terms"].append(search_term)
                        break
                else:
                    variants.append({
                        "metadata": metadata,
                        "search_terms": [search_term],
                    })

                candidate["metadata_conflict"] = len(variants) > 1

    return [candidates[accession] for accession in sorted(candidates)]