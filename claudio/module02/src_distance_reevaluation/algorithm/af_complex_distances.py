"""Calculate distances for strictly mapped endpoints in one AF complex."""

from math import isfinite

import pandas as pd

from claudio.utils.utils import build_xl_dataset
from claudio.module02.src_distance_reevaluation.algorithm.af_complex_mapping import (
    prepare_complex_distance_input,
)
from claudio.module02.src_distance_reevaluation.algorithm.calc_site_distances import (
    calculate_site_dists,
)


def calculate_complex_site_dists(
    data, temp_dir, xl_residues, plddt_cutoff, topolink_bin, verbose_level
):
    cutoff = float(plddt_cutoff)
    if not isfinite(cutoff) or not 0 <= cutoff <= 100:
        raise ValueError("pLDDT cutoff must be between 0 and 100")

    rules = build_xl_dataset(xl_residues)
    result = prepare_complex_distance_input(data, rules)
    for prefix in ("eucl", "topo"):
        result[f"{prefix}_dist_tplk"] = float("nan")
        result[f"{prefix}_dist_lower_bound"] = float("nan")
        result[f"{prefix}_dist_raw"] = pd.Series(
            None, index=result.index, dtype="object",
        )
    result["distance_status"] = "not_evaluable"
    result["distance_reason"] = pd.Series(
        None, index=result.index, dtype="object"
    )

    eligible = []
    for position, (_, row) in enumerate(result.iterrows()):
        if not row["res_criteria_fulfilled"]:
            reason = (
                row["mapping_reason_a"]
                or row["mapping_reason_b"]
                or "unmapped_endpoint"
            )
        elif any(
            not isfinite(value)
            for site in ("a", "b")
            for value in row[f"mapping_coord_{site}"]
        ):
            reason = "invalid_coordinates"
        elif any(
            not isfinite(row[f"pLDDT_{site}"])
            or not 0 <= row[f"pLDDT_{site}"] <= 100
            for site in ("a", "b")
        ):
            reason = "invalid_plddt"
        elif any(row[f"pLDDT_{site}"] < cutoff for site in ("a", "b")):
            reason = "low_plddt"
        elif (
            row["chain_a"] == row["chain_b"]
            and row["mapped_residue_id_a"] == row["mapped_residue_id_b"]
            and row["mapping_atom_a"] == row["mapping_atom_b"]
        ):
            reason = "identical_atom_endpoints"
        else:
            eligible.append(position)
            continue
        result.iat[position, result.columns.get_loc("distance_reason")] = reason

    if not eligible:
        return result

    selected = result.iloc[eligible].copy().reset_index(drop=True)
    # Keep distant pairs within the selected complex so distance
    # violations are not removed by the legacy proximity prefilter.
    selected["is_interfaced"] = True
    computed = calculate_site_dists(
        selected, temp_dir, rules, cutoff, topolink_bin, verbose_level
    )

    distance_columns = [
        f"{prefix}_dist{suffix}"
        for prefix in ("eucl", "topo")
        for suffix in ("_tplk", "_lower_bound", "_raw")
    ]
    for local_position, original_position in enumerate(eligible):
        computed_row = computed.loc[local_position]
        for column in distance_columns:
            result.iat[
                original_position, result.columns.get_loc(column)
            ] = computed_row.get(column, float("nan"))

        available = all(
            pd.notna(computed_row.get(column))
            and isfinite(computed_row[column])
            for column in ("eucl_dist_tplk", "topo_dist_tplk")
        )
        has_bound = any(
            pd.notna(computed_row.get(column))
            and isfinite(computed_row[column])
            for column in ("eucl_dist_lower_bound", "topo_dist_lower_bound")
        )
        if has_bound:
            status, reason = "lower_bound", "topolink_lower_bound"
        elif available:
            status, reason = "computed", None
        else:
            status, reason = "not_evaluable", "distance_unavailable"
        result.iat[
            original_position, result.columns.get_loc("distance_status")
        ] = status
        result.iat[
            original_position, result.columns.get_loc("distance_reason")
        ] = reason

    return result