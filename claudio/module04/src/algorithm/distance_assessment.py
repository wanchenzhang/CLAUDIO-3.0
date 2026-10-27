"""Assess one distance method against an explicitly configured interval."""
from math import isfinite
from numbers import Real
import numpy as np
import pandas as pd


def _finite_nonnegative(value):
    return (
        isinstance(value, Real)
        and not isinstance(value, bool)
        and isfinite(value)
        and value >= 0
    )


def assess_distance_interval(
    exact_distance, lower_bound, *, minimum, maximum, boundary,
):
    """Treat lower_bound as a strict > value, never as an exact distance."""
    if (
        not _finite_nonnegative(minimum)
        or not _finite_nonnegative(maximum)
        or minimum >= maximum
    ):
        raise ValueError("Invalid distance interval")
    if boundary not in ("inclusive", "exclusive"):
        raise ValueError("Unsupported interval boundary policy")

    has_exact = not pd.isna(exact_distance)
    has_bound = not pd.isna(lower_bound)
    values = [
        value
        for value, present in (
            (exact_distance, has_exact), (lower_bound, has_bound),
        )
        if present
    ]

    if any(not _finite_nonnegative(value) for value in values):
        status, reason = "not_evaluable", "invalid_distance_value"
    elif not values:
        status, reason = "not_evaluable", "distance_unavailable"
    elif has_exact and has_bound:
        status, reason = "not_evaluable", "conflicting_distance_values"
    elif has_exact:
        below = exact_distance < minimum or (
            boundary == "exclusive" and exact_distance == minimum
        )
        above = exact_distance > maximum or (
            boundary == "exclusive" and exact_distance == maximum
        )
        if below:
            status, reason = "violated", "exact_below_range"
        elif above:
            status, reason = "violated", "exact_above_range"
        else:
            status, reason = "satisfied", "exact_within_range"
    elif lower_bound >= maximum:
        status, reason = "violated", "lower_bound_above_range"
    else:
        status, reason = "indeterminate", "lower_bound_inconclusive"

    return {
        "status": status,
        "reason": reason,
        "minimum": float(minimum),
        "maximum": float(maximum),
        "boundary": boundary,
        "rule_version": "distance-interval-v1",
    }
def _missing_endpoint_value(value):
    return pd.isna(value) or (
        isinstance(value, str) and value.strip() in ("", "-")
    )


def check_endpoint_eligibility(row, *, structure_source, plddt_cutoff):
    """Check M1 status, M2 endpoint mapping and AF confidence only."""
    if structure_source not in ("pdb", "af_monomer", "af_complex"):
        raise ValueError("Unsupported structure_source")
    if not _finite_nonnegative(plddt_cutoff) or plddt_cutoff > 100:
        raise ValueError("pLDDT cutoff must be between 0 and 100")

    def outcome(reason, site=None):
        upstream = (
            row.get(f"mapping_reason_{site}")
            if reason == "unmapped_endpoint" else None
        )
        return {
            "eligible": reason == "endpoints_evaluable",
            "reason": reason,
            "site": site,
            "upstream_reason": (
                None if _missing_endpoint_value(upstream) else upstream
            ),
            "structure_source": structure_source,
            "plddt_cutoff": float(plddt_cutoff),
            "rule_version": "endpoint-eligibility-v1",
        }

    for site in ("a", "b"):
        status = row.get(f"position_status_{site}")
        if _missing_endpoint_value(status):
            return outcome("missing_position_status", site)
        if status in ("invalid", "ambiguous"):
            return outcome("unresolved_position", site)
        if status not in ("valid", "corrected"):
            return outcome("unrecognized_position_status", site)

    for site in ("a", "b"):
        mapping = row.get(f"mapping_status_{site}")
        if _missing_endpoint_value(mapping):
            if structure_source == "af_complex":
                return outcome("missing_mapping_status", site)
        elif mapping != "mapped":
            return outcome("unmapped_endpoint", site)

        confirmed = row.get(f"res_crit_{site}")
        if not isinstance(confirmed, (bool, np.bool_)) or not confirmed:
            return outcome("unmapped_endpoint", site)
        if _missing_endpoint_value(row.get(f"chain_{site}")):
            return outcome("missing_chain_assignment", site)

        number = row.get(f"pdb_pos_{site}")
        if _missing_endpoint_value(number):
            return outcome("missing_residue_identity", site)
        if (
            not isinstance(number, Real)
            or isinstance(number, bool)
            or not isfinite(number)
            or int(number) != number
        ):
            return outcome("invalid_residue_identity", site)

    if structure_source != "pdb":
        for site in ("a", "b"):
            confidence = row.get(f"pLDDT_{site}")
            if _missing_endpoint_value(confidence):
                return outcome("missing_plddt", site)
            if isinstance(confidence, (bool, np.bool_)):
                return outcome("invalid_plddt", site)
            try:
                confidence = float(confidence)
            except (TypeError, ValueError, OverflowError):
                return outcome("invalid_plddt", site)
            if not isfinite(confidence) or not 0 <= confidence <= 100:
                return outcome("invalid_plddt", site)
            if confidence < plddt_cutoff:
                return outcome("low_plddt", site)

    return outcome("endpoints_evaluable")