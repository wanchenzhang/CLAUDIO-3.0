"""Assess one distance method against an explicitly configured interval."""
from math import isfinite
from numbers import Real

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