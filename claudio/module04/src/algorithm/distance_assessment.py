"""Assess one distance method against an explicitly configured interval."""
from math import isfinite
from numbers import Real
import numpy as np
import pandas as pd
from collections.abc import Mapping
from claudio.module04.src.algorithm.structure_assessment import (
    annotate_structure_relationships,
)
from claudio.module04.src.algorithm.relationship_classification import (
    classify_model_candidates,
)
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
def assess_candidate_distance(
    row, *, metric, minimum, maximum, boundary, plddt_cutoff,
):
    """Assess one method using existing M2 assignments and stage results."""
    prefixes = {"euclidean": "eucl", "topological": "topo"}
    if metric not in prefixes:
        raise ValueError("Unsupported distance metric")

    result = assess_distance_interval(
        None, None, minimum=minimum, maximum=maximum, boundary=boundary,
    )
    if not _finite_nonnegative(plddt_cutoff) or plddt_cutoff > 100:
        raise ValueError("pLDDT cutoff must be between 0 and 100")

    context = (
        classify_model_candidates(pd.DataFrame([row]))
        if "XL_type_scope" in row
        else annotate_structure_relationships(pd.DataFrame([row]))
    ).iloc[0]
    source = context["structure_source"]
    stage_status = row.get("distance_status")
    if _missing_endpoint_value(stage_status):
        stage_status = None

    result.update({
        "metric": metric,
        "structure_source": source,
        "model_chain_relation": context["model_chain_relation"],
        "chain_relation_rule_version": context["chain_relation_rule_version"],
        "plddt_cutoff": float(plddt_cutoff),
        "upstream_distance_status": stage_status,
        "site": None,
        "upstream_reason": None,
        "interval_rule_version": result["rule_version"],
        "endpoint_rule_version": None,
        "relationship_rule_version": context.get("XL_type_rule_version"),
        "rule_version": (
            "candidate-distance-v2" if "XL_type_scope" in row
            else "candidate-distance-v1"
        ),
    })

    def stop(reason, stage, *, site=None, upstream=None):
        result.update(
            status="not_evaluable", reason=reason, assessment_stage=stage,
            site=site, upstream_reason=upstream,
        )
        return result

    if context["model_chain_relation"] == "unknown":
        return stop(context["model_chain_relation_reason"], "structure")

    gate = check_endpoint_eligibility(
        row, structure_source=source, plddt_cutoff=plddt_cutoff,
    )
    result["endpoint_rule_version"] = gate["rule_version"]
    if not gate["eligible"]:
        return stop(
            gate["reason"], "endpoints",
            site=gate["site"], upstream=gate["upstream_reason"],
        )

    if stage_status is None:
        if source == "af_complex":
            return stop("missing_distance_status", "distance")
    elif stage_status == "not_evaluable":
        upstream = row.get("distance_reason")
        return stop(
            "upstream_distance_not_evaluable", "distance",
            upstream=None if _missing_endpoint_value(upstream) else upstream,
        )
    elif stage_status not in ("computed", "lower_bound"):
        return stop("unrecognized_distance_status", "distance")
    if context.get("XL_type") == "undetermined" and "XL_type_scope" in row:
        return stop(context["XL_type_reason"], "relationship")
    bounds = [row.get(f"{prefix}_dist_lower_bound") for prefix in ("eucl", "topo")]
    exact = [row.get(f"{prefix}_dist") for prefix in ("eucl", "topo")]
    if stage_status == "computed" and (
        any(pd.notna(value) for value in bounds)
        or not all(_finite_nonnegative(value) for value in exact)
    ):
        return stop("inconsistent_distance_status", "distance")
    if stage_status == "lower_bound" and not any(
        _finite_nonnegative(value) for value in bounds
    ):
        return stop("inconsistent_distance_status", "distance")

    prefix = prefixes[metric]
    interval = assess_distance_interval(
        row.get(f"{prefix}_dist"), row.get(f"{prefix}_dist_lower_bound"),
        minimum=minimum, maximum=maximum, boundary=boundary,
    )
    result.update(
        status=interval["status"], reason=interval["reason"],
        assessment_stage="interval",
    )
    return result
def annotate_distance_assessments(
    data, *, euclidean_policy, topological_policy, plddt_cutoff,
):
    """Append per-method assessments while retaining every candidate row."""
    methods = (
        ("euclidean", "eucl", euclidean_policy),
        ("topological", "topo", topological_policy),
    )
    required = {"minimum", "maximum", "boundary"}
    for metric, _, policy in methods:
        if not isinstance(policy, Mapping) or set(policy) != required:
            raise ValueError(
                f"{metric} policy requires minimum, maximum and boundary"
            )
        assess_distance_interval(None, None, **policy)
    if not _finite_nonnegative(plddt_cutoff) or plddt_cutoff > 100:
        raise ValueError("pLDDT cutoff must be between 0 and 100")
    if not data.columns.is_unique:
        raise ValueError("Candidate input column names must be unique")

    result = (
        classify_model_candidates(data)
        if "XL_type_scope" in data.columns
        else annotate_structure_relationships(data)
    )
    fields = {
        "status": "status",
        "reason": "reason",
        "stage": "assessment_stage",
        "site": "site",
        "upstream_reason": "upstream_reason",
        "minimum": "minimum",
        "maximum": "maximum",
        "boundary": "boundary",
        "plddt_cutoff": "plddt_cutoff",
        "upstream_distance_status": "upstream_distance_status",
        "metric": "metric",
        "rule_version": "rule_version",
        "interval_rule_version": "interval_rule_version",
        "endpoint_rule_version": "endpoint_rule_version",
        "chain_relation_rule_version": "chain_relation_rule_version",
        "relationship_rule_version": "relationship_rule_version",
    }
    for metric, prefix, policy in methods:
        decisions = [
            assess_candidate_distance(
                row, metric=metric, plddt_cutoff=plddt_cutoff, **policy,
            )
            for _, row in result.iterrows()
        ]
        for column, field in fields.items():
            dtype = (
                "float64"
                if column in ("minimum", "maximum", "plddt_cutoff")
                else "object"
            )
            result[f"{prefix}_assessment_{column}"] = pd.Series(
                [decision[field] for decision in decisions],
                index=result.index,
                dtype=dtype,
            )
    return result