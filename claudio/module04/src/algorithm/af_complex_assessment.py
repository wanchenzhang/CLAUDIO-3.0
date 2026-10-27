"""Describe AF model chain assignments without reclassifying source records."""
import pandas as pd


def _label(value):
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text if text not in ("", "-") else None


def _chain_relationship(row):
    status = _label(row.get("af_match_status"))
    if status is None:
        return "unknown", "missing_af_match_status"
    if status != "matched":
        reason = (
            status if status in ("no_candidate_model", "no_matching_chain_pair")
            else "unusable_af_match_status"
        )
        return "unknown", reason

    if (
        _label(row.get("af_model_id")) is None
        or _label(row.get("af_model_version")) is None
    ):
        return "unknown", "missing_model_identity"
    if _label(row.get("path")) is None:
        return "unknown", "no_usable_structure"

    chain_a = _label(row.get("chain_a"))
    chain_b = _label(row.get("chain_b"))
    if chain_a is None or chain_b is None:
        return "unknown", "missing_chain_assignment"
    if chain_a == chain_b:
        return "same_chain", "same_model_chain"
    return "cross_chain", "different_model_chains"


def annotate_complex_relationships(data):
    """Describe existing AF assignments; retain sources and distance results."""
    result = data.copy()
    relationships = [_chain_relationship(row) for _, row in result.iterrows()]
    result["model_chain_relation"] = pd.Series(
        [relation for relation, _ in relationships],
        index=result.index, dtype="object",
    )
    result["model_chain_relation_reason"] = pd.Series(
        [reason for _, reason in relationships],
        index=result.index, dtype="object",
    )
    result["chain_relation_rule_version"] = "af-chain-v1"
    return result