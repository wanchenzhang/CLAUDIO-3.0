"""Describe selected structure chains independently of source XL classifications."""
import pandas as pd


def _label(value):
    if value is None or pd.isna(value):
        return None
    text = str(value).strip()
    return text if text not in ("", "-") else None


def _infer_source(row):
    if any(
        _label(row.get(column)) is not None
        for column in ("af_match_status", "af_model_id", "af_model_version")
    ):
        return "af_complex"

    methods = [
        _label(row.get(column))
        for column in ("pdb_method", "method_a", "method_b")
    ]
    if any(method and method.casefold() == "alphafold" for method in methods):
        return "af_monomer"

    if (
        _label(row.get("pdb_id")) is not None
        or _label(row.get("pdb_method")) is not None
    ):
        return "pdb"
    return "unknown"


def _chain_relationship(row, source):
    if source == "af_complex":
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

    elif source in ("pdb", "af_monomer"):
        if _label(row.get("pdb_id")) is None:
            return "unknown", "missing_model_identity"
    else:
        return "unknown", "unknown_structure_source"

    if _label(row.get("path")) is None:
        return "unknown", "no_usable_structure"

    chain_a = _label(row.get("chain_a"))
    chain_b = _label(row.get("chain_b"))
    if chain_a is None or chain_b is None:
        return "unknown", "missing_chain_assignment"
    if chain_a == chain_b:
        return "same_chain", "same_model_chain"
    return "cross_chain", "different_model_chains"


def annotate_structure_relationships(data, *, structure_source=None):
    """Annotate existing assignments without changing sources or distances."""
    if structure_source is not None and (
        not isinstance(structure_source, str)
        or structure_source not in ("pdb", "af_monomer", "af_complex", "unknown")
    ):
        raise ValueError("Unsupported structure_source")

    result = data.copy()
    sources = []
    relationships = []
    for _, row in result.iterrows():
        source = (
            structure_source
            if structure_source is not None
            else _infer_source(row)
        )
        sources.append(source)
        relationships.append(_chain_relationship(row, source))

    result["structure_source"] = pd.Series(
        sources, index=result.index, dtype="object",
    )
    result["model_chain_relation"] = pd.Series(
        [relation for relation, _ in relationships],
        index=result.index, dtype="object",
    )
    result["model_chain_relation_reason"] = pd.Series(
        [reason for _, reason in relationships],
        index=result.index, dtype="object",
    )
    result["chain_relation_rule_version"] = "structure-chain-v1"
    return result