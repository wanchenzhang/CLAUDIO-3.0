"""Classify relationships within resolved model candidates, not experimental truth."""
from math import isfinite
from numbers import Real

import numpy as np
import pandas as pd

from claudio.module04.src.algorithm.structure_assessment import (
    annotate_structure_relationships,
)


def _missing(value):
    return pd.isna(value) or (
        isinstance(value, str) and value.strip() in ("", "-")
    )


def _protein_relation(row):
    first, second = row.get("unip_id_a"), row.get("unip_id_b")
    if _missing(first) or _missing(second):
        return "unknown"
    return (
        "same_accession"
        if str(first).strip() == str(second).strip()
        else "different_accessions"
    )


def _classify(row):
    relation = row["model_chain_relation"]
    if relation == "unknown":
        return "undetermined", row["model_chain_relation_reason"], None
    if row["protein_relationship"] == "unknown":
        return "undetermined", "missing_protein_identity", None

    for site in ("a", "b"):
        status = row.get(f"position_status_{site}")
        if _missing(status):
            return "undetermined", "missing_position_status", site
        if status in ("invalid", "ambiguous"):
            return "undetermined", "unresolved_position", site
        if status not in ("valid", "corrected"):
            return "undetermined", "unrecognized_position_status", site

        mapping = row.get(f"mapping_status_{site}")
        if _missing(mapping):
            if row["structure_source"] == "af_complex":
                return "undetermined", "missing_mapping_status", site
        elif mapping != "mapped":
            return "undetermined", "unmapped_endpoint", site

        confirmed = row.get(f"res_crit_{site}")
        if not isinstance(confirmed, (bool, np.bool_)) or not confirmed:
            return "undetermined", "unmapped_endpoint", site

        number = row.get(f"pdb_pos_{site}")
        if _missing(number):
            return "undetermined", "missing_residue_identity", site
        if (
            not isinstance(number, Real) or isinstance(number, bool)
            or not isfinite(number) or int(number) != number
        ):
            return "undetermined", "invalid_residue_identity", site

    reason = row.get("distance_reason")
    if isinstance(reason, str) and reason == "identical_atom_endpoints":
        return "undetermined", "identical_atom_endpoints", None

    if relation == "same_chain":
        if row["protein_relationship"] == "different_accessions":
            return "undetermined", "protein_chain_conflict", None

        if row["pdb_pos_a"] == row["pdb_pos_b"]:
            fields = ("pdb_hetflag", "pdb_icode")
            identity = [
                tuple(row.get(f"{field}_{site}") for field in fields)
                for site in ("a", "b")
            ]
            if any(
                pd.isna(value) or value == "-"
                for key in identity for value in key
            ):
                return "undetermined", "incomplete_residue_identity", None

            identity = [
                tuple(str(value).strip() for value in key)
                for key in identity
            ]
            if identity[0] == identity[1]:
                atom_a = row.get("mapping_atom_a")
                atom_b = row.get("mapping_atom_b")
                if _missing(atom_a) or _missing(atom_b) or atom_a == atom_b:
                    return "undetermined", "coincident_residue_endpoints", None

        return "intra", "mapped_same_chain", None

    return "inter", "mapped_cross_chain", None


def classify_model_candidates(data):
    if not data.columns.is_unique:
        raise ValueError("Candidate input column names must be unique")

    result = annotate_structure_relationships(data)
    result["protein_relationship"] = pd.Series(
        [_protein_relation(row) for _, row in result.iterrows()],
        index=result.index,
        dtype="object",
    )
    decisions = [_classify(row) for _, row in result.iterrows()]

    for position, column in enumerate(
        ("XL_type", "XL_type_reason", "XL_type_site")
    ):
        result[column] = pd.Series(
            [decision[position] for decision in decisions],
            index=result.index,
            dtype="object",
        )

    result["XL_type_scope"] = "model_candidate"
    result["XL_type_rule_version"] = "model-candidate-type-v1"
    return result