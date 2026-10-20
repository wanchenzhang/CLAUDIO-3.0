"""Map protein positions to residues in existing AF complex models."""

from Bio.PDB import MMCIFParser

from claudio.utils.af_complex import (
    ComplexMetadataMismatchError,
    read_complex_chains,
)
from math import isfinite
from numbers import Integral, Real
import json
from pathlib import Path

import pandas as pd
def map_complex_residue(
    structure_path,
    records,
    protein_id,
    protein_sequence,
    chain_id,
    protein_position,
    atom_name,
    position_status,
):
    """Map a position using an exact full-chain sequence match."""

    def not_evaluable(reason):
        return {
            "mapping_status": "not_evaluable",
            "mapping_reason": reason,
            "protein_id": protein_id,
            "protein_position": protein_position,
            "chain_id": chain_id,
            "residue_id": None,
            "residue_name": None,
            "atom_name": atom_name,
            "coord": None,
        }

    if (
        not isinstance(position_status, str)
        or position_status not in ("valid", "corrected")
    ):
        return not_evaluable("unresolved_position")
    if (
        isinstance(protein_position, bool)
        or not isinstance(protein_position, Real)
    ):
        return not_evaluable("invalid_protein_position")

    if not isinstance(protein_position, Integral):
        if (
            not isfinite(protein_position)
            or protein_position != int(protein_position)
        ):
            return not_evaluable("invalid_protein_position")

    protein_position = int(protein_position)
    try:
        chains = read_complex_chains(structure_path, records)
    except ComplexMetadataMismatchError as error:
        result = not_evaluable("structure_metadata_mismatch")
        result["mapping_detail"] = str(error)
        return result

    chain_info = next(
        (chain for chain in chains if chain["chain_id"] == chain_id),
        None,
    )
    if chain_info is None:
        return not_evaluable("chain_not_found")

    if chain_info["protein_id"] != protein_id:
        return not_evaluable("protein_identity_mismatch")

    if chain_info["sequence"] != protein_sequence:
        return not_evaluable("sequence_mismatch")

    structure = MMCIFParser(QUIET=True).get_structure(
        "af_complex", str(structure_path)
    )
    model = next(structure.get_models())
    chain = model[chain_id]

    residues = [
        residue for residue in chain
        if residue.id[0] == " "
    ]

    if not 1 <= protein_position <= len(residues):
        return not_evaluable("position_out_of_range")

    residue = residues[protein_position - 1]

    if atom_name not in residue:
        result = not_evaluable("missing_atom")
        result["residue_id"] = residue.id
        result["residue_name"] = residue.get_resname()
        return result

    atom = residue[atom_name]

    return {
        "mapping_status": "mapped",
        "mapping_reason": None,
        "protein_id": protein_id,
        "protein_position": protein_position,
        "chain_id": chain.id,
        "residue_id": residue.id,
        "residue_name": residue.get_resname(),
        "atom_name": atom.get_name(),
        "coord": tuple(float(value) for value in atom.get_coord()),
    }
def map_complex_record(row, records, atom_a, atom_b):
    """Map both endpoints while preserving the original record."""
    result = row.copy(deep=True)
    af_status = row.get("af_match_status")

    if isinstance(af_status, str) and af_status in (
        "no_candidate_model",
        "no_matching_chain_pair",
    ):
        for site, atom_name in (("a", atom_a), ("b", atom_b)):
            result[f"mapping_status_{site}"] = "not_evaluable"
            result[f"mapping_reason_{site}"] = af_status
            result[f"mapping_detail_{site}"] = None
            result[f"mapped_residue_id_{site}"] = None
            result[f"mapped_residue_name_{site}"] = None
            result[f"mapping_atom_{site}"] = atom_name
            result[f"mapping_coord_{site}"] = None

        return result
    for site, atom_name in (("a", atom_a), ("b", atom_b)):
        mapped = map_complex_residue(
            structure_path=row["path"],
            records=records,
            protein_id=row[f"unip_id_{site}"],
            protein_sequence=row[f"seq_{site}"],
            chain_id=row[f"chain_{site}"],
            protein_position=row[f"pos_{site}"],
            atom_name=atom_name,
            position_status=row.get(f"position_status_{site}"),
        )

        result[f"mapping_status_{site}"] = mapped["mapping_status"]
        result[f"mapping_reason_{site}"] = mapped["mapping_reason"]
        result[f"mapping_detail_{site}"] = mapped.get("mapping_detail")
        result[f"mapped_residue_id_{site}"] = mapped["residue_id"]
        result[f"mapped_residue_name_{site}"] = mapped["residue_name"]
        result[f"mapping_atom_{site}"] = mapped["atom_name"]
        result[f"mapping_coord_{site}"] = mapped["coord"]

    return result
def map_complex_dataset(data, atom_a, atom_b):
    """Map all candidate rows while preserving input fields and indices."""
    result = data.copy(deep=True)

    mapping_fields = (
        "mapping_status",
        "mapping_reason",
        "mapping_detail",
        "mapped_residue_id",
        "mapped_residue_name",
        "mapping_atom",
        "mapping_coord",
    )
    mapping_columns = [
        f"{field}_{site}"
        for site in ("a", "b")
        for field in mapping_fields
    ]

    for column in mapping_columns:
        result[column] = pd.Series(
            None, index=result.index, dtype="object"
        )

    metadata_cache = {}

    for row_number, (_, row) in enumerate(data.iterrows()):
        af_status = row.get("af_match_status")

        if isinstance(af_status, str) and af_status in (
            "no_candidate_model",
            "no_matching_chain_pair",
        ):
            records = []
        else:
            metadata_path = Path(row["af_metadata_path"])

            if metadata_path not in metadata_cache:
                metadata_cache[metadata_path] = json.loads(
                    metadata_path.read_text(encoding="utf-8")
                )

            records = metadata_cache[metadata_path]

        mapped = map_complex_record(
            row=row,
            records=records,
            atom_a=atom_a,
            atom_b=atom_b,
        )

        for column in mapping_columns:
            result.iat[
                row_number, result.columns.get_loc(column)
            ] = mapped[column]

    return result