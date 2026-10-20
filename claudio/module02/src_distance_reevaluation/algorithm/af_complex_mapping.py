"""Map protein positions to residues in existing AF complex models."""

from Bio.PDB import MMCIFParser

from claudio.utils.af_complex import read_complex_chains


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

    chains = read_complex_chains(structure_path, records)

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