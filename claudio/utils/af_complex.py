"""Retrieve and inspect existing AlphaFold protein-complex models."""

from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1


def read_complex_chains(structure_path, records):
    """Read protein chains and associate them with AFDB metadata."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("af_complex", str(structure_path))
    model = next(structure.get_models())

    metadata_by_chain = {
        record["chainId"]: record
        for record in records
    }

    chains = []

    for chain in model:
        residues = [
            residue for residue in chain
            if residue.id[0] == " "
        ]

        sequence = "".join(
            seq1(residue.resname)
            for residue in residues
        )
        
        if chain.id not in metadata_by_chain:
            raise ValueError(
                f"Missing metadata for chain {chain.id}"
            )

        metadata = metadata_by_chain[chain.id]
        if sequence != metadata["sequence"]:
            raise ValueError(
                f"Sequence mismatch for chain {chain.id}: "
                "structure sequence differs from AFDB metadata"
            )

        chains.append({
            "chain_id": chain.id,
            "protein_id": metadata["uniprotAccession"],
            "residue_count": len(residues),
            "sequence": sequence,
        })

    return chains
