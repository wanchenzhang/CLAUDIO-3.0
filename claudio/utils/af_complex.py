"""Retrieve and inspect existing AlphaFold protein-complex models."""

from Bio.PDB import MMCIFParser
from Bio.SeqUtils import seq1
import requests
from pathlib import Path
import gzip
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

def fetch_complex_metadata(model_id):
    """Retrieve AFDB metadata for an existing protein-complex model."""
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{model_id}"

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    records = response.json()
    if records == []:
        raise ValueError(f"No metadata found for model {model_id}")

    for record in records:
        if record.get("isComplex") is not True:
            raise ValueError(f"Model {model_id} is not a complex")

    return records
def download_complex_cif(records, output_dir):
    """Download a complex mmCIF file and return its local path."""
    cif_urls = {record["cifUrl"] for record in records}
    if len(cif_urls) != 1:
        raise ValueError("Inconsistent CIF URLs in complex metadata")

    model_identities = {
        (record["modelEntityId"], record["latestVersion"])
        for record in records
    }
    if len(model_identities) != 1:
        raise ValueError("Inconsistent model identity in complex metadata")
    record = records[0]
    url = record["cifUrl"]

    response = requests.get(url, timeout=30)
    response.raise_for_status()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = (
        f"{record['modelEntityId']}"
        f"-model_v{record['latestVersion']}.cif"
    )
    structure_path = output_dir / filename
    content = response.content
    if content.startswith(b"\x1f\x8b"):
        content = gzip.decompress(content)

    structure_path.write_bytes(content)

    return structure_path

def fetch_complex_candidates(uniprot_id):
    """Return AFDB complex candidates containing a UniProt accession."""
    url = f"https://alphafold.ebi.ac.uk/api/complex/{uniprot_id}"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    return response.json()
