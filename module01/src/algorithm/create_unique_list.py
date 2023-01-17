import sys

import pandas as pd
import requests as r
from io import StringIO
import os


def create_list_of_unique_proteins(data, search_tool, intra_only, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                                   hhsearch_out):
    # create pandas dataframe of unique proteins depending on uniprot ids
    #
    # input data: pd.DataFrame, search_tool: str, intra_only: bool, blast_bin: str/None, blast_db: str,
    # hhsearch_bin: str/None, hhsearch_db: str, hhsearch_out: str
    # return unique_proteins_list: pd.DataFrame

    if intra_only:
        # Collect first rows in dataset of unique proteins
        unique_protein_rows = [data[data.unip_id == protein].iloc[0] for protein in data.unip_id.unique()]

        # Collect uniprot ids, indeces, sequences and counts of unique proteins
        unique_proteins = data.unip_id.unique().tolist()
        unique_protein_indeces = [row.name for row in unique_protein_rows]
        unique_sequences = [row.seq for row in unique_protein_rows]
        unique_protein_counts = [len(data[data.unip_id == protein].index) for protein in unique_proteins]
    else:
        # Collect first rows in dataset of unique proteins
        # print([data[(data.unip_id_a == protein)] for protein in (data.unip_id_a + data.unip_id_b).unique()])
        unique_protein_rows = [data[(data.unip_id_a == protein) | (data.unip_id_b == protein)].iloc[0]
                               for protein in pd.concat([data.unip_id_a, data.unip_id_b]).unique()]

        # Collect uniprot ids, indeces, sequences and counts of unique proteins
        unique_proteins = pd.concat([data.unip_id_a, data.unip_id_b]).unique().tolist()
        unique_protein_indeces = [row.name for row in unique_protein_rows]
        unique_sequences = [unique_protein_rows[i].seq_a
                            if unique_protein_rows[i].unip_id_a == protein else unique_protein_rows[i].seq_b
                            for i, protein in enumerate(unique_proteins)]
        unique_protein_counts = [len(data[(data.unip_id_a == protein) | (data.unip_id_b == protein)].index)
                                 for protein in unique_proteins]

    # Apply uniprot search for information on unique proteins
    infos = search_uniprot_metadata(unique_proteins)
    # Apply blastp or hhsearch search for pdb entries on unique proteins by their sequences
    pdbs = search_pdb_entries(unique_proteins, unique_sequences, search_tool, blast_bin, blast_db, hhsearch_bin,
                              hhsearch_db, hhsearch_out)

    # Collect information of unique proteins for final unique protein list dataset
    unique_proteins_list = pd.DataFrame()
    unique_proteins_list["Index"] = unique_protein_indeces
    for i, head in enumerate(infos[0][0].split('\t')):
        unique_proteins_list[head] = [info[1].split('\t')[i].replace('\"', '') for info in infos]
    unique_proteins_list["Sequence"] = unique_sequences
    unique_proteins_list["PDB"] = pdbs
    unique_proteins_list["Count"] = unique_protein_counts

    return unique_proteins_list


def search_uniprot_metadata(unique_proteins):
    # retrieve full sequences and info from uniprot database by given entries
    #
    # input unique_proteins: list(str)
    # return infos: list(list(str))

    # Create data container lists
    infos = []

    ind = 1
    # Iterate over proteins (proteins = uniprot ids)
    for protein in unique_proteins:
        print(f"\r\t[{round(ind * 100 / len(unique_proteins), 2)}%]", end='')
        ind += 1

        # Retrieve uniprot information on protein
        urllib = f"https://rest.uniprot.org/uniprotkb/search?query={protein}&format=tsv"
        try:
            info = r.get(urllib).text.split('\n')
        except ConnectionError as e:
            print("No connection to UniProt API possible. Please try again later.")
            print(e)
            sys.exit()

        # Add information and sequence to container lists
        infos.append(info)
    print()

    return infos


def search_pdb_entries(proteins, sequences, search_tool, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                       hhsearch_out):
    # use either hhsearch or blastp as search tool on protein sequence in order to retrieve matching pdb id, if no
    # result was found add an alphafold entry id instead (id: af<uniprot_id>_A)
    #
    # input proteins: list(str), sequences: list(str), search_tool: str, blast_bin: str/None, blast_db: str,
    # hhsearch_bin: str/None, hhsearch_db: str, hhsearch_out: str
    # return pdbs: list(str)

    # Create data container list
    pdbs = []

    ind = 1
    # Iterate over proteins (proteins = uniprot ids)
    for i, protein in enumerate(proteins):
        print(f"\r\t[{round(ind * 100 / len(proteins), 2)}%]", end='')
        ind += 1

        # Create temporary fasta file at data/temp/unique_protein_list for commandline application in search tools
        project_path = '/'.join(os.path.abspath(__file__).split('/')[:-4])
        project_path = project_path + '/' if project_path else ""
        temp_path = f"{project_path}data/temp/unique_protein_list/"
        with open(f"{temp_path}tmp.fasta", 'w') as tmp_file:
            tmp_file.write(f">Name\n{sequences[i]}\n")

        # Initialize result as False
        res = False

        # Depending on given string either perform blastp or hhsearch
        if search_tool == "blastp":
            blast_call = "blastp" if blast_bin is None else f"{blast_bin}blastp"
            cmd = f"{blast_call} -query {temp_path}tmp.fasta -db {blast_db}pdbaa -evalue 1e-5 -max_target_seqs 20 " \
                  f"-outfmt \"6 delim=, saccver pident qcovs evalue\""
            res = pd.read_csv(StringIO(os.popen(cmd).read()),
                              sep=',',
                              names=["pdb", "ident", "cov", "eval"],
                              dtype={"pdb": str, "ident": float, "cov": float, "eval": float})
            # Isolate search result to entries with identity of at least 90% and coverage of at least 50%
            res = res[(res["ident"] >= 90) & (res["cov"] >= 50)]
            # If no result remains reset result to False, else take first fitting entry from pdb column
            if res.empty:
                res = False
            else:
                res = res.iloc[0, :]["pdb"]
        elif search_tool == "hsearch":
            hhearch_call = "hhsearch" if hhsearch_bin is None else f"{hhsearch_bin}hhsearch"
            cmd = f"{hhearch_call} -i {temp_path}tmp.fasta -d {hhsearch_db}pdb70 -e 1e-5 -qid 90 -cov 50 -blasttab " \
                  f"{hhsearch_out} -v 0 -cpu 20"
            os.system(cmd)
            # Open hhsearch output (Note: hhsearch outs cannot be retrieved from the commandline, as it is the case with
            # blastp)
            res = open(f"{temp_path}tmp.hhr", 'r').read().split('\n')[9][4:10]
        # If result is not False, append it to container list, else append alphafold entry id instead
        if res:
            pdbs.append(res)
        else:
            pdbs.append(f"af{protein}_A")
    print()

    return pdbs
