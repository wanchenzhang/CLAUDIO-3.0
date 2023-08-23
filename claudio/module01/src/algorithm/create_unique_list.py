import socket
import sys

import requests as r
from io import StringIO

from utils.utils import *


def create_list_of_unique_proteins(data, unique_protein_temp_dir, search_tool, blast_bin, blast_db, hhsearch_bin,
                                   hhsearch_db, verbose_level):
    # create pandas dataframe of unique proteins depending on uniprot ids
    #
    # input data: pd.DataFrame, unique_protein_temp_dir: str, search_tool: str, blast_bin: str/None, blast_db: str,
    # hhsearch_bin: str/None, hhsearch_db: str, verbose_level: int
    # return unique_proteins_list: pd.DataFrame

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
    infos = search_uniprot_metadata(unique_proteins, verbose_level)
    # Apply blastp or hhsearch search for pdb entries on unique proteins by their sequences
    pdbs = search_pdb_entries(unique_proteins, unique_sequences, unique_protein_temp_dir, search_tool, blast_bin,
                              blast_db, hhsearch_bin, hhsearch_db, verbose_level)

    # Collect information of unique proteins for final unique protein list dataset
    unique_proteins_list = pd.DataFrame()
    unique_proteins_list["Index"] = unique_protein_indeces
    for i, head in enumerate(infos[0][0].split('\t')):
        try:
            unique_proteins_list[head] = [info[1].split('\t')[i].replace('\"', '') for info in infos]
        except:
            print(head)
            for info in infos:
                try:
                    _ = info[1].split('\t')[i].replace('\"', '')
                except:
                    print(info)
                    sys.exit()
    unique_proteins_list["Sequence"] = unique_sequences
    unique_proteins_list["PDB"] = pdbs
    unique_proteins_list["Count"] = unique_protein_counts

    return unique_proteins_list


def search_uniprot_metadata(unique_proteins, verbose_level):
    # retrieve full sequences and info from uniprot database by given entries
    #
    # input unique_proteins: list(str), verbose_level: int
    # return infos: list(list(str))

    # Create data container lists
    infos = []

    ind = 1
    # Iterate over proteins (proteins = uniprot ids)
    for protein in unique_proteins:
        verbose_print(f"\r\tMetadata search:[{round_self(ind * 100 / len(unique_proteins), 2)}%]", 1, verbose_level,
                      end='')
        ind += 1

        # Retrieve uniprot information on protein
        urllib = f"https://rest.uniprot.org/uniprotkb/search?query={protein}&format=tsv"
        try:
            info = r.get(urllib).text.split('\n')
        except (r.exceptions.Timeout, ConnectionError, socket.gaierror, r.exceptions.ConnectionError) as e:
            print("No connection to UniProt API possible. Please try again later.")
            print(e)
            sys.exit()

        # Add information and sequence to container lists
        infos.append(info)
    verbose_print("", 1, verbose_level)

    return infos


def search_pdb_entries(proteins, sequences, unique_protein_temp_dir, search_tool, blast_bin, blast_db, hhsearch_bin,
                       hhsearch_db, verbose_level):
    # use either hhsearch or blastp as search tool on protein sequence in order to retrieve matching pdb id, if no
    # result was found add an alphafold entry id instead (id: af<uniprot_id>_A)
    #
    # input proteins: list(str), sequences: list(str), unique_protein_temp_dir: str, search_tool: str,
    # blast_bin: str/None, blast_db: str, hhsearch_bin: str/None, hhsearch_db: str, verbose_level: int
    # return pdbs: list(str)

    # Create data container list
    pdbs = []

    ind = 1
    # Iterate over proteins (proteins = uniprot ids)
    for i, protein in enumerate(proteins):
        verbose_print(f"\r\tStructure search:[{round_self(ind * 100 / len(proteins), 2)}%]", 1, verbose_level, end='')

        # Create temporary fasta file at data/temp/unique_protein_list for commandline application in search tools
        with open(f"{unique_protein_temp_dir}tmp{i}.fasta", 'w') as tmp_file:
            tmp_file.write(f">{protein}\n{sequences[i]}\n")

        # Initialize result as False
        res = False

        # Depending on given string either perform blastp or hhsearch
        if search_tool == "blastp":
            blast_call = "blastp" if blast_bin is None else f"{blast_bin}blastp"
            cmd = f"{blast_call} -query {unique_protein_temp_dir}tmp{i}.fasta -db {blast_db}pdbaa -evalue 1e-5 " \
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
        elif search_tool == "hhsearch":
            hhearch_call = "hhsearch" if hhsearch_bin is None else f"{hhsearch_bin}hhsearch"
            cmd = f"{hhearch_call} -i {unique_protein_temp_dir}tmp{i}.fasta -d {hhsearch_db}pdb70 -e 1e-5 " \
                  f"-blasttab {unique_protein_temp_dir}tmp{i}.hhr -qid 90 -cov 50 -v 0 -cpu 20"
            os.system(cmd)
            # Open hhsearch output (Note: hhsearch outs cannot be retrieved from the commandline, as it is the case with
            # blastp)
            res = [line.split('\t')[1] for line in open(f"{unique_protein_temp_dir}tmp{i}.hhr", 'r').read().split('\n')
                   if line.split('\t')[0] == protein][0]
        # If result is not False, append it to container list, else append alphafold entry id instead
        if res:
            pdbs.append(res)
        else:
            pdbs.append(f"af{protein}_A")
        ind += 1
    verbose_print("", 1, verbose_level)

    return pdbs
