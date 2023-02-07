import os
from io import StringIO
import pandas as pd
import numpy as np


def structure_search(data, filename, search_tool, e_value, query_id, coverage, intra_only, temp_path, blast_bin,
                     blast_db, hhsearch_bin, hhsearch_db, hhsearch_out):
    # Perform either hhsearch or blastp search for any unique uniprot entry in rcsb database and retrieve the best
    # results, and save all entries for which no sufficient result was returned to retrieve these from alphafold
    # database later
    #
    # input data: pd.DataFrame, filename: str, search_tool: str, e_value: float, query_id: float, coverage: float,
    # intra_only: bool, temp_path: str, blast_bin: str/None, blast_db: str, hhsearch_bin: str/None, hhsearch_db: str,
    # hhsearch_out: str
    # return dataset: pd.DataFrame

    # Container lists for search results
    not_found = []
    ind = 0

    # Save each already search uniprot entries' id and result for quick retrieval if reencountered (instead of repeating
    # the same search)
    already_searched = {}

    for i, row in data.iterrows():
        # Perform single uniprot id search for intra dataset
        if intra_only:
            # If encountered a new uniprot entry, do search and save results
            if row["unip_id"] not in already_searched.keys():
                best_results, pdb_id, chain = perform_search(row, '', search_tool, e_value, query_id, coverage,
                                                             temp_path, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                                                             hhsearch_out)
                already_searched[row["unip_id"]] = (pdb_id, chain, ' '.join(best_results))

            if already_searched[row["unip_id"]][2]:
                # Sort results into dataset
                pdb_id, chain, best_results = already_searched[row["unip_id"]]
                data.loc[i, "pdb_id"] = pdb_id
                data.loc[i, "chain"] = chain
                data.loc[i, "all_results"] = best_results
            else:
                not_found.append(i)

        # Else perform separate search for uniprot_id_a and uniprot_id_b
        else:
            # If encountered a new uniprot entry at site_a, do search and save results
            if row["unip_id_a"] not in already_searched.keys():
                best_results, pdb_id, chain = perform_search(row, 'a', search_tool, e_value, query_id, coverage,
                                                             temp_path, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                                                             hhsearch_out)
                already_searched[row["unip_id_a"]] = (pdb_id, chain, ' '.join(best_results))
            # If encountered a new uniprot entry at site_b, do search and save results
            if row["unip_id_b"] not in already_searched.keys():
                best_results, pdb_id, chain = perform_search(row, 'b', search_tool, e_value, query_id, coverage,
                                                             temp_path, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                                                             hhsearch_out)
                already_searched[row["unip_id_b"]] = (pdb_id, chain, ' '.join(best_results))

            if already_searched[row["unip_id_a"]][2] and already_searched[row["unip_id_b"]][2]:
                # Sort results into dataset
                _, _, best_results_a = already_searched[row["unip_id_a"]]
                _, _, best_results_b = already_searched[row["unip_id_b"]]

                # Join search results
                results_a = {result.split('_')[0]: result.split('_')[1] for result in best_results_a.split(' ')}
                results_b = {result.split('_')[0]: result.split('_')[1] for result in best_results_b.split(' ')}
                intersect_results = {key: (results_a[key], results_b[key])
                                     for key in [key for key in results_a.keys() if key in results_b.keys()]}

                # Add results to dataset
                if intersect_results:
                    pdb_id = list(intersect_results.keys())[0]
                    chain_a, chain_b = intersect_results[pdb_id]
                    data.loc[i, "pdb_id"] = pdb_id
                    data.loc[i, "chain_a"] = chain_a
                    data.loc[i, "chain_b"] = chain_b
                    data.loc[i, "all_results"] = ' '.join([f"{key}_{value[0]}_{value[1]}"
                                                           for key, value in intersect_results.items()])
            else:
                not_found.append(i)
        ind += 1
        print(f"\r\t[{round(ind * 100 / len(data.index), 2)}%]", end='')
    print()

    # Save results to temporary save file
    marker_string = f"_{search_tool}_bltmp."
    temp_save_filepath = f"{temp_path}{'.'.join(filename.split('.')[:-1])}{marker_string}{filename.split('.')[-1]}"
    if intra_only:
        data[["pdb_id", "chain", "all_results"]].to_csv(temp_save_filepath, index=False)
    else:
        data[["pdb_id", "chain_a", "chain_b", "all_results"]].to_csv(temp_save_filepath, index=False)

    # Print ids of entries which were not found in rcsb database (will be retrieved from alphafold database instead)
    if intra_only:
        not_found_proteins = data.loc[not_found].unip_id.unique()
        print(f"\tProteins which yielded no results from RCSB database (will be retrieved from AlphaFold "
              f"(n = {len(not_found_proteins)})): {not_found_proteins}")
    elif not data[data.unip_id_a == data.unip_id_b].empty:
        not_found_proteins = pd.concat(data.loc[not_found].unip_id_a, data.loc[not_found].unip_id_b).unique()
        print(f"\tProteins which yielded no results from RCSB database (will be retrieved from AlphaFold "
              f"(n = {len(not_found_proteins)})): {not_found_proteins}")

    return data


def perform_search(data, site, search_tool, e_value, query_id, coverage, temp_path, blast_bin, blast_db, hhsearch_bin,
                   hhsearch_db, hhsearch_out):
    # Perform either hhsearch or blastp search for unique uniprot entry in rcsb database and return possible results
    #
    # input data: pd.Series, site: str, search_tool: str, e_value: float, query_id: float, coverage: float,
    # temp_path: str, blast_bin: str/None, blast_db: str, hhsearch_bin: str/None, hhsearch_db: str, hhsearch_out: str
    # return best_result: list(str), pdb_id: str/None, chain: str/None

    # Save uniprot sequence in a temporary fasta file for search tool commandline call
    # (override before each new search)
    with open(f"{temp_path}tmp.fasta", 'w') as tmp_file:
        tmp_file.write(f">Name\n{data[f'seq_{site}' if site else 'seq']}\n")
        tmp_file.close()
        best_result = []

        # Perform search with either hhsearch or blastp (Note: Watch environmental variables $BLASTDB,
        # $HHDB, $HHOUT to be set according to instructions found in README.md)
        if search_tool == "blastp":
            blast_call = "blastp" if blast_bin is None else f"{blast_bin}blastp"
            command = f"{blast_call} -query {temp_path}tmp.fasta -db {blast_db}pdbaa -evalue {e_value} " \
                      f"-max_target_seqs 20 -outfmt \"6 delim=, saccver pident qcovs evalue\""
            res = pd.read_csv(StringIO(os.popen(command).read()), sep=',', names=["pdb", "ident", "cov", "eval"],
                              dtype={"pdb": str, "ident": float, "cov": float, "eval": float})
            best_result = res[(res["ident"] >= query_id) & (res["cov"] >= coverage)]
            best_result = best_result.loc[:, "pdb"].tolist()
        elif search_tool == "hhsearch":
            hhsearch_call = "hhsearch" if hhsearch_bin is None else f"{hhsearch_bin}hhsearch"
            command = f"{hhsearch_call} -i {temp_path}tmp.fasta -d {hhsearch_db}pdb70 -e {e_value} -qid" \
                      f" {query_id} -cov {coverage} -blasttab {hhsearch_out} -v 0 -cpu 20"
            os.system(command)
            best_result = []
            for j, line in enumerate(open(f"{temp_path}tmp.hhr", 'r').read().split('\n')):
                if j > 8:
                    best_result.append(line[4:10])
        # If search successful, save result to respective container lists, and full result to
        # already_searched_res with specified chain
        if best_result and ('_' in best_result[0]):
            pdb_id = best_result[0].split('_')[0]
            chain = best_result[0].split('_')[1]
        # Elif no chain specified
        elif best_result:
            pdb_id = best_result[0][:4]
            chain = '-'
        # Else, search was unsuccessful
        else:
            pdb_id = None
            chain = None

    return best_result, pdb_id, chain
