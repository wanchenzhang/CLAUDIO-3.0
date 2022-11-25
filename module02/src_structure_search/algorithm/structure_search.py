import os
from io import StringIO
import pandas as pd


def structure_search(data, search_tool, e_value, query_id, coverage, intra_only, temp_path, blast_bin, blast_db,
                     hhsearch_bin, hhsearch_db, hhsearch_out):
    # Perform either hhsearch or blastp search for any unique uniprot entry in rcsb database and retrieve the best
    # results, and save all entries for which no sufficient result was returned to retrieve these from alphafold
    # database later
    #
    # input data: pd.DataFrame, search_tool: str, e_value: float, query_id: float, coverage: float, intra_only: bool,
    # temp_path: str, blast_bin: str/None, blast_db: str, hhsearch_bin: str/None, hhsearch_db: str, hhsearch_out: str
    # return dataset: pd.DataFrame

    dataset, filename = data

    # Container lists for search results
    not_found = []
    ind = 0

    # Save each already search uniprot entries' id and result for quick retrieval if reencountered (instead of repeating
    # the same search)
    already_searched_ids, \
        already_searched_res = ([] for _ in range(2))

    for i, row in dataset.iterrows():
        # If encountered a new uniprot entry perform actual search
        if row["unip_id_a"] not in already_searched_ids:
            # Save uniprot sequence in a temporary fasta file for search tool commandline call
            # (override before each new search)
            with open(f"{temp_path}tmp.fasta", 'w') as tmp_file:
                tmp_file.write(f">Name\n{row['seq_a']}\n")
                tmp_file.close()
                best_result_a = []

                # Perform search with either hhsearch or blastp (Note: Watch environmental variables $BLASTDB, $HHDB,
                # $HHOUT to be set according to instructions found in README.md)
                if search_tool == "blastp":
                    blast_call = "blastp" if blast_bin is None else f"{blast_bin}blastp"
                    command_a = f"{blast_call} -query {temp_path}tmp.fasta -db {blast_db}pdbaa -evalue {e_value} " \
                                f"-max_target_seqs 20 -outfmt \"6 delim=, saccver pident qcovs evalue\""
                    res_a = pd.read_csv(StringIO(os.popen(command_a).read()),
                                        sep=',',
                                        names=["pdb", "ident", "cov", "eval"],
                                        dtype={"pdb": str, "ident": float, "cov": float, "eval": float})
                    best_result_a = res_a[(res_a["ident"] >= query_id) & (res_a["cov"] >= coverage)]
                    best_result_a = best_result_a.loc[:, "pdb"].tolist()
                elif search_tool == "hhsearch":
                    hhsearch_call = "hhsearch" if hhsearch_bin is None else f"{hhsearch_bin}hhsearch"
                    command_a = f"{hhsearch_call} -i {temp_path}tmp.fasta -d {hhsearch_db}pdb70 -e {e_value} -qid " \
                                f"{query_id} -cov {coverage} -blasttab {hhsearch_out} -v 0 -cpu 20"
                    os.system(command_a)
                    best_result_a = []
                    for j, line in enumerate(open(f"{temp_path}tmp.hhr", 'r').read().split('\n')):
                        if j > 8:
                            best_result_a.append(line[4:10])

                # Save uniprot id to already_search_ids
                already_searched_ids.append(row["unip_id_a"])
                # If search successful, save result to respective container lists, and full result to
                # already_searched_res
                if best_result_a:
                    # If '_' in result means a chain is specified in the result
                    if '_' in best_result_a[0]:
                        dataset.loc[i, "pdb_id"] = best_result_a[0].split('_')[0]
                        dataset.loc[i, "chain"] = best_result_a[0].split('_')[1]
                        already_searched_res.append((best_result_a[0].split('_')[0], best_result_a[0].split('_')[1],
                                                     ' '.join(best_result_a)))

                    # Else save '-' inplace instead (means the full pdb (not just a specific chain) will be reviewed
                    # later)
                    else:
                        dataset.loc[i, "pdb_id"] = best_result_a[0][:4]
                        dataset.loc[i, "chain"] = '-'
                        already_searched_res.append((best_result_a[0][:4], '-', ' '.join(best_result_a)))

                    dataset.loc[i, "all_results"] = ' '.join(best_result_a)

                # If search was unsuccessful, save None to respective container lists and already_searched_res;
                # furthermore save entry index to not_found (will be retrieved from alphafold, if possible)
                else:
                    dataset.loc[i, "pdb_id"] = None
                    dataset.loc[i, "chain"] = None
                    already_searched_res.append((None, None, ''))
                    not_found.append(i)

        # If reencountered an uniprot id, copy result from already_searched_res
        else:
            id_a, chain_a, best_result_a = already_searched_res[already_searched_ids.index(row["unip_id_a"])]
            dataset.loc[i, "pdb_id"] = id_a
            dataset.loc[i, "chain"] = chain_a
            dataset.loc[i, "all_results"] = best_result_a

        # If intra_only False, perform same actions independently for site b
        if not intra_only:

            # If encountered a new uniprot entry perform actual search
            if row["unip_id_b"] not in already_searched_ids:

                # Save uniprot sequence in a temporary fasta file for search tool commandline call
                # (override before each new search)
                with open(f"{temp_path}tmp.fasta", 'w') as tmp_file:
                    tmp_file.write(f">Name\n{row['seq_b']}\n")
                    tmp_file.close()
                    best_result_b = []

                    # Perform search with either hhsearch or blastp (Note: Watch environmental variables $BLASTDB,
                    # $HHDB, $HHOUT to be set according to instructions found in README.md)
                    if search_tool == "blastp":
                        command_b = f"{blast_call} -query {temp_path}tmp.fasta -db {blast_db}pdbaa -evalue {e_value} " \
                                    f"-max_target_seqs 20 -outfmt \"6 delim=, saccver pident qcovs evalue\""
                        res_b = pd.read_csv(StringIO(os.popen(command_b).read()),
                                            sep=',',
                                            names=["pdb", "ident", "cov", "eval"],
                                            dtype={"pdb": str, "ident": float, "cov": float, "eval": float})
                        best_result_b = res_b[(res_b["ident"] >= query_id) & (res_b["cov"] >= coverage)]
                        best_result_b = best_result_b.loc[:, "pdb"].tolist()
                    elif search_tool == "hhsearch":
                        command_b = f"{hhsearch_call} -i {temp_path}tmp.fasta -d {hhsearch_db}pdb70 -e {e_value} " \
                                    f"-qid {query_id} -cov {coverage} -blasttab {hhsearch_out} -v 0 -cpu 20"
                        os.system(command_b)
                        best_result_b = []
                        for j, line in enumerate(open(f"{temp_path}tmp.hhr", 'r').read().split('\n')):
                            if j > 8:
                                best_result_b.append(line[4:10])

                    # Save uniprot id to already_search_ids
                    already_searched_ids.append(row["unip_id_b"])
                    # If search successful, save result to respective container lists, and full result to
                    # already_searched_res
                    if best_result_b:
                        # If '_' in result means a chain is specified in the result
                        if '_' in best_result_b[0]:
                            dataset.loc[i, "pdb_id_b"] = best_result_b[0].split('_')[0]
                            dataset.loc[i, "chain_b"] = best_result_b[0].split('_')[1]
                            already_searched_res.append((best_result_b[0].split('_')[0], best_result_b[0].split('_')[1],
                                                         ' '.join(best_result_b)))

                        # Else save '-' inplace instead (means the full pdb (not just a specific chain) will be reviewed
                        # later)
                        else:
                            dataset.loc[i, "pdb_id_b"] = best_result_b[0][:4]
                            dataset.loc[i, "chain_b"] = '-'
                            already_searched_res.append((best_result_b[0][:4], '-', ' '.join(best_result_b)))

                        dataset.loc[i, "all_results_b"] = ' '.join(best_result_b)

                    # If search was unsuccessful, save None to respective container lists and already_searched_res;
                    # furthermore save entry index to not_found (will be retrieved from alphafold, if possible)
                    else:
                        dataset.loc[i, "pdb_id_b"] = None
                        dataset.loc[i, "chain_b"] = None
                        not_found.append(i)
                        already_searched_res.append((None, None, ''))

            # If reencountered an uniprot id, copy result from already_searched_res
            else:
                id_b, chain_b, best_result_b = already_searched_res[already_searched_ids.index(row["unip_id_b"])]
                dataset.loc[i, "pdb_id_b"] = id_b
                dataset.loc[i, "chain_b"] = chain_b
                dataset.loc[i, "all_results_b"] = best_result_b

        ind += 1
        print(f"\r\t[{round(ind * 100 / len(dataset.index), 2)}%]", end='')
    print()

    # Save results to temporary save file
    marker_string = f"_{search_tool}_bltmp."
    temp_save_filepath = f"{temp_path}{'.'.join(filename.split('.')[:-1])}{marker_string}{filename.split('.')[-1]}"
    if intra_only:
        dataset[["pdb_id", "chain", "all_results"]].to_csv(temp_save_filepath, index=False)
    else:
        dataset[["pdb_id", "pdb_id_b", "chain", "chain_b", "all_results", "all_results_b"]].to_csv(temp_save_filepath, index=False)

    # Print ids of entries which were not found in rcsb database (will be retrieved from alphafold database instead)
    not_found_proteins = dataset.iloc[not_found].unip_id_a.unique()
    print(f"\tProteins which yielded no results from RCSB database (will be retrieved from AlphaFold "
          f"(n = {len(not_found_proteins)})): {not_found_proteins}")

    return dataset
