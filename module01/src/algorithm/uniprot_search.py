import sys
import os

import requests as r
import pandas as pd


def do_uniprot_search(data, intra_only, filename):
    # Retrieve full uniprot sequences and IDs if not given
    #
    # input data: pd.DataFrame, intra_only: bool, filename: str
    # return data: pd.DataFrame

    # retrieve sequences from uniprot entries
    _, seqs = search_uniprot(data, intra_only)

    data["seq_a"] = seqs[0]

    # if only intra type interactions copy site a results to site b, else get other results
    data["seq_b"] = seqs[0] if intra_only else seqs[1]

    # save results in temporary save file (can be used on rerun, instead of searching results again)
    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-4])
    project_path = project_path + '/' if project_path else ""
    temp_save_filepath = f"{project_path}data/temp/uniprot_search/" \
                         f"{'.'.join(filename.split('.')[:-1])}_srtmp.{filename.split('.')[-1]}"
    data[["seq_a", "seq_b"]].to_csv(temp_save_filepath, index=False)

    return data


def search_uniprot(data, intra_only):
    # search uniprot database for sequences
    #
    # input data: pd.DataFrame, intra_only: bool
    # return ids: list(list(str)), seqs: list(list(str))

    ids, \
        seqs = ([[], []] for _ in range(2))

    ids[0] = list(data.unip_id_a)
    ids[1] = list(data.unip_id_b)

    # create list of unique uniprot ids
    unip_ids = list(pd.concat([data.unip_id_a, data.unip_id_b]).unique())

    # retrieve all possible uniprot sequences
    unip_search_results = []
    for id in unip_ids:
        url = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query={id}"
        try:
            result = [''.join(x.split('\n')[1:]) for x in r.get(url).text.split('>') if x]
        except ConnectionError as e:
            print("No connection to UniProt API possible. Please try again later.")
            print(e)
            sys.exit()
        unip_search_results.append(result)

    ind = 0
    for _, row in data.iterrows():
        id = row["unip_id_a"]
        ind += 1
        print(f"\r\t[{round(ind * 100 / len(data.index), 2) if intra_only else round((ind * 100 / len(data.index)) / 2, 2)}%]", end='')
        if pd.isna(id):
            seqs[0].append(float('nan'))
            if intra_only:
                seqs[1].append(float('nan'))
        else:
            result = unip_search_results[unip_ids.index(id)]
            fitting_seq_found = False
            seq = ''
            if len(result) > 1:
                for seq in result:
                    peptide_arg = ((row["pep_a"] in seq) and (row["pep_b"] in seq)) \
                        if intra_only else (row["pep_a"] in seq)
                    if peptide_arg:
                        fitting_seq_found = True
                        break
            if not fitting_seq_found:
                seq = result[0]
            seqs[0].append(seq)

    # Do separate search for site b, if not intra types only
    if not intra_only:
        ind = 0
        for i, row in data.iterrows():
            ind += 1
            print(f"\r\t[{round((ind * 100 / len(data.index)) / 2, 2) + 50.0:.1f}%]", end='')
            # If XL-type inter do separate search, else copy results from site a
            if row.XL_type != "intra":
                id = row["unip_id_b"]
                if pd.isna(id):
                    seqs[1].append(float('nan'))
                else:
                    result = unip_search_results[unip_ids.index(id)]
                    fitting_seq_found = False
                    seq = ''
                    if len(result) > 1:
                        for seq in result:
                            if row["pep_b"] in seq:
                                fitting_seq_found = True
                                break
                    if not fitting_seq_found:
                        seq = result[0]
                    seqs[1].append(seq)
            else:
                seqs[1].append(seqs[0][i])
    print()

    return ids, seqs