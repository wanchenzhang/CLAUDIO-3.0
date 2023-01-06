import sys
import os

import requests as r
import pandas as pd


def do_uniprot_search(data, filename):
    # Retrieve full uniprot sequences and IDs if not given
    #
    # input data: pd.DataFrame, filename: str
    # return data: pd.DataFrame

    # retrieve sequences from uniprot entries
    data["seq"] = search_uniprot(data)

    # save results in temporary save file (can be used on rerun, instead of searching results again)
    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-4])
    project_path = project_path + '/' if project_path else ""
    temp_save_filepath = f"{project_path}data/temp/uniprot_search/" \
                         f"{'.'.join(filename.split('.')[:-1])}_srtmp.{filename.split('.')[-1]}"
    data[["seq"]].to_csv(temp_save_filepath, index=False)

    return data


def search_uniprot(data):
    # search uniprot database for sequences
    #
    # input data: pd.DataFrame
    # return seqs: list(str)

    seqs = []

    # create list of unique uniprot ids
    unip_ids = data.unip_id_a.unique().tolist()

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
        id = row.unip_id_a
        ind += 1
        print(f"\r\t[{round(ind * 100 / len(data.index), 2)}%]", end='')
        if pd.isna(id):
            seqs.append(float('nan'))
        else:
            result = unip_search_results[unip_ids.index(id)]
            fitting_seq_found = False
            seq = ''
            if len(result) > 1:
                for seq in result:
                    peptide_arg = (row["pep_a"] in seq) and (row["pep_b"] in seq)
                    if peptide_arg:
                        fitting_seq_found = True
                        break
            if not fitting_seq_found:
                seq = result[0]
            seqs.append(seq)
    print()

    return seqs
