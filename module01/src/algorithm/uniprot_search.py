import socket
import sys

import requests as r

from utils.utils import *


def do_uniprot_search(data, uniprot_search_temp_dir, filename, verbose_level):
    # Retrieve full uniprot sequences and IDs if not given
    #
    # input data: pd.DataFrame, uniprot_search_temp_dir: str, filename: str, verbose_level: int
    # return data: pd.DataFrame

    # retrieve sequences from uniprot entries
    data["seq_a"], search_result_dict = search_uniprot(data, verbose_level, site='a')
    data["seq_b"], _ = search_uniprot(data, verbose_level, already_searched=search_result_dict, site='b')

    # save results in temporary save file (can be used on rerun, instead of searching results again)
    temp_save_filepath = f"{uniprot_search_temp_dir}{filename}_srtmp.csv"
    data[["seq_a", "seq_b"]].to_csv(temp_save_filepath, index=False)

    return data


def search_uniprot(data, verbose_level, already_searched={}, site='a'):
    # search uniprot database for sequences
    #
    # input data: pd.DataFrame, verbose_level: int, already_searched: dict{str: list(str)}, site: str
    # return seqs: list(str), already_searched: dict{str: list(str)}

    seqs = []

    # create list of unique uniprot ids
    unip_ids = data[f"unip_id_{site}"].unique().tolist()

    # retrieve all unique uniprot sequences
    for id in unip_ids:
        # if search for unip id has not been performed yet, do so, and add it to already_searched dictionary
        if id not in already_searched.keys() or already_searched[id] is None:
            url = f"https://rest.uniprot.org/uniprotkb/search?format=fasta&query={id}"
            try:
                url_return_text = r.get(url).text
                return_failed = "Error encountered when streaming data. Please try again later." in url_return_text
                # if successful continue
                if not return_failed:
                    result = [''.join(x.split('\n')[1:]) for x in url_return_text.split('>') if x]
                    already_searched[id] = result
                # else print error message and raise ValueError
                else:
                    verbose_print(f"\tWarning! UniProt API call failed for UniProt_ID={id}.\n\tReturned message: "
                                  f"{url_return_text}", 2, verbose_level)
                    already_searched[id] = None
            except (r.exceptions.Timeout, ConnectionError, socket.gaierror, r.exceptions.ConnectionError) as e:
                print("No connection to UniProt API possible. Please try again later.")
                print(e)
                sys.exit()
            except ValueError:
                print("Error! Encountered at least one faulty return from the UniProt database.")
                sys.exit()

    ind = 0
    for _, row in data.iterrows():
        id = row[f"unip_id_{site}"]
        ind += 1
        verbose_print(f"\r\tSite_{site}:[{round_self(ind * 100 / len(data.index), 2)}%]", 1, verbose_level, end='')
        if pd.isna(id) or already_searched[id] is None or not already_searched[id]:
            seqs.append(float('nan'))
        else:
            result = already_searched[id]
            fitting_seq_found = False
            seq = ''
            if len(result) > 1:
                for seq in result:
                    peptide_arg = (row["pep_a"] in seq) and (row["pep_b"] in seq) \
                        if row["unip_id_a"] == row["unip_id_b"] else (row[f"pep_{site}"] in seq)
                    if peptide_arg:
                        fitting_seq_found = True
                        break
            if not fitting_seq_found:
                seq = result[0]
            seqs.append(seq)
    verbose_print("", 1, verbose_level)

    return seqs, already_searched
