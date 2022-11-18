import sys

import pandas as pd


def read_temp_file(data, search_tool, intra_only, temp_path):
    # read temporary save file containing earlier hhsearch or blastp search results and concatenate them to the dataset,
    # for quick reruns
    #
    # input data: pd.DataFrame, search_tool: str, intra_only: bool, temp_path: str
    # return read_temp_file: pd.DataFrame

    dataset, filename = data
    tmp_filepath = f"{temp_path}{'.'.join(filename.split('.')[:-1])}_{search_tool}_bltmp.{filename.split('.')[-1]}"
    tmp_data = pd.read_csv(tmp_filepath)

    dataset["pdb_id"] = tmp_data["pdb_id"]
    dataset["chain"] = tmp_data["chain"]
    dataset["all_results"] = tmp_data["all_results"]
    if not intra_only:
        dataset["pdb_id_b"] = tmp_data["pdb_id_b"]
        dataset["chain_b"] = tmp_data["chain_b"]
        dataset["all_results_b"] = tmp_data["all_results_b"]
    return dataset.fillna('')
