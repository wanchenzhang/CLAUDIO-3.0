import pandas as pd


def read_temp_file(data, filename, search_tool, temp_path):
    # read temporary save file containing earlier hhsearch or blastp search results and concatenate them to the dataset,
    # for quick reruns
    #
    # input data: pd.DataFrame, filename: str, search_tool: str, temp_path: str
    # return read_temp_file: pd.DataFrame

    tmp_filepath = f"{temp_path}{'.'.join(filename.split('.')[:-1])}_{search_tool}_bltmp.{filename.split('.')[-1]}"
    tmp_data = pd.read_csv(tmp_filepath)

    data["pdb_id"] = tmp_data["pdb_id"]
    data["all_results"] = tmp_data["all_results"]
    data["chain_a"] = tmp_data["chain_a"]
    data["chain_b"] = tmp_data["chain_b"]

    return data.fillna('')
