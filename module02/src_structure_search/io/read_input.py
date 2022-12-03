import pandas as pd


def read_in(input_filepath):
    # read inputfile (.sqcs file), and columns for results and for control
    #
    # input input_filepath: str
    # return data: pd.DataFrame, filename: str

    # Read input
    data = pd.read_csv(input_filepath, index_col=0)
    filename = input_filepath.split('/')[-1]

    # Add result and control columns
    data["all_results"] = ''
    data[["pdb_id", "chain", "path", "pdb_method", "pdb_resolution", "best_res_pdb_method",
          "best_res_pdb_resolution"]] = '-'

    return data, filename
