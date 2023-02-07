import pandas as pd


def read_in(input_filepath):
    # read inputfile (.sqcs file), and columns for results and for control
    #
    # input input_filepath: str
    # return data: pd.DataFrame, filename: str, intra_only: bool

    # Read input
    data = pd.read_csv(input_filepath, index_col=0)
    filename = input_filepath.split('/')[-1]
    intra_only = all([x not in data.columns for x in ["unip_id_a", "unip_id_b", "seq_a", "seq_b"]])

    # Add result and control columns
    data["all_results"] = ''
    data[["pdb_id", "path", "pdb_method", "pdb_resolution", "best_res_pdb_method", "best_res_pdb_resolution"]] = '-'
    if intra_only:
        data["chain"] = '-'
    else:
        data[["chain_a", "chain_b"]] = '-'

    return data, filename, intra_only
