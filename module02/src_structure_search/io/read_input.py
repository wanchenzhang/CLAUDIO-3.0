import pandas as pd


def read_in(input_filepath):
    # read inputfile (.sqcs file), and columns for results and for control
    #
    # input input_filepath: str
    # return tuple(data: pd.DataFrame, filename: str)

    # Read input
    data = pd.read_csv(input_filepath, index_col=0)
    filename = input_filepath.split('/')[-1]
    intra_only = all(data["XL_type"] == "intra")

    # Add result and control columns
    if intra_only:
        data["all_results"] = ''
        data[["pdb_id", "chain", "path", "pdb_method",
              "pdb_resolution", "best_res_pdb_method", "best_res_pdb_resolution"]] = '-'
    else:
        data[["all_results", "all_results_b"]] = ''
        data[["pdb_id", "pdb_id_b", "chain", "chain_b", "path", "path_b", "pdb_method", "pdb_method_b",
              "pdb_resolution", "pdb_resolution_b", "best_res_pdb_method", "best_res_pdb_method_b",
              "best_res_pdb_resolution", "best_res_pdb_resolution_b"]] = '-'

    return (data, filename), intra_only
