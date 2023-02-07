import pandas as pd


def read_unipsearch_out(input_filename):
    # read output of uniprot_search (intra XLs only)
    #
    # input input_filename: str
    # return data: pd.DataFrame, intra_only: bool

    data = pd.read_csv(input_filename, index_col=0)
    intra_only = all([x not in data.columns for x in ["unip_id_a", "unip_id_b", "seq_a", "seq_b"]])

    return data, intra_only
