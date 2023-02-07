import pandas as pd


def read_in(input_filepath):
    # read outfile from uniprot search and isolate intra interactions only
    #
    # input input_filepath: str
    # return data: pd.DataFrame, intra_only: bool

    data = pd.read_csv(input_filepath, index_col=0)
    intra_only = all([x not in data.columns for x in ["unip_id_a", "unip_id_b", "seq_a", "seq_b"]])

    return data, intra_only
