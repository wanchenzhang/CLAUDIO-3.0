import pandas as pd


def read_inputs(file1, file2):
    # read results of both reevaluations, return combined dataset
    #
    # input file1: str, file2: str
    # return data: pd.DataFrame, intra_only: bool

    data1 = pd.read_csv(file1, index_col=0)
    data2 = pd.read_csv(file2, index_col=0)
    data = data1.merge(data2)
    intra_only = all([x not in data.columns for x in ["unip_id_a", "unip_id_b", "seq_a", "seq_b"]])

    return data, intra_only
