import pandas as pd


def read_unipsearch_out(input_filename):
    # read output of uniprot_search (intra XLs only)
    #
    # input input_filename: str
    # return data: pd.DataFrame

    data = pd.read_csv(input_filename, index_col=0)
    return data[data.XL_type == "intra"]
