import pandas as pd


def read_inputs(file1, file2):
    # read results of both reevaluations, return combined dataset
    #
    # input file1: str, file2: str
    # return data: pd.DataFrame

    data1 = pd.read_csv(file1, index_col=0)
    data2 = pd.read_csv(file2, index_col=0)
    data = data1.merge(data2)

    return data
