import pandas as pd


def read_inputs(file1, file2):
    # read results of both reevaluations, return combined dataset
    #
    # input file1: str, file2: str
    # return data: pd.DataFrame, intra_only: bool

    data1 = pd.read_csv(file1, index_col=0)
    data2 = pd.read_csv(file2, index_col=0)
    data = merge_datasets(data1, data2)
    intra_only = all([x not in data.columns for x in ["unip_id_a", "unip_id_b", "seq_a", "seq_b"]])

    return data, intra_only


def merge_datasets(df1, df2):
    # merge ops analysis dataset into structural distance analysis dataset. Furthermore,
    # ensure that results of ops analysis are integrated for new datapoints of structural distance analysis
    #
    # input df1: pd.DataFrame, df2: pd.DataFrame
    # return df1: pd.DataFrame
    for column in df2.columns:
        if column not in df1.columns:
            df2_col_content = df2[column].to_list()
            if type(df1.index[0]) == str:
                for i in df1.index:
                    if '_' in i:
                        df2_col_content.append(df2.loc[int(i.split('_')[0]), column])
            df1[column] = df2_col_content
    return df1
