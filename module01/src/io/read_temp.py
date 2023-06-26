import os
import pandas as pd


def read_temp_search_save(data, tmp_filepath):
    # read temporary save file containing earlier uniprot search results, for quick reruns
    #
    # input data: pd.DataFrame, tmp_filepath: str
    # return read_temp_search_save: pd.DataFrame

    return pd.concat([data, pd.read_csv(tmp_filepath)], axis=1)
