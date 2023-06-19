import os
import pandas as pd


def read_temp_search_save(data, uniprot_search_temp_dir, filename):
    # read temporary save file containing earlier uniprot search results, for quick reruns
    #
    # input data: pd.DataFrame, input_temppath: str, filename: str
    # return read_temp_search_save: pd.DataFrame

    tmp_filepath = f"{uniprot_search_temp_dir}{filename}_srtmp.csv"
    tmp_data = pd.read_csv(tmp_filepath)

    return pd.concat([data, tmp_data], axis=1)
