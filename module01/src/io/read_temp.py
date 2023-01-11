import os
import pandas as pd


def read_temp_search_save(data, filename):
    # read temporary save file containing earlier uniprot search results, for quick reruns
    #
    # input data: pd.DataFrame, filename: str
    # return read_temp_search_save: pd.DataFrame

    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-4])
    project_path = project_path + '/' if project_path else ""
    tmp_filepath = f"{project_path}data/temp/uniprot_search/" \
                   f"{filename}_srtmp.csv"
    tmp_data = pd.read_csv(tmp_filepath)

    return pd.concat([data, tmp_data], axis=1)
