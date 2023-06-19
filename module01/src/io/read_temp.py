import os
import pandas as pd


def read_temp_search_save(data, input_temppath, filename):
    # read temporary save file containing earlier uniprot search results, for quick reruns
    #
    # input data: pd.DataFrame, input_temppath: str, filename: str
    # return read_temp_search_save: pd.DataFrame

    if input_temppath == '/':
        project_path = '/'.join(os.path.abspath(__file__).split('/')[:-4])
        project_path = project_path + '/' if project_path else ""
        input_temppath = project_path
    tmp_filepath = f"{input_temppath}uniprot_search/{filename}_srtmp.csv"
    tmp_data = pd.read_csv(tmp_filepath)

    return pd.concat([data, tmp_data], axis=1)
