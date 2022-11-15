import os
import pandas as pd


def read_temp_search_save(data, filepath):
    # read temporary save file containing earlier uniprot search results, for quick reruns
    #
    # input data: pd.DataFrame, filepath: str
    # return read_temp_search_save: pd.DataFrame

    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-4])
    tmp_filepath = f"{project_path}/data/temp/uniprot_search/" \
                   f"{'.'.join(filepath.split('.')[:-1])}_srtmp.{filepath.split('.')[-1]}"
    tmp_data = pd.read_csv(tmp_filepath)[["seq_a", "seq_b"]]
    return pd.concat([data, tmp_data], axis=1)
