import pandas as pd
import sys


def read_inputfile(input_filepath, projections):
    # read inputfile and use projections to map column names to unified naming
    #
    # input input_filepath: str, projections: dict
    # return data: pd.DataFrame

    data = pd.read_csv(input_filepath)
    try:
        data.rename(columns=projections, inplace=True)
    except:
        print("Error! Given projection of column names failed! Check given comma-separated list!")
        sys.exit()

    return data
