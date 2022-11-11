import pandas as pd
import sys


def read_inputfile(input_filepath, projections):
    # read inputfile and use projections to map column names to unified naming
    #
    # input input_filepath: str, projections: dict
    # return data: pd.DataFrame

    data = pd.read_csv(input_filepath)[list(projections.keys())]
    try:
        data = data.rename(columns=projections)
    except:
        print("Error! Given projection of column names failed! Check given dictionary!")
        sys.exit()

    return data
