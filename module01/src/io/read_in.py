import pandas as pd
import sys


def read_inputfile(input_filepath, projections):
    # read inputfile and use projections to map column names to unified naming
    #
    # input input_filepath: str, projections: dict
    # return data: pd.DataFrame, intra_only: bool

    data = pd.read_csv(input_filepath)[list(projections.keys())]
    try:
        data = data.rename(columns=projections)
    except:
        print("Error! Given projection of column names failed! Check given dictionary!")
        sys.exit()

    # set primitive starting point for XL types
    data["XL_type"] = data.apply(lambda x: "intra" if x.unip_id_a == x.unip_id_b else "inter", axis=1)

    # Intra_only is True if all interactions are intra type
    intra_only = all([x == "intra" for x in data["XL_type"]])
    print(f"\t{'Data contained intra interactions only.' if intra_only else 'Data contained inter interactions.'}")

    return data, intra_only
