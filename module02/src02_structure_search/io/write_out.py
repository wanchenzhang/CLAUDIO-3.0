
def write_output(data, filepath):
    # overwrite dataset to input filepath as csv
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # no return

    data.to_csv(filepath, index=True)
