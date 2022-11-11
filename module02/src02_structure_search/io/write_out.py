
def write_output(data, filepath):
    # write dataset to output file as csv
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # no return

    data.to_csv(filepath, index=True)
