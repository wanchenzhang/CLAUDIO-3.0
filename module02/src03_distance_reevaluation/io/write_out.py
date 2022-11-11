
def write_output(data, filename, output_directory):
    # write dataset to output file as csv
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # no return

    out_filepath = f"{output_directory}{filename}.csv"
    data.to_csv(out_filepath, index=True)
