
def write_output(data, filename, output_directory):
    # overwrite dataset to input filepath as csv
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # no return

    data.to_csv(f"{output_directory}{filename}_structdi.csv", index=True)
