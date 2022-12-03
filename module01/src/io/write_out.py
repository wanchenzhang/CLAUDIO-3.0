
def write_outputs(data, unique_proteins_list, filename, output_directory):
    # write outputs to files: One csv-file for the users information on unique proteins and one
    # csv-file marked with the extension .sqcs, signaling CLAUDIO that this dataset is processed for its uses
    #
    # input data: pd.DataFrame, unique_proteins_list: pd.DataFrame, filename: str, output_directory: str
    # no return

    # write list of unique proteins
    print("\tWrite list of unique proteins")
    output_path = f"{output_directory}{filename}_unique_proteins.csv"
    unique_proteins_list.to_csv(output_path, index=True)

    print("\tWrite SQCS-file (a CSV-file) for further steps in CLAUDIO")
    output_path = f"{output_directory}{filename}.sqcs"
    data.to_csv(output_path, index=True)
