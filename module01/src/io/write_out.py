
def write_output(unique_pair_list, unique_proteins_list, filename, output_directory):
    # write outputs to files
    #
    # input unique_pair_list: pd.DataFrame, unique_proteins_list: pd.DataFrame, filename: str, output_directory: str
    # no return

    # write list of unique interaction pairs
    print(unique_pair_list.columns)
    print(unique_pair_list.head)
    output_path = f"{output_directory}{filename}_unique_pairs.csv"
    unique_pair_list.to_csv(output_path, index=True)

    # write list of unique proteins
    print(unique_proteins_list.columns)
    print(unique_proteins_list.head)
    output_path = f"{output_directory}{filename}_unique_proteins.csv"
    unique_proteins_list.to_csv(output_path, index=True)
