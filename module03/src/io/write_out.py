
def write_output(data, filename, output_directory):
    # write dataset to output file as csv
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # no return

    out_filepath = f"{output_directory}{filename}_ops.csv"
    data[["pep_a", "pep_b", "pos_a", "pos_b", "unip_id_a", "unip_id_b", "pub", "XL_type", "homo_adjacency",
          "homo_int_overl", "homo_pep_overl"]].to_csv(out_filepath, index=True)
