import os


def write_outputs(data, filename, output_directory):
    # write outputs
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # no return

    # save list of already written results by uniprot ids
    already_written = []

    for i, row in data.iterrows():
        unip_id = row.unip_id

        # only write new result files if results for uniprot entry were not written yet
        if (unip_id not in already_written) and (row.XL_type == "inter"):
            # create subdirectories, if not already in output directory
            output_subdirectory = f"{output_directory}homomers"
            if not os.path.exists(output_subdirectory):
                os.mkdir(output_subdirectory)
            output_subdirectory = f"{output_subdirectory}/{unip_id}"
            if not os.path.exists(output_subdirectory):
                os.mkdir(output_subdirectory)

            # if valid oligo-states were found, create fasta with respective number of sequence copies
            # (ex.: 3mer -> 3 seqs)
            if row.swiss_model_homology:
                for oligo_state in row.swiss_model_homology.split('_'):
                    num_mer = int(oligo_state.replace("mer", ''))
                    with open(f"{output_subdirectory}/{unip_id}_{oligo_state}.fasta", 'w') as f:
                        content = f">{unip_id}\n{row.seq}\n"
                        for _ in range(num_mer):
                            f.write(content)
            # else assume unknown multimer, with two sequences for now, but leave oligo-state field in filename empty
            # (user may decide here)
            else:
                with open(f"{output_subdirectory}/{unip_id}_new_multimer.fasta", 'w') as f:
                    content = f">{unip_id}\n{row.seq}\n"
                    f.write(content + content)

            # write all interactions for current uniprot id into a csv
            data[data.unip_id == unip_id][['pos_a', 'pos_b', 'swiss_model_homology', 'XL_type']]\
                .to_csv(f"{output_subdirectory}/{unip_id}_{row.swiss_model_homology}.csv", index=False)

            # add uniprot id to list of already finished results
            already_written.append(unip_id)

    # write full output dataset containing all given and computed information
    data.to_csv(f"{output_directory}{filename}_final.csv")
