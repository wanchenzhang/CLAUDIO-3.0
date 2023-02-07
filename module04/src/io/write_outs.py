import os


def write_outputs(data, intra_only, filename, output_directory):
    # write outputs
    #
    # input data: pd.DataFrame, intra_only: bool, filename: str, output_directory: str
    # no return

    # save list of already written results by uniprot ids
    already_written = []

    # create subdirectories, if not already in output directory
    output_homo_dir = f"{output_directory}homomers"
    output_hetero_dir = f"{output_directory}heteromers"
    if (data.XL_type == "intra").any() and not os.path.exists(output_homo_dir):
        os.mkdir(output_homo_dir)
    if ((data.XL_type == "inter") & (~intra_only)).any() and not os.path.exists(output_hetero_dir):
        os.mkdir(output_hetero_dir)

    for i, row in data.iterrows():
        unip_id = row.unip_id if intra_only else \
            (f"{row.unip_id_a}_{row.unip_id_b}" if row.unip_id_a != row.unip_id_b else row.unip_id_a)

        # only write new result files if results for uniprot entry were not written yet
        if (unip_id not in already_written) and (row.XL_type == "inter"):
            is_intra = intra_only or row.unip_id_a == row.unip_id_b

            output_subdirectory = f"{output_homo_dir}/{unip_id}" if is_intra else f"{output_hetero_dir}/{unip_id}"
            if not os.path.exists(output_subdirectory):
                os.mkdir(output_subdirectory)

            # If intra crosslink write fasta with subsequent number of sequences
            if is_intra:
                # if valid oligo-states were found, create fasta with respective number of sequence copies
                # (ex.: 3mer -> 3 seqs)
                if row.swiss_model_homology:
                    for oligo_state in row.swiss_model_homology.split('_'):
                        num_mer = int(oligo_state.replace("homo", '').replace("mer", ''))
                        with open(f"{output_subdirectory}/{unip_id}_{oligo_state}.fasta", 'w') as f:
                            content = f">{unip_id}\n{row.seq if intra_only else row.seq_a}\n"
                            for _ in range(num_mer):
                                f.write(content)
                # else assume unknown homomer, with two sequences for now, but leave oligo-state field in filename empty
                # (user may decide here)
                else:
                    with open(f"{output_subdirectory}/{unip_id}_new_homomer.fasta", 'w') as f:
                        content = f">{unip_id}\n{row.seq if intra_only else row.seq_a}\n"
                        f.write(content + content)
            # Else write multi fasta with involved sequences for heteromer complex
            else:
                # If heteromeric state given by SWISS-MODEL write title plainly, else write identifier for new heteromer
                multifasta_filename = f"{unip_id}_heteromer.fasta" \
                    if "hetero" in row.swiss_model_homology else f"{unip_id}_new_heteromer.fasta"
                with open(f"{output_subdirectory}/{multifasta_filename}", 'w') as f:
                    content = f">{row.unip_id_a}\n{row.seq_a}\n\n>{row.unip_id_b}\n{row.seq_b}"
                    f.write(content)

            # write all interactions for current uniprot id into a csv
            if intra_only:
                unip_arg = data.unip_id == unip_id
            elif row.unip_id_a == row.unip_id_b:
                unip_arg = data.unip_id_a == unip_id
            else:
                unip_arg = ((data.unip_id_a == unip_id.split('_')[0]) & (data.unip_id_b == unip_id.split('_')[1])) | \
                           ((data.unip_id_a == unip_id.split('_')[1]) & (data.unip_id_b == unip_id.split('_')[0]))
            data[unip_arg][['pos_a', 'pos_b', 'swiss_model_homology', 'XL_type']]\
                .to_csv(f"{output_subdirectory}/{unip_id}_{row.swiss_model_homology}.csv", index=False)

            # add uniprot id to list of already finished results
            if is_intra:
                already_written.append(unip_id)
            else:
                already_written.append(unip_id)
                already_written.append(f"{unip_id.split('_')[1]}_{unip_id.split('_')[0]}")

    # write full output dataset containing all given and computed information
    data.to_csv(f"{output_directory}{filename}_final.csv")
