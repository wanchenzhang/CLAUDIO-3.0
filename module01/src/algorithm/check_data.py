import sys

import pandas as pd
from Bio.Align import PairwiseAligner


def double_check_data(data, filename, output_directory):
    # Check given datapoints and correct it, if possible/needed
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # return data: pd.DataFrame

    log_text = ''
    ind = 0
    data_len = len(data.index)

    # Save known interactions (drop duplicates)
    known_inters = []

    # Save newly created datapoints (in order to ensure that they aren't duplicated multiple times)
    new_datapoints = []

    for i, row in data.iterrows():
        print(f"\r\t[{round(ind * 100 / data_len, 2)}%]", end='')

        # SUCCESS: Remove empty entry
        if type(row["seq_a"]) != str:
            log_text += f"{i}: empty entry. will be removed\n"
            data.drop(i, inplace=True)
            log_text += "\tSUCCESS\n"
            continue

        # FAIL: Neither K nor M is in the peptide (wrong peptide, or faulty interaction found)
        if 'K' not in row["pep_a"] and 'M' not in row["pep_a"]:
            log_text += f"{i}_a: pep_a does not contain K or M!\n"
            log_text += "\tFAIL\n"
        if 'K' not in row["pep_b"] and 'M' not in row["pep_b"]:
            log_text += f"{i}_b: pep_b does not contain K or M!\n"
            log_text += "\tFAIL\n"

        # REPLACE: Replace methionineoxide with methionine in peptide
        if "Mo" in row["pep_a"]:
            row["pep_a"] = row["pep_a"].replace("Mo", "M")
            log_text += f"{i}_a: Found 'Mo' in pep_a.\n"
            log_text += "\tREPLACE: Mo -> M\n"
        if "Mo" in row["pep_b"]:
            row["pep_b"] = row["pep_b"].replace("Mo", "M")
            log_text += f"{i}_b: Found 'Mo' in pep_b.\n"
            log_text += "\tREPLACE: Mo -> M\n"

        # ISSUE: Check if peptide is in full sequence as is (possible insertions are not considered)
        if row["pep_a"] not in row["seq_a"]:
            log_text += f"{i}_a: pep_a is not completely in full seq_a\n"
            log_text += f"\tISSUE\n"
        if row["pep_b"] not in row["seq_b"]:
            log_text += f"{i}_b: pep_b is not completely in full seq_b\n"
            log_text += f"\tISSUE\n"

        # Recompute lysin positions for site a and b
        row, log_text = check_pep_pos(i, row, 'a', log_text)
        row, log_text = check_pep_pos(i, row, 'b', log_text)
        # Update row
        data.loc[i, :] = row

        # ISSUE: Check if peptide occurs multiple times in full sequence as is (possible insertions are not considered)
        pep_found_multiple_times = False
        if row["seq_a"].count(row["pep_a"]) > 1:
            log_text += f"{i}_a: pep_a was found more than once in full seq_a\n"
            log_text += f"\tISSUE\n"
            pep_found_multiple_times = True
        if row["seq_b"].count(row["pep_b"]) > 1:
            log_text += f"{i}_b: pep_b was found more than once in full seq_b\n"
            log_text += f"\tISSUE\n"
            pep_found_multiple_times = True

        # If any peptide was found multiple times, create possible variations (add them later (*))
        if pep_found_multiple_times:
            datapoints, log_text = create_duplicates(row, log_text)
            new_datapoints.append(datapoints)

        # If not in known interactions, save unip ids and positions in known interactions, else drop datapoint
        interaction_info1 = [row["pos_a"], row["unip_id_a"], row["pos_b"], row["unip_id_b"]]
        interaction_info2 = [row["pos_b"], row["unip_id_b"], row["pos_a"], row["unip_id_a"]]
        if not(interaction_info1 in known_inters) and not(interaction_info2 in known_inters):
            known_inters.append(interaction_info1)
            known_inters.append(interaction_info2)
        else:
            for known_inter in known_inters:
                if (interaction_info1 == known_inter) or (interaction_info2 == known_inter):
                    known_i = known_inter[0]
                    log_text += f"{i}: entry is a DUPLICATE (see entry {known_i}). will be removed\n"
                    data.drop(i, inplace=True)
                    log_text += "\tSUCCESS\n"
                    break

        ind += 1
        print(f"\r\t[{round(ind * 100 / data_len, 2)}%]", end='')
    print()

    # Add aforementioned possible variations to the end of the dataset (*)
    for dp in new_datapoints:
        data = pd.concat([data, dp]).drop_duplicates().reset_index(drop=True)

    log_text += f"\n======================================================================\nFinal result:\n" \
                f"\tFAILS: {log_text.count('FAIL')}\n" \
                f"\tSUCCESSES: {log_text.count('SUCCESS')}\n" \
                f"\tISSUES: {log_text.count('ISSUE')}\n" \
                f"\tDUPLICATES: {log_text.count('DUPLICATE')}"

    # Save logfile
    output_path = f"{output_directory}{filename}.log"
    with open(output_path, 'w') as log_f:
        log_f.write(log_text)

    return data


def check_pep_pos(i, row, site, log_text):
    # Check given positions and correct it, if possible/needed
    #
    # input i: int, row: pd.Series, site: str, log_text: str
    # return row: pd.Series, log_text: str

    seq = row[f"seq_{site}"]
    pep_pos = int(row[f"pos_{site}"])
    peptide = row[f"pep_{site}"]

    wrong_pos = False
    try:
        # Check whether aminoacid at given position is either K or M
        if seq[pep_pos - 1] not in ['K', 'M']:
            log_text += f"{i}_{site}: residue in peptide at pos_{site} is not K or M.\n"
            wrong_pos = True
        # Check if residue is M, whether it is in position 1 of chain
        elif (seq[pep_pos - 1] == 'M') and (pep_pos != 1):
            log_text += f"{i}_{site}: residue in peptide at pos_{site} is M, but " \
                        f"should only be possible in pos=1.\n"
            wrong_pos = True
    except IndexError:
        # Catch given position out of bounds exception
        log_text += f"{i}_{site}: given pos_{site} is out of bounds of sequence.\n"
        wrong_pos = True

    # Recompute position if acid at position is not K or M, or if position out of bounds
    if wrong_pos:
        # Count peptide-as-is-occurences in full sequence
        pep_count_in_seq = seq.count(peptide)
        # Peptide occured exactly once in sequence (uniquely identified peptide)
        if pep_count_in_seq == 1:
            log_text += f"\tfound pep_{site} exactly once in full sequence\n"

            # Count Lysins
            k_count = peptide.count('K')
            if k_count == 0:
                log_text += f"\tno K found in pep_{site} (try M)\n"
                m_count = peptide.count('K')

                # FAIL: If no K and no M
                if m_count == 0:
                    log_text += "\tno M found either\n"
                    log_text += "\tFAIL\n"
                else:
                    # SUCCESS: M should be in first position, if true replace given position with 1 (Verify, FAIL if wrong)
                    if seq.find('M') != 0:
                        log_text += "\tno M found in expected first position of full sequence\n"
                        log_text += "\tFAIL\n"
                    else:
                        log_text += f"\tfound M at first position in pep_{site}\n"
                        log_text += f"\tREPLACE: pos_{site}: {pep_pos} -> 1\n"
                        row[f"pos_{site}"] = 1
                        log_text += f"\t\tVERIFY: new residue is M: {seq[0] == 'M'}\n"
                        if seq[0] == 'M':
                            log_text += f"\tSUCCESS\n"
                        else:
                            log_text += f"\tFAIL\n"

            # SUCCESS: Uniquely identified K, replace given position with its (Verify, FAIL if wrong)
            elif k_count == 1:
                log_text += f"\tfound K exactly once in pep_{site}\n"
                new_pos = seq.find(peptide) + peptide.find('K') + 1
                log_text += f"\tREPLACE: pos_{site}: {pep_pos} -> {new_pos}\n"
                row[f"pos_{site}"] = new_pos
                log_text += f"\t\tVERIFY: new residue is K: {seq[new_pos - 1] == 'K'}\n"
                if seq[new_pos - 1] == 'K':
                    log_text += f"\tSUCCESS\n"
                else:
                    log_text += f"\tFAIL\n"

            # More than one K in peptide (impossible to identify correct one, if not specified by data)
            else:
                log_text += f"\tfound K more than once in pep_{site}\n"
                # SUCCESS: Try k_pos specification (only for data by Liu et al.), if so replace given position with its (Verify, FAIL if wrong)
                try:
                    k_pos = int(row[f"k_pos_{site}"])
                    log_text += f"\tk_pos_{site} was given\n"
                    new_pos = seq.find(peptide) + k_pos
                    log_text += f"\tREPLACE: pos_{site}: {pep_pos} -> {new_pos}\n"
                    row[f"pos_{site}"] = new_pos
                    log_text += f"\t\tVERIFY: new residue is K: {seq[new_pos - 1] == 'K'}\n"
                    if seq[new_pos - 1] == 'K':
                        log_text += f"\tSUCCESS\n"
                    else:
                        log_text += f"\tFAIL\n"
                # FAIL: k_pos not given
                except:
                    log_text += f"\tk_pos_{site} was not given\n"
                    log_text += "\tFAIL\n"

        # Align peptide to sequence, if peptide found more than once or not found
        else:
            try:
                log_text += f"\tk_pos_{site} was given\n"
                log_text += "\talign peptide to sequence\n"
                k_pos = int(row[f"k_pos_{site}"])
                new_pos = realign_pep_to_seq(seq, peptide, k_pos)
                log_text += f"\tREPLACE: pos_{site}: {pep_pos} -> {new_pos}\n"
                row[f"pos_{site}"] = new_pos
                log_text += f"\t\tVERIFY: new residue is K: {seq[new_pos - 1] == 'K'}\n"
                if seq[new_pos - 1] == 'K':
                    log_text += f"\tSUCCESS\n"
                else:
                    log_text += f"\tFAIL\n"
            except:
                log_text += f"\tfound pep_{site} more than once or not even once in sequence\n"
                log_text += "\talignment attempt failed\n"
                log_text += "\tFAIL\n"

    return row, log_text


def realign_pep_to_seq(seq, peptide, k_pos):
    # Align peptide to sequence and based on k_pos in peptide compute k_pos in sequence, return None if this fails
    #
    # input seq: str, peptide: str, k_pos: int
    # return new_k_pos: int

    # Set alignment parameters
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -5
    aligner.extend_gap_score = 0
    aligner.match_score = 1
    aligner.mismatch_score = 0
    # aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    # Try aligning, if it fails return None
    try:
        alignment = aligner.align(seq, peptide)[0]
    except ValueError:
        return None
    print(f"\tAlignment:\n{alignment}")
    aligned_seq = str(alignment).split('\n')[0]
    aligned_pep = str(alignment).split('\n')[2]

    # Find k_pos in aligned peptide
    aligned_pep_i = -1
    aligned_pep_pos = 0
    for i, c in enumerate(aligned_pep):
        aligned_pep_i = i
        if c not in [' ', '-']:
            aligned_pep_pos += 1
        if aligned_pep_pos == k_pos:
            break

    # Find k_pos in aligned sequence and compute according new_k_pos in unaligned sequence, return new_k_pos if
    # successful, else return None
    try:
        _ = aligned_seq[aligned_pep_i + 1]
        new_k_pos = len(aligned_seq[:aligned_pep_i + 1].replace(' ', '').replace('-', ''))
        if new_k_pos > 0:
            al_seq = str(alignment).split('\n')[1]
            print(f"\tAligned site (index: {aligned_pep_i}):\n"
                  f"\t\t{aligned_seq[aligned_pep_i]}\n"
                  f"\t\t{al_seq[aligned_pep_i]}\n"
                  f"\t\t{aligned_pep[aligned_pep_i]}")
            return new_k_pos
        else:
            return None
    except IndexError:
        return None


def create_duplicates(row, log_text):
    # Check if peptides in datapoint occur multiple times in sequence, if so create dataset with possible permutations
    #
    # input row: pd.Series, log_text: str
    # return new_datapoints: pd.DataFrame, log_text: str

    new_datapoints = []

    shift_a, \
        shift_b = ('' for _ in range(2))

    # Make sure that the crosslinked residue can be uniquely identified for pep_a
    if ((row.pep_a.count('K') == 1) or row.k_pos_a.replace(' ', '')) and row.seq_a.count(row.pep_a) != 1:
        # Find all positions of pep_a
        pep_a_positions = [i for i in range(len(row.seq_a)) if row.seq_a[i:].startswith(row.pep_a)]
        # Compute shift for pos_a
        shift_a = row.pep_a.find('K') if row.pep_a.count('K') == 1 else int(row.k_pos_a) - 1
        # Add old pos_a if not already in position set
        old_pos_a = row.pos_a - shift_a - 1
        if old_pos_a not in pep_a_positions:
            pep_a_positions.append(old_pos_a)
        pep_a_positions = sorted(pep_a_positions)

    # Make sure that the crosslinked residue can be uniquely identified for pep_b
    if ((row.pep_b.count('K') == 1) or row.k_pos_b) and row.seq_b.count(row.pep_b) != 1:
        # Find all positions of pep_b
        pep_b_positions = [i for i in range(len(row.seq_b)) if row.seq_b[i:].startswith(row.pep_b)]
        # Compute shift for pos_b
        shift_b = row.pep_b.find('K') if row.pep_b.count('K') == 1 else int(row.k_pos_b) - 1
        # Add old pos_b if not already in position set
        old_pos_b = row.pos_b - shift_b - 1
        if old_pos_b not in pep_b_positions:
            pep_b_positions.append(old_pos_b)
        pep_b_positions = sorted(pep_b_positions)

    # Iterate through all possible peptide position permutations
    ind = 1
    if shift_a and shift_b:
        total_count = len(pep_a_positions) * len(pep_b_positions)
        for pos_a in pep_a_positions:
            for pos_b in pep_b_positions:
                # Create new datapoint and set up shifts
                new_datapoint = row.copy()
                new_pos_a = pos_a + shift_a + 1
                new_pos_b = pos_b + shift_b + 1

                # If the new datapoint is not equal to the existing one
                if (not ((row.pos_a == new_pos_a) and (row.pos_b == new_pos_b))) or \
                        ((row.pep_a == row.pep_b) and (row.pos_b == new_pos_a) and (row.pos_a == new_pos_b)):
                    log_text += f"\tcreate new datapoint ({ind} of {total_count}):\n"
                    ind += 1

                    # Verify new_pos_a
                    log_text += f"\t\tpos_a: {row.pos_a} -> {new_pos_a}\n" \
                                f"\t\tVERIFY: new residue is K: {row.seq_a[new_pos_a - 1] == 'K'}\n"
                    if row.seq_a[new_pos_a - 1] == 'K':
                        log_text += f"\t\t\tSUCCESS\n"
                    else:
                        log_text += f"\t\t\tFAIL\n"
                    new_datapoint.pos_a = new_pos_a

                    # Verify new_pos_b
                    log_text += f"\t\tpos_b: {row.pos_b} -> {new_pos_b}\n" \
                                f"\t\tVERIFY: new residue is K: {row.seq_b[new_pos_b - 1] == 'K'}\n"
                    if row.seq_b[new_pos_b - 1] == 'K':
                        log_text += f"\t\t\tSUCCESS\n"
                    else:
                        log_text += f"\t\t\tFAIL\n"
                    new_datapoint.pos_b = new_pos_b

                    new_datapoints.append(new_datapoint)
                else:
                    log_text += f"\tdatapoint already in dataset ({ind} of {total_count})\n"
                    ind += 1

    elif shift_a:
        total_count = len(pep_a_positions)
        for pos_a in pep_a_positions:
            # Create new datapoint and set up shifts
            new_datapoint = row.copy()
            new_pos_a = pos_a + shift_a + 1

            # If the new datapoint is not equal to the existing one
            if not (row.pos_a == new_pos_a):
                log_text += f"\tcreate new datapoint ({ind} of {total_count}):\n"
                ind += 1

                # Verify new_pos_a
                log_text += f"\t\tpos_a: {row.pos_a} -> {new_pos_a}\n" \
                            f"\t\tVERIFY: new residue is K: {row.seq_a[new_pos_a - 1] == 'K'}\n"
                if row.seq_a[new_pos_a - 1] == 'K':
                    log_text += f"\t\t\tSUCCESS\n"
                else:
                    log_text += f"\t\t\tFAIL\n"
                new_datapoint.pos_a = new_pos_a

                new_datapoints.append(new_datapoint)
            else:
                log_text += f"\tdatapoint already in dataset ({ind} of {total_count})\n"
                ind += 1

    elif shift_b:
        total_count = len(pep_b_positions)
        for pos_b in pep_b_positions:
            # Create new datapoint and set up shifts
            new_datapoint = row.copy()
            new_pos_b = pos_b + shift_b + 1

            # If the new datapoint is not equal to the existing one
            if not (row.pos_b == new_pos_b):
                log_text += f"\tcreate new datapoint ({ind} of {total_count}):\n"
                ind += 1

                # Verify new_pos_b
                log_text += f"\t\tpos_b: {row.pos_b} -> {new_pos_b}\n" \
                            f"\t\tVERIFY: new residue is K: {row.seq_b[new_pos_b - 1] == 'K'}\n"
                if row.seq_b[new_pos_b - 1] == 'K':
                    log_text += f"\t\t\tSUCCESS\n"
                else:
                    log_text += f"\t\t\tFAIL\n"
                new_datapoint.pos_b = new_pos_b

                new_datapoints.append(new_datapoint)
            else:
                log_text += f"\tdatapoint already in dataset ({ind} of {total_count})\n"
                ind += 1

    else:
        log_text += f"\tattempt to create new datapoints for peptide copies in full protein sequence failed\n" \
                    f"\tFAIL\n"

    return pd.concat(new_datapoints, axis=1).transpose(), log_text
