import pandas as pd


def verbose_print(print_string, threshold, verbose_level, end='\n'):
    # print given string, if verbose_level is higher than threshold
    #
    # input print_string: str, threshold: int, verbose_level: int, end: str
    # no return

    if verbose_level > threshold:
        print(print_string, end=end)


def build_xl_dataset(xl_residues):
    # build residue dataset from comma-separated xl_residues input string, specifying residue, atom type, and position
    #
    # input xl_residues: str
    # return df_xl_res: pd.DataFrame

    res_list, pos_list, atom_list = ([] for _ in range(3))

    for s in xl_residues.replace(';', ',').split(','):
        if ':' in s:
            if s.count(':') != 2:
                print(f"Error! Found ':' in one xl_res input, but less or more than two-times. If you wish "
                      f"to specify either the position or the atom type make sure you always add two ':' "
                      f"in the input (specific: {s}, full: {xl_residues}).")
                return None
            res_list.append(s.split(':')[0])
            atom_list.append(s.split(':')[1] if s.split(':')[1] else "CB")
            pos_list.append(int(s.split(':')[2]) if s.split(':')[2] else 0)
        else:
            res_list.append(s)
            atom_list.append("CB")
            pos_list.append(0)
    df_xl_res = pd.DataFrame()
    df_xl_res["res"] = res_list
    df_xl_res["atom"] = atom_list
    df_xl_res["pos"] = pos_list

    return df_xl_res


def round_self(value, decimals):
    # simple decimal rounding function (python by itself has a tendency to round fragmented with the built-in function)
    #
    # input value: float, decimals: int
    # return rounded_value: float/int

    # If decimal less than 1, the resulting value will be an integer
    if pd.isna(value):
        return float("Nan")
    if decimals < 1:
        rounded_value = int(f"{int((value * (10 ** decimals)) + .5) / (10 ** decimals):.{decimals}f}")
        return rounded_value
    # Else, the resulting value will be a float
    else:
        rounded_value = float(f"{int((value * (10 ** decimals)) + .5) / (10 ** decimals):.{decimals}f}")
        return rounded_value


def clean_dataset(data):
    # Clean structure data dataset for final outputs
    #
    # input data: pd.DataFrame
    # return data: pd.DataFrame

    # Drop datapoints with ident index, if data values are all identical
    if "pdb_id" in data.columns:
        drop_indeces,\
            already_checked = ([] for _ in range(2))
        for i, row in data.iterrows():
            if ('_' in str(i)) and (i not in already_checked):
                drop_criteria = (data.pos_a == row.pos_a) & (data.pos_b == row.pos_b) & \
                                (data.pep_a == row.pep_a) & (data.pep_b == row.pep_b) & \
                                ((data.res_pos_a == row.res_pos_a) | (pd.isna(data.res_pos_a) & pd.isna(row.res_pos_a))) & \
                                ((data.res_pos_b == row.res_pos_b) | (pd.isna(data.res_pos_b) & pd.isna(row.res_pos_b)))
                if all([rename_col in data.columns for rename_col in ["chain_a", "chain_b", "pdb_id"]]):
                    drop_criteria = drop_criteria & (data.chain_a == row.chain_a) & (data.chain_b == row.chain_b) & \
                                    (data.pdb_id == row.pdb_id)
                if all([rename_col in data.columns for rename_col in ["evidence"]]):
                    drop_criteria = drop_criteria & (data.evidence == row.evidence)

                if len(data[drop_criteria].index) > 1:
                    drop_indeces.extend(list(data[drop_criteria].index)[1:])
                already_checked.extend(list(data[drop_criteria].index))
        data = data.drop(index=drop_indeces)

    # Drop specified data columns
    for drop_col in ["all_results", "path", "best_res_pdb_method", "best_res_pdb_resolution", "eucl_dist",
                     "res_criteria_fulfilled", "res_crit_a", "res_crit_b", "method_a", "method_b"]:
        if drop_col in data.columns:
            data = data.drop(drop_col, axis=1)

    # Rename certain result columns
    if all([rename_col in data.columns for rename_col in ["eucl_dist_tplk", "topo_dist_tplk"]]):
        data = data.rename(columns={"eucl_dist_tplk": "eucl_dist", "topo_dist_tplk": "topo_dist"})

    # Ascertain data types
    data = data.astype({"pos_a": int, "pos_b": int, "pep_a": str, "pep_b": str, "res_pos_a": int, "res_pos_b": int},
                       errors="ignore")
    if "unip_id" in data.columns:
        data = data.astype({"unip_id": str, "seq": str}, errors="ignore")
    else:
        data = data.astype({"unip_id_a": str, "unip_id_b": str, "seq_a": str, "seq_b": str}, errors="ignore")
    if "pdb_id" in data.columns:
        data = data.astype({"pdb_id": str, "pdb_method": str, "pdb_resolution": str,
                            "pdb_pos_a": int, "pdb_pos_b": int, "pLDDT_a": float, "pLDDT_b": float,
                            "topo_dist": float}, errors="ignore")
        if "unip_id" in data.columns:
            data.astype({"chain": str}, errors="ignore")
        else:
            data.astype({"chain_a": str, "chain_b": str}, errors="ignore")
    if "homo_pep_overl" in data.columns:
        data = data.astype({"homo_adjacency": float, "homo_int_overl": float, "homo_pep_overl": bool}, errors="ignore")
    if "evidence" in data.columns:
        data = data.astype({"evidence": str, "XL_type": str, "swiss_model_homology": str}, errors="ignore")

    # Reduce multichain examples to fitting ones (if present)
    if "evidence" in data.columns:
        drop_indeces,\
            already_checked = ([] for _ in range(2))
        for i, row in data.iterrows():
            if ('_' in str(i)) and (i not in already_checked):
                multi_chain_criteria = (data.pos_a == row.pos_a) & (data.pos_b == row.pos_b) & \
                                       (data.pep_a == row.pep_a) & (data.pep_b == row.pep_b) & \
                                       ((data.res_pos_a == row.res_pos_a) |
                                        (pd.isna(data.res_pos_a) & pd.isna(row.res_pos_a))) & \
                                       ((data.res_pos_b == row.res_pos_b) |
                                        (pd.isna(data.res_pos_b) & pd.isna(row.res_pos_b)))
                multi_chain_set = data[multi_chain_criteria]
                if not multi_chain_set[data.evidence == ''].empty:
                    drop_indeces.extend(list(multi_chain_set[data.evidence != ''].index))
                already_checked.extend(list(multi_chain_set.index))
        data = data.drop(index=drop_indeces)

    # Sort rows
    if any([type(i) == str for i in data.index]):
        data["sort_index_i"] = [int(i.split('_')[0]) for i in data.index]
        data["sort_index_j"] = [int(i.split('_')[1]) if '_' in i else 1 for i in data.index]
        data = data.sort_values(["sort_index_i", "sort_index_j"]).drop("sort_index_i", axis=1).drop("sort_index_j", axis=1)

    return data
