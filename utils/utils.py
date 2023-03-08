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
