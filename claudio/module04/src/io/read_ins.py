"""Reads in the outputs of Module02 and Module03 and merges them into a single DataFrame."""
import pandas as pd


def read_inputs(file1: str, file2: str):
    """
    read results of both reevaluations, return combined dataset
    
    Parameters
    ----------
    file1 : str,
    file2 : str

    Returns
    -------
    data : pd.DataFrame
    """

    data1 = pd.read_csv(
        file1, index_col=0, low_memory=False,
        dtype={"eucl_dist_raw": str, "topo_dist_raw": str},
    )
    data2 = pd.read_csv(file2, index_col=0)
    data = merge_datasets(data1, data2)
    different_chains = data["chain_a"] != data["chain_b"]
    overlap_evaluated = data["homo_pep_overl"].notna()

    data.loc[
        different_chains & overlap_evaluated,
        "homo_pep_overl",
    ] = False

    return data


def merge_datasets(df1: pd.DataFrame, df2: pd.DataFrame):
    """
    merge ops analysis dataset into structural distance analysis dataset. 
    Furthermore, ensure that results of ops analysis are integrated for new
    datapoints of structural distance analysis

    Parameters
    ----------
    df1 : pd.DataFrame,
    df2 : pd.DataFrame

    Returns
    -------
    df1 : pd.DataFrame
    """
    if "af_match_status" in df1.columns:
        return _merge_af_candidates(df1, df2)
    for column in df2.columns:
        if column not in df1.columns:
            if isinstance(df1.index[0],str):
                df1[column] = False
                for i in df1.index:
                    if '_' in str(i):
                        df1.loc[i, column] = df2.loc[int(str(i).split('_',maxsplit=1)[0]),
                                                     column]
                    else:
                        df1.loc[i, column] = df2.loc[int(i), column]
            else:
                df1[column] = df2[column].tolist()
    return df1

def _merge_af_candidates(structural, ops):
    """Attach M3 results using the original M1 row behind each AF candidate."""
    if not ops.index.is_unique:
        raise ValueError("M3 input indices must be unique")

    source_positions = {
        str(index): position
        for position, index in enumerate(ops.index)
    }
    positions = []
    for index in structural.index:
        source_index = str(index).split("_", maxsplit=1)[0]
        if source_index not in source_positions:
            raise ValueError(
                f"No M3 row found for structural candidate {index}"
            )
        positions.append(source_positions[source_index])

    aligned = ops.iloc[positions].copy()
    aligned.index = structural.index

    identity_columns = [
        "source_record_id", "candidate_id", "unip_id_a", "unip_id_b",
        "pos_a", "pos_b", "pep_a", "pep_b", "res_pos_a", "res_pos_b",
        "seq_a", "seq_b",
    ]
    for column in identity_columns:
        if column in structural.columns and column in aligned.columns:
            matches = structural[column].eq(aligned[column]).fillna(False)
            matches |= structural[column].isna() & aligned[column].isna()
            if not matches.all():
                raise ValueError(
                    f"M3 input does not match structural {column}"
                )

    result = structural.copy()
    ops_columns = {
        "homo_adjacency", "homo_int_overl", "homo_pep_overl",
    }
    for column in aligned.columns:
        if column not in result.columns or column in ops_columns:
            result[column] = aligned[column].to_numpy()

    return result