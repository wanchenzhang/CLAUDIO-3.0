import os
import socket
import sys
import time
import requests as r


from utils.utils import *


def download_pdbs(dataset, search_tool, intra_only, res_cutoff, output_directory, verbose_level):
    # Download pdb files either from RCSB or AlphaFold database (depending on earlier hhsearch or blastp search) into
    # output directory
    #
    # input dataset: pd.DataFrame, search_tool: str, intra_only: bool, res_cutoff: float, output_directory: str,
    # verbose_level: int
    # return dataset: pd.DataFrame

    # clear output directory of old pdb file results
    clear_output_dir(search_tool, output_directory)

    ind = 0
    # Download pdb files for each datapoint
    for i, row in dataset.iterrows():
        # Iterate over results
        for j, res in enumerate((row["all_results"] + ' ').split(' ')):
            # If an entry was found in the rcsb database, download from there
            pdb_id = res.split('_')[0]
            chain = res.split('_')[1] if pdb_id else '-'
            if not intra_only:
                chain_b = res.split('_')[2] if pdb_id else '-'
            pdb_file = ''
            if pdb_id:
                # Create custom pdb filename
                filename = f"{output_directory}{search_tool}_{pdb_id}.pdb"

                # If no similar pdb was already downloaded, then download
                if filename not in dataset["path"]:
                    # Download pdb as str text
                    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    pdb_file = download_pdb_from_db(url, 0, 5)
                    if pdb_file is None:
                        pdb_file = ''
                        filename = '-'
                        pdb_id = '-'
                        chain = '-'
                        if not intra_only:
                            chain_b = '-'
            # If no entry was found, attempt download from alphafold database instead, if it is an intra crosslink
            else:
                if intra_only or row["unip_id_a"] == row["unip_id_b"]:
                    # Create custom alphafold pdb filename
                    pdb_id = f"af{row['unip_id' if intra_only else 'unip_id_a']}"
                    chain = 'A'
                    if not intra_only:
                        chain_b = 'A'
                    filename = f"{output_directory}{search_tool}_{pdb_id}.pdb"

                    # If no similar pdb was already downloaded, then download
                    if filename not in dataset["path"]:
                        # Download pdb as str text
                        url = f"https://alphafold.ebi.ac.uk/files/AF-"\
                              f"{row['unip_id' if intra_only else 'unip_id_a']}-F1-model_v1.pdb"
                        pdb_file = download_pdb_from_db(url, 0, 5)

            # Check whether method and resolution are accepted, return respective bool, method and resolution
            method_a_accepted, method_a, resolution_a = accept_resolution_method(pdb_file, pdb_id, res_cutoff)

            # Save method and resolution for best structure search result
            if j == 0:
                # Add method and resolution to dataset
                dataset.loc[i, "best_res_pdb_method"] = method_a
                dataset.loc[i, "best_res_pdb_resolution"] = resolution_a

            # Stop Iteration of results if result gets accepted
            if method_a_accepted:
                # Update pdb_id and chain in dataset
                dataset.loc[i, "pdb_id"] = pdb_id
                if intra_only:
                    dataset.loc[i, "chain"] = chain
                else:
                    dataset.loc[i, "chain_a"] = chain
                    dataset.loc[i, "chain_b"] = chain_b
                # Add filename to paths
                dataset.loc[i, "path"] = filename
                # Add method and resolution to dataset
                dataset.loc[i, "pdb_method"] = method_a
                dataset.loc[i, "pdb_resolution"] = resolution_a

                # Save pdb text to new pdb file with custom name
                if not os.path.exists(filename):
                    with open(filename, 'w') as f:
                        f.write(pdb_file)
                break

        ind += 1
        verbose_print(f"\r\t[{round_self(ind * 100 / len(dataset.index), 2)}%]", 1, verbose_level, end='')
    verbose_print("", 1, verbose_level)

    return dataset


def clear_output_dir(search_tool, output_directory):
    # Clear output directory of pdb files starting with the specified search tool's name, to ensure that no old results
    # or other files intervene in later parsing of the results
    #
    # input search_tool: str, output_directory: str
    # no return

    pdb_files = [x for x in os.listdir(output_directory) if x.endswith(".pdb") and x.startswith(search_tool)]
    for f in pdb_files:
        os.remove(f"{output_directory}{f}")


def download_pdb_from_db(url, i_try, max_try):
    # Attempt pdb download from online database, either as .pdb- or .cif-file. If this fails retry until max retry is
    # reached.
    #
    # input url: str, i_try: int, max_try: int
    # output pdb_file: str/None
    try:
        if url.startswith("https://files.rcsb.org/"):
            # Attempt regular .pdb call from RCSB database
            pdb_file = ''.join(r.post(url).text)
            # If ordinary download call fails attempt .cif call (for mmCIF file)
            if pdb_file.startswith("<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML 2.0//EN\">"):
                cif_url = f"{'.'.join(url.split('.')[:-1])}.cif"
                pdb_file = ''.join(r.post(cif_url, timeout=1).text)
            return pdb_file
        else:
            # Attempt .pdb call from AlphaFold database
            pdb_file = r.get(url).text
            return pdb_file
    # Retry on timeout if not reached max_try already, else return None
    except (r.exceptions.Timeout, TimeoutError):
        if i_try >= max_try:
            return None
        else:
            time.sleep(1)
            return download_pdb_from_db(url, i_try + 1, max_try)
    # Break execution if no connection to database possible
    except (ConnectionError, socket.gaierror, r.exceptions.ConnectionError) as e:
        if i_try == max_try:
            print(f"No connection to {'RCSB' if url.startswith('https://files.rcsb.org/') else 'AlphaFold'} API possible. "
                  f"Please try again later.")
            print(e)
        return download_pdb_from_db(url, i_try + 1, max_try)


def accept_resolution_method(pdb, pdb_id, res_cutoff):
    # decide whether a pdb should be accepted based on the used method and its resolution
    #
    # input pdb: str, pdb_id: str, res_cutoff: float
    # return accept_resolution_method: boolean, method: str, resolution: str/float

    method, resolution = ('-' for _ in range(2))

    # If given pdb-file is empty return False
    if not pdb:
        return False, method, resolution
    # If pdb_id has length equal or higher than 5, it is an alphafold entry, return True
    elif len(pdb_id) >= 5:
        return True, "ALPHAFOLD", "ALPHAFOLD"
    # In other cases check whether experimental method for structure determination was diffraction or microscopy method,
    # and whether the resolution is below or equal to the threshhold, if so return True
    else:
        all_pdb_methods = ["ELECTRON DIFFRACTION", "FIBER DIFFRACTION", "FLUORESCENCE TRANSFER", "NEUTRON DIFFRACTION",
                           "SOLUTION NMR", "NMR", "THEORETICAL MODEL", "X-RAY DIFFRACTION", "X-RAY CRYSTALLOGRAPHY",
                           "ELECTRON CRYSTALLOGRAPHY", "ELECTRON MICROSCOPY"]
        accepted_pdb_methods = ["ELECTRON DIFFRACTION", "FIBER DIFFRACTION", "NEUTRON DIFFRACTION", "X-RAY DIFFRACTION",
                                "X-RAY CRYSTALLOGRAPHY", "ELECTRON CRYSTALLOGRAPHY", "ELECTRON MICROSCOPY"]
        accept_method = False
        accept_resolution = False

        for line in pdb.split('\n'):
            # If line startswith EXPDTA, it contains the information of the experimental method used for structure
            # determination
            if line.startswith("EXPDTA"):
                method = ' '.join([w for w in line.replace('  ', ' ').split() if w][1:])
                accept_method = method in accepted_pdb_methods
            # If line contains ANGSTROMS. and RESOLUTION. it contains the float value of the resolution, accept it if it
            # is below or equal to the threshhold
            elif ("ANGSTROMS." in line) and ("RESOLUTION." in line):
                resolution = float([w for w in line.replace('  ', ' ').split() if w][-2])
                accept_resolution = resolution <= res_cutoff
                break

        return accept_method and accept_resolution, method, resolution
