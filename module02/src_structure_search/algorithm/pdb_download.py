import os
import sys

import requests as r


def download_pdbs(dataset, search_tool, res_cutoff, output_directory):
    # Download pdb files either from RCSB or AlphaFold database (depending on earlier hhsearch or blastp search) into
    # output directory
    #
    # input dataset: pd.DataFrame, search_tool: str, res_cutoff: float, output_directory:str
    # return dataset: pd.DataFrame

    # clear output directory of old pdb file results
    clear_output_dir(search_tool, output_directory)

    ind = 0
    # Download pdb files for each datapoint
    for i, row in dataset.iterrows():
        # Iterate over results
        for j, res in enumerate((row['all_results'] + ' ').split(' ')):
            # If an entry was found in the rcsb database, download from there
            pdb_id = res.split('_')[0]
            chain = res.split('_')[1] if pdb_id else '-'
            a_pdb = ''
            if pdb_id:
                # Create custom pdb filename
                filename = f"{output_directory}{search_tool}_{pdb_id}.pdb"

                # If no similar pdb was already downloaded, then download
                if filename not in dataset["path"]:
                    # Download pdb as str text
                    a_URL = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    a_pdb = ''.join(r.post(a_URL).text)
                    # If ordinary download call fails attempt .cif call (for mmCIF file)
                    if a_pdb.startswith("<!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML 2.0//EN\">"):
                        a_URL = f"https://files.rcsb.org/download/{pdb_id}.cif"
                        a_pdb = ''.join(r.post(a_URL).text)

            # If no entry was found, attempt download from alphafold database instead
            else:
                # Create custom alphafold pdb filename
                pdb_id = f"af{row['unip_id']}"
                chain = 'A'
                filename = f"{output_directory}{search_tool}_{pdb_id}.pdb"

                # If no similar pdb was already downloaded, then download
                if filename not in dataset["path"]:
                    # Download pdb as str text
                    a_URL = f"https://alphafold.ebi.ac.uk/files/AF-{row['unip_id']}-F1-model_v1.pdb"
                    try:
                        a_pdb = r.get(a_URL).text
                    except ConnectionError as e:
                        print("No connection to AlphaFold API possible. Please try again later.")
                        print(e)
                        sys.exit()
                    # If no pdb retrieved from alphafold set filename to None, else save pdb with custom filename
                    if "NoSuchKey" in a_pdb:
                        filename = '-'
                        pdb_id = '-'
                        chain = '-'

            # Check whether method and resolution are accepted, return respective bool, method and resolution
            method_a_accepted, method_a, resolution_a = accept_resolution_method(a_pdb, pdb_id, res_cutoff)

            # Save method and resolution for best structure search result
            if j == 0:
                # Add method and resolution to dataset
                dataset.loc[i, "best_res_pdb_method"] = method_a
                dataset.loc[i, "best_res_pdb_resolution"] = resolution_a

            # Stop Iteration of results if result gets accepted
            if method_a_accepted:
                # Update pdb_id and chain in dataset
                dataset.loc[i, "pdb_id"] = pdb_id
                dataset.loc[i, "chain"] = chain
                # Add filename to paths
                dataset.loc[i, "path"] = filename
                # Add method and resolution to dataset
                dataset.loc[i, "pdb_method"] = method_a
                dataset.loc[i, "pdb_resolution"] = resolution_a

                # Save pdb text to new pdb file with custom name
                with open(filename, 'w') as f:
                    f.write(a_pdb)
                break

        ind += 1
        print(f"\r\t[{round(ind * 100 / len(dataset.index), 2)}%]", end='')
    print()

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
