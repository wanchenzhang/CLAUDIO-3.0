import click
import sys
import os
import time

from module01.src.io.read_in import read_inputfile
from module01.src.io.read_temp import read_temp_search_save
from module01.src.algorithm.uniprot_search import do_uniprot_search
from module01.src.algorithm.check_data import double_check_data
from module01.src.algorithm.create_unique_list import create_list_of_unique_proteins
from module01.src.io.write_out import write_outputs

from utils.utils import *


@click.command()
@click.option("-i", "--input-filepath", default="data/in/sample_data.csv")
@click.option("-it", "--input-temppath", default="")
@click.option("-p", "--projections", default="peptide1,peptide2,position1,position2,k_pos1,k_pos2,entry1,entry2")
@click.option("-s", "--uniprot-search", default=True)
@click.option("-x", "--xl-residues", default="K,M:N:1")
@click.option("-t", "--search-tool", default="blastp")
@click.option("-o", "--output-directory", default="data/out/unique_protein_list")
@click.option("-bl", "--blast-bin", default=None)
@click.option("-bldb", "--blast-db", default="$BLASTDB")
@click.option("-hh", "--hhsearch-bin", default=None)
@click.option("-hhdb", "--hhsearch-db", default="$HHDB")
@click.option("-v", "--verbose-level", default=2)
def main(input_filepath, input_temppath, projections, uniprot_search, xl_residues, search_tool, output_directory,
         blast_bin, blast_db, hhsearch_bin, hhsearch_db, verbose_level):
    verbose_print("Start Unique Protein List Tool", 0, verbose_level)
    start_time = time.time()

    input_filepath, input_temppath, output_directory, \
        blast_bin, blast_db, hhsearch_bin, hhsearch_db = \
        (path.replace('\\', '/') for path in
         [input_filepath, input_temppath, output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db])

    filename = '.'.join(input_filepath.split('/')[-1].split('.')[:-1])

    # Check output directory
    output_directory = create_out_path(output_directory, input_filepath)

    # Create temporary dirs
    uniprot_search_temp_dir = create_temp_dir(input_temppath, "uniprot_search")
    unique_protein_temp_dir = create_temp_dir(input_temppath, "unique_protein_list")

    # Convert directory paths to literals if None
    if blast_bin == "None":
        blast_bin = None
    if hhsearch_bin == "None":
        hhsearch_bin = None

    # Add '/' to end of directory paths if not there
    if not output_directory.endswith('/'):
        output_directory += '/'
    if (blast_bin is not None) and (not blast_bin.endswith('/')):
        blast_bin += '/'
    if not blast_db.endswith('/'):
        blast_db += '/'
    if (hhsearch_bin is not None) and (not hhsearch_bin.endswith('/')):
        hhsearch_bin += '/'
    if not hhsearch_db.endswith('/'):
        hhsearch_db += '/'

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, uniprot_search_temp_dir, unique_protein_temp_dir, projections, uniprot_search,
                    xl_residues, search_tool, output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                    verbose_level):
        # Use projections to apply unified column names to input dataset
        # (for example see module01/src/dict/default_projections.py)
        new_keys = ["pep_a", "pep_b", "pos_a", "pos_b", "res_pos_a", "res_pos_b", "unip_id_a", "unip_id_b"]
        projections = {projections.split(',')[i]: new_keys[i] for i in range(len(new_keys))}

        # Define dataset for crosslink residues including possible positions and atom types
        df_xl_res = build_xl_dataset(xl_residues)

        # Read input file
        verbose_print("Read input", 0, verbose_level)
        data = read_inputfile(input_filepath, projections)

        # uniprot_search parameter is True actually perform a new search, else try to retrieve previous results
        # from temporary save file
        tmp_filepath = f"{uniprot_search_temp_dir}{filename}_srtmp.csv"
        if (not uniprot_search) and os.path.exists(tmp_filepath):
            verbose_print("Retrieve UniProt sequences from temporary save", 0, verbose_level)
            read_temp_search_save(data, tmp_filepath)
        else:
            verbose_print("Retrieve UniProt sequences from UniProtKB", 0, verbose_level)
            do_uniprot_search(data, tmp_filepath, verbose_level)

        # Check datapoints for inconsistencies and correct them if possible (creates logfile in the process)
        verbose_print("Check datapoints for inconsistencies", 0, verbose_level)
        data = double_check_data(data, filename, df_xl_res, output_directory, verbose_level)
        verbose_print("Changes made to dataset written to log-file", 0, verbose_level)

        # Write list of unique protein pairs and unique proteins overall
        verbose_print("Create unique protein list", 0, verbose_level)
        unique_proteins_list = create_list_of_unique_proteins(data, unique_protein_temp_dir, search_tool, blast_bin,
                                                              blast_db, hhsearch_bin, hhsearch_db, verbose_level)

        # Clean dataset for output
        data = clean_dataset(data)

        # Write ouput csv
        verbose_print("Write output", 0, verbose_level)
        write_outputs(data, unique_proteins_list, filename, output_directory, verbose_level)

    verbose_print(f"\nEnd script (Elapsed time: {round_self(time.time() - start_time, 2)}s)", 0, verbose_level)
    verbose_print("===================================", 0, verbose_level)
    sys.exit()


def inputs_valid(input_filepath, uniprot_search_temp_dir, unique_protein_temp_dir, projections, uniprot_search,
                 xl_residues, search_tool, output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                 verbose_level):
    # check validity of inputted parameters
    #
    # input input_filepath: str, uniprot_search_temp_dir: str, unique_protein_temp_dir: str, projections: str,
    # uniprot_search: bool, xl_residues: str, search_tool: str, output_directory: str, blast_bin: str/None,
    # blast_db: str, hhsearch_bin: str/None, hhsearch_db: str, verbose_level: int
    # return inputs_valid: bool

    filename = input_filepath.split('/')[-1]
    # check whether an inputfile is specified
    if input_filepath:
        # check whether the number of comma-separated values is acceptable
        if len(projections.split(',')) == 8:
            # if uniprot_search False then check whether temporary save file exists, return False if it fails,
            # else continue
            if not uniprot_search:
                try:
                    pd.read_csv(f"{uniprot_search_temp_dir}{'.'.join(filename.split('.')[:-1])}_srtmp."
                                f"{filename.split('.')[-1]}")
                except FileNotFoundError:
                    print(f"Error! No temporary save file was found. Run the program with \"-s True\" to perform an "
                          f"actual search first (given: {uniprot_search}).")
                    return False
            # check whether xl_residues can be turned into a proper DataFrame, else return False
            try:
                if build_xl_dataset(xl_residues) is None:
                    return False

                # check whether specified structure search tool is either blastp or hhsearch
                if search_tool in ["blastp", "hhsearch"]:
                    return True
                else:
                    print(f"Error! Given search tool is neither blastp or hhsearch (given: {search_tool}).")
            except:
                print(f"Error! Could not properly parse xl_residues for accepted crosslinked residues "
                      f"(given: {xl_residues}).")
        else:
            print(f"Error! Could not find all necessary keys in \"projections\" parameter (given: "
                  f"{projections}).")
    else:
        print(f"Error! The parameter \"input-filepath\" was not given (given: {input_filepath}).")
    return False
