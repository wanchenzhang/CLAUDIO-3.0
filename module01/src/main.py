import click
import sys
import os
import time
import ast
import pandas as pd

from module01.src.dict.default_projections import *
from module01.src.io.read_in import read_inputfile
from module01.src.io.read_temp import read_temp_search_save
from module01.src.algorithm.uniprot_search import do_uniprot_search
from module01.src.algorithm.check_data import double_check_data
from module01.src.algorithm.create_unique_list import create_list_of_unique_proteins
from module01.src.io.write_out import write_outputs


@click.command()
@click.option("-i", "--input-filepath", default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv")
@click.option("-p", "--projections", default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant))
@click.option("-s", "--uniprot-search", default=True)
@click.option("-t", "--search-tool", default="blastp")
@click.option("-o", "--output-directory", default="data/out/unique_protein_list")
@click.option("-bl", "--blast-bin", default=None)
@click.option("-bldb", "--blast-db", default="$BLASTDB")
@click.option("-hh", "--hhsearch-bin", default=None)
@click.option("-hhdb", "--hhsearch-db", default="$HHDB")
@click.option("-hhout", "--hhsearch-out", default="$HHOUT")
def main(input_filepath, projections, uniprot_search, search_tool, output_directory, blast_bin, blast_db, hhsearch_bin,
         hhsearch_db, hhsearch_out):
    print("Start Unique Protein List Tool")
    start_time = time.time()

    filename = '.'.join(input_filepath.split('/')[-1].split('.')[:-1])
    output_directory = output_directory if output_directory else '/'.join(input_filepath.split('/')[:-1])

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
    if not hhsearch_out.endswith('/'):
        hhsearch_out += '/'

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, projections, uniprot_search, search_tool, output_directory, blast_bin, blast_db,
                    hhsearch_bin, hhsearch_db, hhsearch_out):
        # Use projections to apply unified column names to input dataset
        # (for example see module01/src/dict/default_projections.py)
        projections = ast.literal_eval(projections)

        # Read input file: extract dataset, whether all datapoints are intra type interactions
        print("Read input")
        data, intra_only = read_inputfile(input_filepath, projections)

        # uniprot_search parameter is True actually perform a new search, else try to retrieve previous results
        # from temporary save file
        print("Retrieve UniProt sequences" if uniprot_search else "Retrieve UniProt sequences from temporary save")
        data = do_uniprot_search(data, intra_only, filename) if uniprot_search \
            else read_temp_search_save(data, filename)

        print("Check datapoints for inconsistencies")
        # Check datapoints for inconsistencies and correct them if possible (creates logfile in the process)
        data = double_check_data(data, filename, output_directory)
        print("Changes made to dataset written to log-file")

        # Write list of unique protein pairs and unique proteins overall
        print("Write unique protein and protein pairs lists")
        unique_pair_list, unique_proteins_list = create_list_of_unique_proteins(data, search_tool, blast_bin, blast_db,
                                                                                hhsearch_bin, hhsearch_db, hhsearch_out)

        # Write ouput csv
        print("Write output")
        write_outputs(data, unique_pair_list, unique_proteins_list, intra_only, filename, output_directory)

    print(f"\nEnd script (Elapsed time: {round(time.time() - start_time, 2)}s)")
    print("===================================")
    sys.exit()


def inputs_valid(input_filepath, projections, uniprot_search, search_tool, output_directory, blast_bin, blast_db,
                 hhsearch_bin, hhsearch_db, hhsearch_out):
    # check validity of inputted parameters
    #
    # input input_filepath: str, projections: str, uniprot_search: bool, search_tool: str, output_directory: str,
    # blast_bin: str/None, blast_db: str, hhsearch_bin: str/None, hhsearch_db: str, hhsearch_out: str
    # return inputs_valid: bool

    filename = input_filepath.split('/')[-1]
    # check whether an inputfile is specified
    if input_filepath:
        try:
            # check whether given dict can be read
            ast.literal_eval(projections)
            # if uniprot_search False then check whether temporary save file exists, return False if it fails,
            # else continue
            if not uniprot_search:
                try:
                    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-3])
                    project_path = project_path + '/' if project_path else ""
                    pd.read_csv(f"{project_path}data/temp/uniprot_search/{'.'.join(filename.split('.')[:-1])}_srtmp."
                                f"{filename.split('.')[-1]}")
                except FileNotFoundError:
                    print(f"Error! No temporary save file was found. Run the program with \"-s True\" to perform an "
                          f"actual search first (given: {uniprot_search}).")
                    return False
            if search_tool in ["blastp", "hhsearch"]:
                return True
            else:
                print(f"Error! Given search tool is neither blastp or hhsearch (given: {search_tool}).")
        except ValueError:
            print(f"Error! Could not construct dictionary from value given for \"projections\" parameter (given: "
                  f"{projections}).")
    else:
        print(f"Error! The parameter \"input-filepath\" was not given (given: {input_filepath}).")
    return False
