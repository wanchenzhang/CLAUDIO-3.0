import click
import sys
import os
import time
import ast
import pandas as pd

from module02.src01_uniprot_search.dict.default_projections import *
from module02.src01_uniprot_search.io.read_in import read_inputfile
from module02.src01_uniprot_search.algorithm.uniprot_search import do_uniprot_search
from module02.src01_uniprot_search.io.read_temp import read_temp_search_save
from module02.src01_uniprot_search.algorithm.check_data import double_check_data
from module02.src01_uniprot_search.io.write_out import write_sqcs_output


@click.command()
@click.option("-i", "--input-filepath", default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv")
@click.option("-p", "--projections", default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant))
@click.option("-s", "--uniprot-search", default=True)
@click.option("-o", "--output-directory", default="data/out/uniprot_search")
def main(input_filepath, projections, uniprot_search, output_directory):
    print("Start uniprot search")
    start_time = time.time()

    output_directory = output_directory if output_directory else '/'.join(input_filepath.split('/')[:-1])

    if not output_directory.endswith('/'):
        output_directory += '/'

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, projections, uniprot_search, output_directory):
        filename = input_filepath.split('/')[-1]

        print("Read input")
        # create projections dictionary from string input
        projections = ast.literal_eval(projections)

        # read input file: extract dataset, whether all datapoints are intra type interactions, whether a uniprot search
        # was already performed (e.g. uniprot IDs given), whether the uniprot sequences are already given, and otherwise
        # a list of search queries
        data, intra_only = read_inputfile(input_filepath, projections)

        # uniprot_search parameter is True actually perform a new search, else try to retrieve previous results
        # from temporary save file
        data = do_uniprot_search(data, intra_only, filename) if uniprot_search \
            else read_temp_search_save(data, filename)

        print("Check datapoints for inconsistencies")
        # Check datapoints for inconsistencies and correct them if possible (creates logfile in the process)
        data = double_check_data(data, intra_only, filename, output_directory)
        print("Changes made to dataset written to log-file")

        # Write output
        print("Write output")
        write_sqcs_output(data, intra_only, filename, output_directory)

    print(f"\nEnd script (Elapsed time: {round(time.time() - start_time, 2)}s)")
    print("===================================")
    sys.exit()


def inputs_valid(input_filepath, projections, uniprot_search, output_directory):
    # check validity of inputted parameters
    #
    # input input_filepath: str, projections: str, uniprot_search: bool, output_directory: str
    # return inputs_valid: bool

    filename = input_filepath.split('/')[-1]
    # check whether an inputfile is specified
    if input_filepath:
        try:
            # check whether given dict can be read
            ast.literal_eval(projections)
            if not uniprot_search:
                # if uniprot_search False then check whether temporary save file exists
                try:
                    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-3])
                    pd.read_csv(f"{project_path}/data/temp/uniprot_search/{'.'.join(filename.split('.')[:-1])}_srtmp."
                                f"{filename.split('.')[-1]}")
                    return True
                except FileNotFoundError:
                    print(f"Error! No temporary save file was found. Run the program with \"-s True\" to perform an "
                          f"actual search first (given: {uniprot_search}).")
            else:
                return True
        except ValueError:
            print(f"Error! Could not construct dictionary from value given for \"projections\" parameter "
                  f"(given: {projections}).")
    else:
        print(f"Error! The parameter \"input-filepath\" was not given (given: {input_filepath}).")
    return False
