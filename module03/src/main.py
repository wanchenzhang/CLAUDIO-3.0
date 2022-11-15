import click
import os
import sys
import time
import ast
import pandas as pd

from module02.src01_uniprot_search.dict.default_projections import *
from module02.src01_uniprot_search.main import main as run_unip_search
from module03.src.io.read_in import read_in
from module03.src.algorithm.signal_analysis import analyse_homo_signals
from module03.src.algorithm.create_histograms import create_homo_signal_histograms
from module03.src.io.write_out import write_output


@click.command()
@click.option("-i", "--input-filepath", default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv")
@click.option("-p", "--projections", default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant))
@click.option("-s", "--uniprot-search", default=True)
@click.option("-o", "--output-directory", default="data/out/homo_signal/")
def main(input_filepath, projections, uniprot_search, output_directory):
    print("Start Homo-signal analysis")
    start_time = time.time()

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, projections, uniprot_search, output_directory):
        output_directory = output_directory if output_directory else '/'.join(input_filepath.split('/')[:-1])
        if not output_directory.endswith('/'):
            output_directory += '/'

        # If no uniprot search run was performed already, do so now
        sqcs_path = f"{output_directory}{input_filepath.split('/')[-1]}.sqcs"
        if not os.path.exists(sqcs_path):
            try:
                run_unip_search(["-i", input_filepath, "-p", projections, "-s", uniprot_search, "-o", output_directory])
            except SystemExit:
                pass

        # Read dataset and add columns for results
        print("Read input")
        data = read_in(sqcs_path)

        # Analyse homo signals
        print("Analyse homo signals")
        data = analyse_homo_signals(data)

        # Create Homo-signal statistic histograms
        print("Create homo-signal histograms")
        create_homo_signal_histograms(data, input_filepath.split('/')[-1], output_directory)

        # Write Output
        print("Write output")
        write_output(data, input_filepath.split('/')[-1], output_directory)

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
