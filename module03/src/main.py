import click
import os
import sys
import time
import ast
import pandas as pd

from module03.src.io.read_in import read_in
from module03.src.algorithm.signal_analysis import analyse_homo_signals
from module03.src.algorithm.create_histograms import create_homo_signal_histograms
from module03.src.io.write_out import write_output


@click.command()
@click.option("-i", "--input-filepath", default="data/out/unique_protein_list/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.sqcs")
@click.option("-o", "--output-directory", default="data/out/homo_signal/")
def main(input_filepath, output_directory):
    print("Start Homo-signal analysis")
    start_time = time.time()

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, output_directory):
        output_directory = output_directory if output_directory else '/'.join(input_filepath.split('/')[:-1])
        if not output_directory.endswith('/'):
            output_directory += '/'

        # Read dataset and add columns for results
        print("Read input")
        data = read_in(input_filepath)

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


def inputs_valid(input_filepath, output_directory):
    # check validity of inputted parameters
    #
    # input input_filepath: str, output_directory: str
    # return inputs_valid: bool

    # check whether an inputfile with the extension .sqcs is specified
    if input_filepath and input_filepath.endswith(".sqcs"):
        return True
    else:
        print(f"Error! The parameter \"input-filepath\" was not given or was not ending with the '.sqcs' extension "
              f"(given: {input_filepath}).")
    return False
