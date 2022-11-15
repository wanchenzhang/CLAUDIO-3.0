import click
import ast
import time
import sys

from module01.src.dict.default_projections import *
from module01.src.io.read_in import read_inputfile
from module01.src.algorithm.create_unique_list import create_list_of_unique_proteins
from module01.src.io.write_out import write_output


@click.command()
@click.option("-i", "--input-filepath", default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv")
@click.option("-p", "--projections", default=str(liu18_schweppe17_linked_residues_intra_homo_2370_nonredundant_unique))
@click.option("-t", "--search-tool", default="blastp")
@click.option("-o", "--output-directory", default="data/out/unique_protein_list")
def main(input_filepath, projections, search_tool, output_directory):
    print("Start Unique Protein List Tool")
    start_time = time.time()

    filename = ''.join(input_filepath.split('/')[-1].split('.')[:-1])
    output_directory = output_directory if output_directory else '/'.join(input_filepath.split('/')[:-1])

    if not output_directory.endswith('/'):
        output_directory += '/'

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, projections, search_tool, output_directory):
        # Use projections to apply unified column names to input dataset
        # (for example see module01/src/dict/default_projections.py)
        projections = ast.literal_eval(projections)

        # Read inputfile
        print("Read input")
        data = read_inputfile(input_filepath, projections)

        # Write list of unique protein pairs and unique proteins overall
        print("Write unique protein and protein pairs lists")
        unique_pair_list, unique_proteins_list = create_list_of_unique_proteins(data, search_tool)

        # Write ouput csv
        print("Write output")
        write_output(unique_pair_list, unique_proteins_list, filename, output_directory)

    print(f"\nEnd script (Elapsed time: {round(time.time() - start_time, 2)}s)")
    print("===================================")
    sys.exit()


def inputs_valid(input_filepath, projections, search_tool, output_directory):
    # check validity of inputted parameters
    #
    # input input_filepath: str, projections: str, search_tool: str, output_directory: str
    # return inputs_valid: bool

    filename = input_filepath.split('/')[-1]
    # check whether an inputfile is specified
    if input_filepath:
        try:
            # check whether given dict can be read
            ast.literal_eval(projections)
            if search_tool in ["blastp", "hhsearch"]:
                return True
            else:
                print(f"Error! Given search tool is neither blastp or hhsearch (given: {search_tool}).")
        except ValueError:
            print(f"Error! Could not construct dictionary from value given for \"projections\" parameter (given: {projections}).")
    else:
        print(f"Error! The parameter \"input-filepath\" was not given (given: {input_filepath}).")
    return False

