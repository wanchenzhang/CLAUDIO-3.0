import click
import sys
import time
import os

from module02.src_structure_search.io.read_input import read_in
from module02.src_structure_search.algorithm.structure_search import structure_search
from module02.src_structure_search.io.read_temp import read_temp_file
from module02.src_structure_search.algorithm.pdb_download import download_pdbs
from module02.src_structure_search.algorithm.chain_copies import create_ident_chain_copies
from module02.src_structure_search.io.write_out import write_output

from utils.utils import *


@click.command()
@click.option("-i", "--input-filepath", default="data/out/unique_protein_list/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.sqcs")
@click.option("-s", "--do-structure-search", default=True)
@click.option("-t", "--search-tool", default="blastp")
@click.option("-e", "--e-value", default=1e-5)
@click.option("-q", "--query-id", default=90.0)
@click.option("-c", "--coverage", default=50.0)
@click.option("-r", "--res-cutoff", default=6.5)
@click.option("-o", "--output-directory", default="data/out/structure_search")
@click.option("-bl", "--blast-bin", default=None)
@click.option("-bldb", "--blast-db", default="$BLASTDB")
@click.option("-hh", "--hhsearch-bin", default=None)
@click.option("-hhdb", "--hhsearch-db", default="$HHDB")
@click.option("-v", "--verbose-level", default=3)
def main(input_filepath, do_structure_search, search_tool, e_value, query_id, coverage, res_cutoff, output_directory,
         blast_bin, blast_db, hhsearch_bin, hhsearch_db, verbose_level):
    verbose_print("Start structure search", 0, verbose_level)
    start_time = time.time()

    project_path = '/'.join(os.path.abspath(__file__).split('/')[:-3])
    project_path = project_path + '/' if project_path else ""
    temp_save_search_path = f"{project_path}data/temp/structure_search/"
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

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, do_structure_search, search_tool, e_value, query_id, coverage, res_cutoff,
                    output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db, verbose_level):
        e_value = float(e_value)
        query_id = float(query_id)
        coverage = float(coverage)
        res_cutoff = float(res_cutoff)

        # Read dataset and add columns for results
        verbose_print("Read input", 0, verbose_level)
        data, filename = read_in(input_filepath)

        # If given variable do_structure_search is True perform new search, else retrieve results from earlier temporary
        # save file
        if do_structure_search:
            verbose_print(f"Perform {search_tool} search", 0, verbose_level)
            data = structure_search(data, filename, search_tool, e_value, query_id, coverage, temp_save_search_path,
                                    blast_bin, blast_db, hhsearch_bin, hhsearch_db, verbose_level)
        else:
            verbose_print("Read from temporary save file", 0, verbose_level)
            data = read_temp_file(data, filename, search_tool, temp_save_search_path)

        # Download structure files from RCSB database into structures subdirectory, if search tool found a proper
        # result, else use uniprot ID in order to attempt retrieval of matching AlphaFold entry
        verbose_print("Download structures from RCSB database, or AlphaFold database if not found there", 0,
                      verbose_level)
        if not os.path.exists(f"{output_directory}structures"):
            os.mkdir(f"{output_directory}structures")
        data = download_pdbs(data, search_tool, res_cutoff, f"{output_directory}structures/", verbose_level)

        # Create copy datapoints with marked indeces for alternative but identical chains
        data = create_ident_chain_copies(data)

        # Write new output
        write_output(data, filename, output_directory)

    verbose_print(f"\nEnd script (Elapsed time: {round_self(time.time() - start_time, 2)}s)", 0, verbose_level)
    verbose_print("===================================", 0, verbose_level)
    sys.exit()


def inputs_valid(input_filepath, do_structure_search, search_tool, e_value, query_id, coverage, res_cutoff,
                 output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db, verbose_level):
    # check validity of inputted parameters
    #
    # input input_filepath: str, do_structure_search: bool, search_tool: str, e_value: float, query_id: float,
    # coverage: float, res_cutoff: float, output_directory: str, blast_bin: str/None, blast_db: str,
    # hhsearch_bin: str/None, hhsearch_db: str, verbose_level: int
    # return inputs_valid: bool

    # check whether an inputfile with the extension .sqcs is specified
    if input_filepath and input_filepath.endswith(".sqcs"):
        # check whether inputted search tool is either hhsearch or blastp
        if search_tool in ["blastp", "hhsearch"]:
            # check whether value given for e-value can be turned into a float variable
            try:
                e_value = float(e_value)
                # check whether e-value is a value between 0 and 1
                if 0 < e_value < 1:
                    # check whether value given for query identity can be turned into a float variable
                    try:
                        query_id = float(query_id)
                        # check whether query identity is a value in [0,100]
                        if 0 <= query_id <= 100:
                            # check whether value given for coverage can be turned into a float variable
                            try:
                                coverage = float(coverage)
                                # check whether coverage is a value in [0,100]
                                if 0 <= coverage <= 100:
                                    try:
                                        # check whether value given for resolution cutoff can be turned into a float
                                        # variable
                                        res_cutoff = float(res_cutoff)
                                        return True
                                    except:
                                        print(f"Error! Given value for resolution cutoff could not be turned into a"
                                              f" float variable (given: {res_cutoff}")
                                else:
                                    print(f"Error! Given coverage is not within allowed interval of [0,100] "
                                          f"(given: {coverage}).")
                            except:
                                print(f"Error! Given value for coverage could not be turned into a float variable "
                                      f"(given: {coverage}")
                        else:
                            print(f"Error! Given query identity is not within allowed interval of [0,100] "
                                  f"(given: {query_id}).")
                    except:
                        print(f"Error! Given value for query identity could not be turned into a float variable "
                              f"(given: {query_id}")
                else:
                    print(f"Error! Given e-value is not within allowed interval of (0,1) (given: {e_value}).")
            except:
                print(f"Error! Given value for e-value could not be turned into a float variable (given: {e_value}")
        else:
            print(f"Error! Invalid search tool selected (given: {search_tool}).")
    else:
        print(f"Error! The parameter \"input-filepath\" was not given or was not ending with the '.sqcs' extension "
              f"(given: {input_filepath}).")
    return False
