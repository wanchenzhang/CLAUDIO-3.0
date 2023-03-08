import click
import time
import sys
import os

from module02.src_distance_reevaluation.io.read_uniprot_search_out import read_unipsearch_out
from module02.src_distance_reevaluation.algorithm.search_pdb_pos import search_site_pos_in_pdb
from module02.src_distance_reevaluation.algorithm.calc_site_distances import calculate_site_dists
from module02.src_distance_reevaluation.io.write_out import write_output
from module02.src_distance_reevaluation.algorithm.create_plots import create_histogram

from utils.utils import *


@click.command()
@click.option("-i", "--input-directory", default="data/out/structure_search")
@click.option("-i2", "--input-filepath", default="data/out/unique_protein_list/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.sqcs_structdi.csv")
@click.option("-t", "--search-tool", default="blastp")
@click.option("-x", "--xl-residues", default="K,M:N:1")
@click.option("-p", "--plddt-cutoff", default=70.0)
@click.option("-o", "--output-directory", default="data/out/dist_reeval")
@click.option("-tl", "--topolink-bin", default=None)
@click.option("-v", "--verbose-level", default=3)
def main(input_directory, input_filepath, search_tool, xl_residues, plddt_cutoff, output_directory, topolink_bin,
         verbose_level):
    verbose_print("Start intra interaction check", 0, verbose_level)
    start_time = time.time()

    output_directory = output_directory if output_directory else '/'.join(input_filepath.split('/')[:-1])

    # Convert directory paths to literals if None
    if topolink_bin == "None":
        topolink_bin = None

    # Add '/' to end of directory paths if not there
    if not input_directory.endswith('/'):
        input_directory += '/'
    if not output_directory.endswith('/'):
        output_directory += '/'
    if (topolink_bin is not None) and (not topolink_bin.endswith('/')):
        topolink_bin += '/'

    # If parameters inputted by user valid
    if inputs_valid(input_directory, input_filepath, search_tool, xl_residues, plddt_cutoff, output_directory,
                    topolink_bin, verbose_level):
        # Define dataset for crosslink residues including possible positions and atom types
        df_xl_res = build_xl_dataset(xl_residues)

        # Read result from uniprot_search, e.g. sqcs-file
        verbose_print("Read peptide information from uniprot search results", 0, verbose_level)
        data, intra_only = read_unipsearch_out(input_filepath)

        # Search for site positions in pdb files (replace rcsb pdb with alphafold, if not able to find it there)
        verbose_print("Search site pos in pdb files (replace rcsb-pdb with alphafold-pdb if needed)", 0, verbose_level)
        data = search_site_pos_in_pdb(data, df_xl_res, intra_only, verbose_level)

        # Compute distances of sites, and if distance calculation successful compute new xl_type
        verbose_print("Calculate presumed interaction site distances and evaluate interaction likelihood", 0,
                      verbose_level)
        data = calculate_site_dists(data, df_xl_res, plddt_cutoff, intra_only, topolink_bin, verbose_level)

        # Plot histograms of distances
        verbose_print("Create distance histograms", 0, verbose_level)
        create_histogram(data, input_filepath.split('/')[-1], output_directory)

        # Drop temporary result columns
        data.drop(["path", "pdb_method", "pdb_resolution", "eucl_dist", "res_criteria_fulfilled", "res_crit_a",
                   "res_crit_b", "method_a", "method_b"], axis=1, inplace=True)
        # Rename certain result columns
        data.rename(columns={"best_res_pdb_method": "pdb_method", "best_res_pdb_resolution": "pdb_resolution",
                             "eucl_dist_tplk": "eucl_dist", "topo_dist_tplk": "topo_dist"}, inplace=True)

        # Overwrite previous outputfile of module02
        verbose_print("Overwrite outputfile", 0, verbose_level)
        write_output(data, input_filepath)

    verbose_print(f"\nEnd script (Elapsed time: {round_self(time.time() - start_time, 2)}s)", 0, verbose_level)
    verbose_print("===================================", 0, verbose_level)
    sys.exit()


def inputs_valid(input_directory, input_filename, search_tool, xl_residues, plddt_cutoff, output_directory,
                 topolink_bin, verbose_level):
    # check validity of inputted parameters
    #
    # input input_directory: str, input_filename: str, search_tool: str, xl_residues: str, plddt_cutoff: float,
    # output_directory: str, topolink_bin: str/None, verbose_level: int
    # return inputs_valid: bool

    if any([".pdb" in filename for filename in os.listdir(input_directory)]):
        if input_filename.endswith(".sqcs_structdi.csv"):
            if search_tool in ["blastp", "hhsearch"]:
                # check whether xl_residues can be turned into a proper DataFrame, else return False
                try:
                    if build_xl_dataset(xl_residues) is None:
                        return False

                    # check whether plddt cutoff has valid value
                    try:
                        plddt_cutoff = float(plddt_cutoff)
                        if 0 <= plddt_cutoff <= 100:
                            return True
                        else:
                            print(f"pLDDT cutoff value should be in [0, 100] (given: {plddt_cutoff}).")
                    except:
                        print(
                            f"Value given for pLDDT cutoff should be possible to turn into a float "
                            f"(given: {plddt_cutoff}).")
                except:
                    print(f"Error! Could not properly parse xl_residues for accepted crosslinked residues "
                          f"(given: {xl_residues}).")
            else:
                print(f"Error! Invalid search tool! (given: {search_tool})")
        else:
            print(f"Error! Given inputfile (-i2) extension hints at incorrect output of uniprot-search "
                  f"(given: {input_filename})!")
    else:
        print(f"Error! No structure files created by structure_search tool found (given: {input_directory})!")
    return False
