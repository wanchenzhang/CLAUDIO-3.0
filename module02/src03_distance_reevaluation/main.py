import click
import time
import sys
import os

from module02.src03_distance_reevaluation.io.read_uniprot_search_out import read_unipsearch_out
from module02.src03_distance_reevaluation.algorithm.search_pdb_pos import search_site_pos_in_pdb
from module02.src03_distance_reevaluation.algorithm.calc_site_distances import calculate_site_dists
from module02.src03_distance_reevaluation.io.write_out import write_output
from module02.src03_distance_reevaluation.algorithm.create_plots import create_histogram


@click.command()
@click.option("-i", "--input-directory", default="data/out/structure_search")
@click.option("-i2", "--input-filepath", default="data/out/uniprot_search/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv.sqcs")
@click.option("-t", "--search-tool", default="blastp")
@click.option("-p", "--plddt-cutoff", default=70.0)
@click.option("-o", "--output-directory", default="data/out/dist_reeval")
@click.option("-tl", "--topolink-bin", default=None)
def main(input_directory, input_filepath, search_tool, plddt_cutoff, output_directory, topolink_bin):
    print("Start intra interaction check\n")
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
    if inputs_valid(input_directory, input_filepath, search_tool, plddt_cutoff, output_directory, topolink_bin):
        # Read result from uniprot_search, e.g. sqcs-file
        print("Read peptide information from uniprot search results")
        data = read_unipsearch_out(input_filepath)

        # Search for site positions in pdb files (replace rcsb pdb with alphafold, if not able to find it there)
        print("Search site pos in pdb files (replace rcsb-pdb with alphafold-pdb if needed)")
        data = search_site_pos_in_pdb(data)

        # Compute distances of sites, and if distance calculation successful compute new xl_type
        print("Calculate presumed interaction site distances and evaluate interaction likelihood")
        data = calculate_site_dists(data, plddt_cutoff, topolink_bin)

        # Plot histograms of distances
        print("Create distance histograms")
        create_histogram(data, input_filepath.split('/')[-1], output_directory)

        # Write outputfile as csv
        print("Write outputfile")
        write_output(data, input_filepath.split('/')[-1], output_directory)

    print(f"\nEnd script (Elapsed time: {round(time.time() - start_time, 2)}s)")
    print("===================================")
    sys.exit()


def inputs_valid(input_directory, input_filename, search_tool, plddt_cutoff, output_directory, topolink_bin):
    # check validity of inputted parameters
    #
    # input input_directory: str, input_filename: str, search_tool: str, plddt_cutoff: float, output_directory: str,
    # topolink_bin: str/None
    # return inputs_valid: bool

    if any([".pdb" in filename for filename in os.listdir(input_directory)]):
        if input_filename.endswith(".csv.sqcs"):
            if search_tool in ["blastp", "hhsearch"]:
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
            else:
                print(f"Error! Invalid search tool! (given: {search_tool})")
        else:
            print(f"Error! Given inputfile (-i2) extension hints at incorrect output of uniprot-search "
                  f"(given: {input_filename})!")
    else:
        print(f"Error! No structure files created by structure_search tool found (given: {input_directory})!")
    return False
