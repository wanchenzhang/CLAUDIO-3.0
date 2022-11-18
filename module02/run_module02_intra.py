import click
import os
import sys

from module02.src_structure_search.main import main as run_structure_search
from module02.src_distance_reevaluation.main import main as run_distance_analysis


@click.command()
@click.option("-i", "--input-filepath", default="data/out/unique_protein_list/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.sqcs")
@click.option("-rt", "--read-temps", default=False)
@click.option("-t", "--search-tool", default="blastp")
@click.option("-pc", "--plddt-cutoff", default=70.0)
@click.option("-e", "--e-value", default=1e-5)
@click.option("-qi", "--query-id", default=90.0)
@click.option("-c", "--coverage", default=50.0)
@click.option("-r", "--res-cutoff", default=6.5)
@click.option("-o", "--output-directory", default="data/out/module02")
@click.option("-bl", "--blast-bin", default=None)
@click.option("-bldb", "--blast-db", default="$BLASTDB")
@click.option("-hh", "--hhsearch-bin", default=None)
@click.option("-hhdb", "--hhsearch-db", default="$HHDB")
@click.option("-hhout", "--hhsearch-out", default="$HHOUT")
@click.option("-tl", "--topolink-bin", default=None)
def main(input_filepath, read_temps, search_tool, plddt_cutoff, e_value, query_id, coverage, res_cutoff,
         output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db, hhsearch_out, topolink_bin):

    if not output_directory.endswith('/'):
        output_directory += '/'

    filename = input_filepath.split('/')[-1]
    try:
        run_structure_search(["-i", input_filepath, "-s", not read_temps, "-t", search_tool,
                              "-e", e_value, "-q", query_id, "-c", coverage, "-r", res_cutoff, "-o", output_directory,
                              "-bl", blast_bin, "-bldb", blast_db, "-hh", hhsearch_bin, "-hhdb", hhsearch_db, "-hhout",
                              hhsearch_out])
    except SystemExit:
        pass
    try:
        run_distance_analysis(["-i", f"{output_directory}structures",
                               "-i2", f"{output_directory}{filename}_structdi.csv", "-t", search_tool,
                               "-p", plddt_cutoff, "-o", output_directory, "-tl", topolink_bin])
    except SystemExit:
        pass

    sys.exit()
