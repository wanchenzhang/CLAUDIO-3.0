import click
import os
import sys
import ast
import pandas as pd

from module02.src01_uniprot_search.dict.default_projections import *


@click.command()
@click.option("-i", "--input-filepath", default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv")
@click.option("-p", "--projections", default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant))
@click.option("-rt", "--read-temps", default=False)
@click.option("-t", "--search-tool", default="blastp")
@click.option("-pc", "--plddt-cutoff", default=70.0)
@click.option("-e", "--e-value", default=1e-5)
@click.option("-qi", "--query-id", default=90.0)
@click.option("-c", "--coverage", default=50.0)
@click.option("-r", "--res-cutoff", default=6.5)
@click.option("-o", "--output-directory", default="data/out/module02")
def main(input_filepath, projections, read_temps, search_tool, plddt_cutoff, e_value, query_id, coverage,
         res_cutoff, output_directory):

    if not output_directory.endswith('/'):
        output_directory += '/'

    command = f"python3 claudio_mod02_01.py -i {input_filepath} -p \"{projections}\" -s {not read_temps} " \
              f"-o {output_directory}"
    command += " && "

    filename = input_filepath.split('/')[-1]
    if not os.path.exists(f"{output_directory}structures"):
        os.mkdir(f"{output_directory}structures")
    command += f"python3 claudio_mod02_02.py -i {f'{output_directory}{filename}.sqcs'} " \
               f"-s {not read_temps} -t {search_tool} -e {e_value} -q {query_id} -c {coverage} " \
               f"-r {res_cutoff} -o {f'{output_directory}structures'}"
    command += " && "

    command += f"python3 claudio_mod02_03.py -i {f'{output_directory}structures'} " \
               f"-i2 {f'{output_directory}{filename}.sqcs'} -t {search_tool} -p {plddt_cutoff} " \
               f"-o {output_directory}"
    command += " && "

    command += f"rm {output_directory}{filename}.sqcs"

    os.system(command)
    sys.exit()
