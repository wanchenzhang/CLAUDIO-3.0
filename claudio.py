import click
import os
import time
import sys
import ast
import pandas as pd

from module02.src01_uniprot_search.dict.default_projections import *


@click.command()
@click.option("-i", "--input-filepath", default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv")
@click.option("-p", "--projections", default=str(liu18_schweppe17_linked_residues_intra_homo_2672_nonredundant))
@click.option("-rt", "--read-temps", default=False)
@click.option("-t", "--search-tool", default="blastp")
@click.option("-e", "--e-value", default=1e-5)
@click.option("-qi", "--query-id", default=90.0)
@click.option("-cv", "--coverage", default=50.0)
@click.option("-r", "--res-cutoff", default=6.5)
@click.option("-pc", "--plddt-cutoff", default=70.0)
@click.option("-lmin", "--linker-minimum", default=0.0)
@click.option("-lmax", "--linker-maximum", default=35.0)
@click.option("-es", "--euclidean-strictness", default=5.0)
@click.option("-dm", "--distance-maximum", default=50.0)
@click.option("-ct", "--cutoff", default=0.0)
@click.option("-o", "--output-directory", default="data/out/full")
@click.option("-c", "--config", default='')
def main(input_filepath, projections, read_temps, search_tool, e_value, query_id, coverage, res_cutoff, plddt_cutoff,
         linker_minimum, linker_maximum, euclidean_strictness, distance_maximum, cutoff, output_directory, config):
    print("Start full CLAUDIO pipeline")
    print("===================================")
    start_time = time.time()

    # If configuration file given, ignore(/overwrite) all other parameters
    if config:
        print("Configuration file given, ignore other inputs")
        input_filepath, projections, read_temps, search_tool, e_value, query_id, coverage, res_cutoff, plddt_cutoff,\
            linker_minimum, linker_maximum, euclidean_strictness, distance_maximum, cutoff, \
            output_directory = read_config(config)

    if not output_directory.endswith('/'):
        output_directory += '/'

    command = f"python3 claudio_lists.py -i {input_filepath} -p \"{projections}\" " \
              f"-t {search_tool} -o {output_directory}"
    command += " && "

    command += f"python3 claudio_ops.py -i {input_filepath} -p \"{projections}\" "\
               f"-s {not read_temps} -o {output_directory}"
    command += " && "

    command += f"python3 claudio_structdi.py -i {input_filepath} -p \"{projections}\" -rt {read_temps} " \
               f"-t {search_tool} -pc {plddt_cutoff} -e {e_value} -qi {query_id} -c {coverage} -r {res_cutoff} " \
               f"-o {output_directory}"
    command += " && "

    filename = input_filepath.split('/')[-1]
    command += f"python3 claudio_xl.py -i {f'{output_directory}{filename}.sqcs.csv'} " \
               f"-i2 {f'{output_directory}{filename}_homosig.csv'} -p {plddt_cutoff} -lmin {linker_minimum} " \
               f"-lmax {linker_maximum} -es {euclidean_strictness} -dm {distance_maximum} -c {cutoff} " \
               f"-o {output_directory}"
    command += " && "

    command += f"rm {output_directory}{filename}.sqcs.csv {output_directory}{filename}_homosig.csv"

    os.system(command)

    print(f"\nEnd full CLAUDIO pipeline execution (Total elapsed time: {round(time.time() - start_time, 2)}s)")
    print("===================================")
    sys.exit()


def read_config(path):
    # read configuration file and set values for other input parameters
    #
    # input path: str
    # return input_filepath: str, projections: str, read_temps: bool, search_tool: str, e_value: float,
    # query_id: float, coverage: float, res_cutoff: float, plddt_cutoff: float, linker_minimum: float,
    # linker_maximum: float, euclidean_strictness: float, distance_maximum: float, cutoff: float, output_directory: str

    with open(path, 'r') as f:
        config_content = f.read()
        input_lines = [l.replace('"', "'").replace(" ", "")
                       for l in config_content.split('\n') if l and not l.startswith('#')]

        return input_lines[0], input_lines[1], input_lines[2] == "True", input_lines[3], \
            float(input_lines[4]), float(input_lines[5]), float(input_lines[6]), float(input_lines[7]), \
            float(input_lines[8]), float(input_lines[9]), float(input_lines[10]), float(input_lines[11]),\
            float(input_lines[12]), float(input_lines[13]), input_lines[14]


if __name__ == "__main__":
    main()
