import click
import os
import time
import sys
import ast
import pandas as pd

from module02.src01_uniprot_search.dict.default_projections import *
from module01.src.main import main as run_claudio_lists
from module02.run_module02_intra import main as run_claudio_structdi
from module03.src.main import main as run_claudio_ops
from module04.src.main import main as run_claudio_xl


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

    try:
        run_claudio_lists(["-i", input_filepath, "-p", projections, "-t", search_tool, "-o", output_directory])
    except SystemExit:
        pass
    try:
        run_claudio_ops(["-i", input_filepath, "-p", projections, "-s", not read_temps, "-o", output_directory])
    except SystemExit:
        pass
    try:
        run_claudio_structdi(["-i", input_filepath, "-p", projections, "-rt", read_temps, "-t", search_tool,
                              "-pc", plddt_cutoff, "-e", e_value, "-qi", query_id, "-c", coverage, "-r", res_cutoff,
                              "-o", output_directory])
    except SystemExit:
        pass

    filename = input_filepath.split('/')[-1]
    try:
        run_claudio_xl(["-i", f"{output_directory}{filename}.sqcs.csv",
                        "-i2", f"{output_directory}{filename}_homosig.csv", "-p", plddt_cutoff, "-lmin", linker_minimum,
                        "-lmax", linker_maximum, "-es", euclidean_strictness, "-dm", distance_maximum, "-c", cutoff,
                        "-o", output_directory])
    except SystemExit:
        pass

    os.remove(f"{output_directory}{filename}.sqcs.csv")
    os.remove(f"{output_directory}{filename}_homosig.csv")

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
        line_markers = ["input_filepath=", "projections=", "read_temps=", "search_tool=", "e_value=", "query_id=",
                        "coverage=", "res_cutoff=", "plddt_cutoff=", "linker_minimum=", "linker_maximum=",
                        "euclidean_strictness=", "distance_maximum=", "cutoff=", "output_directory="]
        params = [l[len(marker):] for marker in line_markers for l in input_lines if l.startswith(marker)]

        # Check parameters
        if len(params) != len(line_markers):
            print(f"Error! Number of parameters in configuration file do not match the expected number. Check whether "
                  f"you are either missing a parameter, or one is duplicated.\n"
                  f"\tExpected: {len(line_markers)}\n"
                  f"\tReceived: {len(list(params))}")
            sys.exit()
        try:
            params[2] = params[2] == "True"
            params[4] = float(params[4])
            params[5] = float(params[5])
            params[6] = float(params[6])
            params[7] = float(params[7])
            params[8] = float(params[8])
            params[9] = float(params[9])
            params[10] = float(params[10])
            params[11] = float(params[11])
            params[12] = float(params[12])
            params[13] = float(params[13])
        except:
            print("Error! Could not change type of one or multiple parameters as intended.")
            sys.exit()
        return tuple(params)


if __name__ == "__main__":
    main()
