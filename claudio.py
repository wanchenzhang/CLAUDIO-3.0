import click
import os
import time
import sys

from module01.src.main import main as run_claudio_lists
from module02.run_module02_intra import main as run_claudio_structdi
from module03.src.main import main as run_claudio_ops
from module04.src.main import main as run_claudio_xl

from utils.utils import *

_DEFAULT_OPTIONS = ["data/in/sample_data.csv", "", "peptide1,peptide2,position1,position2,k_pos1,k_pos2,entry1,entry2",
                    False, "K,M:N:1", "blastp", 1e-5, 90.0, 50.0, 6.5, 70.0, 0.0, 35.0, None, 50.0, 0.0,
                    "data/out/full", None, "$HOME/BLAST/db", None, "$HOME/HHSUITE/db", None, False, 3]


@click.command()
@click.option("-i", "--input-filepath", default=_DEFAULT_OPTIONS[0])
@click.option("-it", "--input-temppath", default=_DEFAULT_OPTIONS[1])
@click.option("-p", "--projections", default=_DEFAULT_OPTIONS[2])
@click.option("-rt", "--read-temps", default=_DEFAULT_OPTIONS[3])
@click.option("-x", "--xl-residues", default=_DEFAULT_OPTIONS[4])
@click.option("-t", "--search-tool", default=_DEFAULT_OPTIONS[5])
@click.option("-e", "--e-value", default=_DEFAULT_OPTIONS[6])
@click.option("-qi", "--query-id", default=_DEFAULT_OPTIONS[7])
@click.option("-cv", "--coverage", default=_DEFAULT_OPTIONS[8])
@click.option("-r", "--res-cutoff", default=_DEFAULT_OPTIONS[9])
@click.option("-pc", "--plddt-cutoff", default=_DEFAULT_OPTIONS[10])
@click.option("-lmin", "--linker-minimum", default=_DEFAULT_OPTIONS[11])
@click.option("-lmax", "--linker-maximum", default=_DEFAULT_OPTIONS[12])
@click.option("-es", "--euclidean-strictness", default=_DEFAULT_OPTIONS[13])
@click.option("-dm", "--distance-maximum", default=_DEFAULT_OPTIONS[14])
@click.option("-ct", "--cutoff", default=_DEFAULT_OPTIONS[15])
@click.option("-o", "--output-directory", default=_DEFAULT_OPTIONS[16])
@click.option("-bl", "--blast-bin", default=_DEFAULT_OPTIONS[17])
@click.option("-bldb", "--blast-db", default=_DEFAULT_OPTIONS[18])
@click.option("-hh", "--hhsearch-bin", default=_DEFAULT_OPTIONS[19])
@click.option("-hhdb", "--hhsearch-db", default=_DEFAULT_OPTIONS[20])
@click.option("-tl", "--topolink-bin", default=_DEFAULT_OPTIONS[21])
@click.option("-s", "--compute-scoring", default=_DEFAULT_OPTIONS[22])
@click.option("-v", "--verbose-level", default=_DEFAULT_OPTIONS[23])
@click.option("-c", "--config", default="")
def main(input_filepath, input_temppath, projections, read_temps, xl_residues, search_tool, e_value, query_id, coverage,
         res_cutoff, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness, distance_maximum, cutoff,
         output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db, topolink_bin, compute_scoring, verbose_level,
         config):
    verbose_print("Start full CLAUDIO pipeline", 0, verbose_level)
    verbose_print("===================================", 0, verbose_level)
    start_time = time.time()

    # If configuration file given, ignore(/overwrite) all other parameters
    if config:
        verbose_print("Configuration file given", 0, verbose_level)
        args = [input_filepath, input_temppath, projections, read_temps, xl_residues, search_tool, e_value, query_id,
                coverage, res_cutoff, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness,
                distance_maximum, cutoff, output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db,
                topolink_bin, compute_scoring, verbose_level]
        input_filepath, input_temppath, projections, read_temps, xl_residues, search_tool, e_value, query_id, \
            coverage, res_cutoff, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness, \
            distance_maximum, cutoff, output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db, topolink_bin, \
            compute_scoring, verbose_level \
            = read_config(config, args)

    if not output_directory.endswith('/'):
        output_directory += '/'
    filename = '.'.join(input_filepath.split('/')[-1].split('.')[:-1])
    verbose_level = int(verbose_level)

    # Run Module01
    try:
        run_claudio_lists(["-i", input_filepath, "-it", input_temppath, "-p", projections, "-s", not read_temps,
                           "-x", xl_residues, "-t", search_tool, "-o", output_directory, "-bl", blast_bin,
                           "-bldb", blast_db, "-hh", hhsearch_bin, "-hhdb", hhsearch_db, "-v", verbose_level])
    except SystemExit:
        pass
    if not os.path.exists(f"{output_directory}{filename}.sqcs"):
        sys.exit()

    # Run Module03
    try:
        run_claudio_ops(["-i", f"{output_directory}{filename}.sqcs", "-o", output_directory, "-v", verbose_level])
    except SystemExit:
        pass

    # Run Module02
    try:
        run_claudio_structdi(["-i", f"{output_directory}{filename}.sqcs", "-it", input_temppath, "-rt", read_temps,
                              "-t", search_tool, "-x", xl_residues, "-pc", plddt_cutoff, "-lmin", linker_minimum,
                              "-lmax", linker_maximum, "-e", e_value, "-qi", query_id, "-c", coverage, "-r", res_cutoff,
                              "-o", output_directory, "-bl", blast_bin, "-bldb", blast_db, "-hh", hhsearch_bin,
                              "-hhdb", hhsearch_db, "-tl", topolink_bin, "-v", verbose_level])
    except SystemExit:
        pass

    # Run Module04
    try:
        run_claudio_xl(["-i", f"{output_directory}{filename}.sqcs_structdi.csv",
                        "-i2", f"{output_directory}{filename}.sqcs_ops.csv", "-p", plddt_cutoff, "-lmin", linker_minimum,
                        "-lmax", linker_maximum, "-es", euclidean_strictness, "-dm", distance_maximum, "-c", cutoff,
                        "-o", output_directory, "-s", compute_scoring, "-v", verbose_level])
    except SystemExit:
        if os.path.exists(f"{output_directory}{filename}_final.csv"):
            os.remove(f"{output_directory}{filename}.sqcs")
            os.remove(f"{output_directory}{filename}.sqcs_structdi.csv")
            os.remove(f"{output_directory}{filename}.sqcs_ops.csv")
        pass

    verbose_print(f"\nEnd full CLAUDIO pipeline execution (Total elapsed time: "
                  f"{round_self(time.time() - start_time, 2)}s)",
                  0, verbose_level)
    verbose_print("===================================", 0, verbose_level)
    sys.exit()


def read_config(path, args):
    # read configuration file from path, compare already set values to defaults thus checking whether params were given
    # additionally to the config-file and set values for other input parameters
    #
    # input path: str, args: list(str)
    # return input_filepath: str, input_temppath: str, projections: str, read_temps: bool, xl_residues: str,
    # search_tool: str, e_value: float, query_id: float, coverage: float, res_cutoff: float, plddt_cutoff: float,
    # linker_minimum: float, linker_maximum: float, euclidean_strictness: float/None, distance_maximum: float,
    # cutoff: float, output_directory: str, blast_bin: str/None, blast_db: str, hhsearch_bin: str/None,
    # hhsearch_db: str, topolink_bin: str, compute_scoring: bool, verbose_level: int

    with open(path, 'r') as f:
        config_content = f.read()
        input_lines = [l.replace('"', "'").replace(" ", "")
                       for l in config_content.split('\n') if l and not l.startswith('#')]
        line_markers = ["input_filepath=", "input_temppath=", "projections=", "read_temps=", "xl_residues=",
                        "search_tool=", "e_value=", "query_id=", "coverage=", "res_cutoff=", "plddt_cutoff=",
                        "linker_minimum=", "linker_maximum=", "euclidean_strictness=", "distance_maximum=", "cutoff=",
                        "output_directory=", "blast_bin=", "blast_db=", "hhsearch_bin=", "hhsearch_db=",
                        "topolink_bin=", "compute_scoring=", "verbose_level="]
        config_params = {marker: l[len(marker):] for marker in line_markers for l in input_lines
                         if l.startswith(marker)}

        # Check number of not empty parameters, if number is as expected convert boolean params from string into boolean
        if len([x for x in config_params.values() if x]) != len(line_markers):
            print(f"Error! Number of parameters in configuration file do not match the expected number. Check whether "
                  f"you are either missing a parameter, or one is duplicated.\n"
                  f"\tExpected number: {len(line_markers)}\n"
                  f"\tReceived number: {len([x for x in config_params.values() if x])}")
            sys.exit()
        else:
            # Check whether boolean params can be correctly converted
            if config_params["read_temps="] == "True":
                config_params["read_temps="] = True
            elif config_params["read_temps="] == "False":
                config_params["read_temps="] = False
            else:
                raise ValueError(f"Error! Could not change type of read_temp to boolean (given:{config_params[2]}).")

            if config_params["compute_scoring="] == "True":
                config_params["compute_scoring="] = True
            elif config_params["compute_scoring="] == "False":
                config_params["compute_scoring="] = False
            else:
                raise ValueError(f"Error! Could not change type of compute_scoring to boolean "
                                 f"(given:{config_params[-2]}).")
            # Set input_temppath to empty string, if given None
            if config_params["input_temppath="] == "None":
                config_params["input_temppath="] = ""

        # Check whether already given args have default value
        defaults = _DEFAULT_OPTIONS.copy()

        # define dicts using line markers as keys, for arguments and default values
        defaults = {marker: defaults[i] for i, marker in enumerate(line_markers)}
        args = {marker: args[i] for i, marker in enumerate(line_markers)}

        # define dict containing booleans for each line marker, if already given argument has default value
        arg_is_default = {marker: args[marker] == defaults[marker] for marker in line_markers}

        # define final params, if already given argument is not default, it takes precedence over param value in
        # config-file
        final_params = [config_params[marker] if arg_is_default[marker] else args[marker]
                        for marker in line_markers]

        return tuple(final_params)


if __name__ == "__main__":
    main()
