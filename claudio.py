import click
import os
import time
import sys

from module01.src.main import main as run_claudio_lists
from module02.run_module02_intra import main as run_claudio_structdi
from module03.src.main import main as run_claudio_ops
from module04.src.main import main as run_claudio_xl


@click.command()
@click.option("-i", "--input-filepath", default="data/in/liu18_schweppe17_linked_residues_intra-homo_2370_nonredundant.csv")
@click.option("-p", "--projections", default="peptide1,peptide2,position1,position2,k_pos1,k_pos2,entry1,entry2")
@click.option("-rt", "--read-temps", default=False)
@click.option("-x", "--xl-residues", default="K,M:1")
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
@click.option("-bl", "--blast-bin", default=None)
@click.option("-bldb", "--blast-db", default="$BLASTDB")
@click.option("-hh", "--hhsearch-bin", default=None)
@click.option("-hhdb", "--hhsearch-db", default="$HHDB")
@click.option("-hhout", "--hhsearch-out", default="$HHOUT")
@click.option("-tl", "--topolink-bin", default=None)
@click.option("-c", "--config", default='')
def main(input_filepath, projections, read_temps, xl_residues, search_tool, e_value, query_id, coverage, res_cutoff,
         plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness, distance_maximum, cutoff, output_directory,
         blast_bin, blast_db, hhsearch_bin, hhsearch_db, hhsearch_out, topolink_bin, config):
    print("Start full CLAUDIO pipeline")
    print("===================================")
    start_time = time.time()

    # If configuration file given, ignore(/overwrite) all other parameters
    if config:
        print("Configuration file given, ignore other inputs")
        input_filepath, projections, read_temps, search_tool, e_value, query_id, coverage, res_cutoff, plddt_cutoff,\
            linker_minimum, linker_maximum, euclidean_strictness, distance_maximum, cutoff, \
            output_directory, blast_bin, blast_db, hhsearch_bin, hhsearch_db, hhsearch_out, topolink_bin \
            = read_config(config)

    if not output_directory.endswith('/'):
        output_directory += '/'
    filename = '.'.join(input_filepath.split('/')[-1].split('.')[:-1])

    # Run Module01
    try:
        run_claudio_lists(["-i", input_filepath, "-p", projections, "-s", not read_temps, "-x", xl_residues,
                           "-t", search_tool, "-o", output_directory, "-bl", blast_bin, "-bldb", blast_db,
                           "-hh", hhsearch_bin, "-hhdb", hhsearch_db, "-hhout", hhsearch_out])
    except SystemExit:
        pass

    # Run Module03
    try:
        run_claudio_ops(["-i", f"{output_directory}{filename}.sqcs", "-o", output_directory])
    except SystemExit:
        pass

    # Run Module02
    try:
        run_claudio_structdi(["-i", f"{output_directory}{filename}.sqcs", "-rt", read_temps, "-t", search_tool,
                              "-pc", plddt_cutoff, "-e", e_value, "-qi", query_id, "-c", coverage, "-r", res_cutoff,
                              "-o", output_directory, "-bl", blast_bin, "-bldb", blast_db, "-hh", hhsearch_bin,
                              "-hhdb", hhsearch_db, "-hhout", hhsearch_out, "-tl", topolink_bin])
    except SystemExit:
        pass

    # Run Module04
    try:
        run_claudio_xl(["-i", f"{output_directory}{filename}.sqcs_structdi.csv",
                        "-i2", f"{output_directory}{filename}.sqcs_ops.csv", "-p", plddt_cutoff, "-lmin", linker_minimum,
                        "-lmax", linker_maximum, "-es", euclidean_strictness, "-dm", distance_maximum, "-c", cutoff,
                        "-o", output_directory])
    except SystemExit:
        pass

    os.remove(f"{output_directory}{filename}.sqcs")
    os.remove(f"{output_directory}{filename}.sqcs_structdi.csv")
    os.remove(f"{output_directory}{filename}.sqcs_ops.csv")

    print(f"\nEnd full CLAUDIO pipeline execution (Total elapsed time: {round(time.time() - start_time, 2)}s)")
    print("===================================")
    sys.exit()


def read_config(path):
    # read configuration file and set values for other input parameters
    #
    # input path: str
    # return input_filepath: str, projections: str, read_temps: bool, search_tool: str, e_value: float,
    # query_id: float, coverage: float, res_cutoff: float, plddt_cutoff: float, linker_minimum: float,
    # linker_maximum: float, euclidean_strictness: float/None, distance_maximum: float, cutoff: float,
    # output_directory: str, blast_bin: str/None, blast_db: str, hhsearch_bin: str/None, hhsearch_db: str,
    # hhsearch_out: str, topolink_bin: str

    with open(path, 'r') as f:
        config_content = f.read()
        input_lines = [l.replace('"', "'").replace(" ", "")
                       for l in config_content.split('\n') if l and not l.startswith('#')]
        line_markers = ["input_filepath=", "projections=", "read_temps=", "search_tool=", "e_value=", "query_id=",
                        "coverage=", "res_cutoff=", "plddt_cutoff=", "linker_minimum=", "linker_maximum=",
                        "euclidean_strictness=", "distance_maximum=", "cutoff=", "output_directory=",
                        "blast_bin=", "blast_db=", "hhsearch_bin=", "hhsearch_db=", "hhsearch_out=", "topolink_bin="]
        params = [l[len(marker):] for marker in line_markers for l in input_lines if l.startswith(marker)]

        # Check number of parameters
        if len(params) != len(line_markers):
            print(f"Error! Number of parameters in configuration file do not match the expected number. Check whether "
                  f"you are either missing a parameter, or one is duplicated.\n"
                  f"\tExpected number: {len(line_markers)}\n"
                  f"\tReceived number: {len(list(params))}")
            sys.exit()

        # Check whether read_temp can be correctly converted too boolean
        try:
            if params[2] == "True":
                params[2] = True
            elif params[2] == "False":
                params[2] = False
            else:
                raise ValueError
        except ValueError:
            print(f"Error! Could not change type of read_temp to boolean (given:{params[2]}).")
            sys.exit()

        return tuple(params)


if __name__ == "__main__":
    main()
