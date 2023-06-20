import click
import sys
import time
import ast

from module04.src.io.read_ins import read_inputs
from module04.src.algorithm.combine_reevals import combine_inter_reevaluations
from module04.src.algorithm.retrieve_oligo_state import retrieve_oligomeric_states
from module04.src.algorithm.create_histo import create_histograms
from module04.src.io.write_outs import write_outputs

from utils.utils import *


@click.command()
@click.option("-i", "--input-filepath", default="data/out/dist_reeval/sample_data.sqcs_structdi.csv")
@click.option("-i2", "--input-filepath2", default="data/out/homo_signal/sample_data.sqcs_ops.csv")
@click.option("-p", "--plddt-cutoff", default=70.0)
@click.option("-lmin", "--linker-minimum", default=5.0)
@click.option("-lmax", "--linker-maximum", default=35.0)
@click.option("-es", "--euclidean-strictness", default=None)
@click.option("-dm", "--distance-maximum", default=50.0)
@click.option("-c", "--cutoff", default=0.0)
@click.option("-o", "--output-directory", default="data/out/new_inter/")
@click.option("-s", "--compute-scoring", default=False)
@click.option("-v", "--verbose-level", default=3)
def main(input_filepath, input_filepath2, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness,
         distance_maximum, cutoff, output_directory, compute_scoring, verbose_level):
    verbose_print("Start New Inter Interaction Analysis", 0, verbose_level)
    start_time = time.time()

    output_directory = output_directory if output_directory else '/'.join(input_filepath.split('/')[:-1])
    if not output_directory.endswith('/'):
        output_directory += '/'

    # If parameters inputted by user valid
    if inputs_valid(input_filepath, input_filepath2, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness,
                    distance_maximum, cutoff, output_directory, compute_scoring, verbose_level):
        plddt_cutoff = float(plddt_cutoff)
        linker_minimum = float(linker_minimum)
        linker_maximum = float(linker_maximum)
        euclidean_strictness = float(euclidean_strictness) \
            if euclidean_strictness != "None" and euclidean_strictness is not None else None
        distance_maximum = float(distance_maximum)
        cutoff = float(cutoff)
        filename = input_filepath.split('/')[-1].split('.')[0]

        # Read inputs
        verbose_print("Read inputs", 0, verbose_level)
        data = read_inputs(input_filepath, input_filepath2)

        # Combine results of both reevaluations
        verbose_print("Combine results of both reevaluations", 0, verbose_level)
        data = combine_inter_reevaluations(data, plddt_cutoff, linker_minimum, linker_maximum,
                                           euclidean_strictness, distance_maximum, cutoff, compute_scoring)

        # Retrieve known oligomeric states from SWISS-MODEL
        verbose_print("Retrieve known oligomeric states from SWISS-MODEL", 0, verbose_level)
        data = retrieve_oligomeric_states(data, verbose_level)

        # Clean dataset for output
        data = clean_dataset(data)

        # Create inter score histogram
        verbose_print("Create score histogram", 0, verbose_level)
        create_histograms(data, filename, cutoff, compute_scoring, output_directory)

        # Write final csv containing all computed information, fastas for alphafold and protein-specific csv with
        # interaction restraints
        verbose_print("Write outputs", 0, verbose_level)
        write_outputs(data, filename, compute_scoring, output_directory)

    verbose_print(f"\nEnd script (Elapsed time: {round_self(time.time() - start_time, 2)}s)", 0, verbose_level)
    verbose_print("===================================", 0, verbose_level)
    sys.exit()


def inputs_valid(input_filepath, input_filepath2, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness,
                 distance_maximum, cutoff, output_directory, compute_scoring, verbose_level):
    # check validity of inputted parameters
    #
    # input input_filepath: str, input_filepath2: str, plddt_cutoff: float, linker_minimum: float,
    # linker_maximum: float, euclidean_strictness: float, distance_maximum: float, cutoff: float, output_directory: str,
    # compute_scoring: bool, verbose_level: int
    # return inputs_valid: bool

    # check whether outputfile from distance-based reevauation is specified
    if input_filepath.endswith(".sqcs_structdi.csv"):
        # check whether outputfile from homo-signal-based reevauation is specified
        if input_filepath2.endswith(".sqcs_ops.csv"):
            # check whether plddt cutoff has valid value
            try:
                plddt_cutoff = float(plddt_cutoff)
                if 0 <= plddt_cutoff <= 100:
                    # check whether minimum crosslinker distance has valid value
                    try:
                        linker_minimum = float(linker_minimum)
                        # check whether maximum crosslinker distance has valid value
                        try:
                            linker_maximum = float(linker_maximum)
                            # check whether euclidean strictness has valid value (either float or None)
                            try:
                                if euclidean_strictness is not None:
                                    euclidean_strictness = float(euclidean_strictness)
                            except ValueError:
                                try:
                                    if euclidean_strictness is not None:
                                        euclidean_strictness = ast.literal_eval(euclidean_strictness)
                                except ValueError:
                                    print(f"Error! Could not change type of given euclidean strictness to either float "
                                          f"or None (given: {euclidean_strictness}).")
                                    return False
                            # check whether maximum distance has valid value
                            try:
                                distance_maximum = float(distance_maximum)
                                # check whether reevaluation cutoff has valid value
                                try:
                                    cutoff = float(cutoff)
                                    if 0 <= cutoff <= 1:
                                        return True
                                    else:
                                        print(f"Cutoff value for reclassification should be in [0, 1] "
                                              f"(given: {cutoff}).")
                                except:
                                    print(f"Value given for reclassification cutoff should be possible to turn into a"
                                          f" float (given: {cutoff}).")
                            except:
                                print(f"Value given for maximum distance value should be possible to turn into a "
                                      f"float (given: {distance_maximum}).")
                        except:
                            print(f"Value given for crosslinker maximum should be possible to turn into a float "
                                  f"(given: {linker_maximum}).")
                    except:
                        print(f"Value given for crosslinker minimum should be possible to turn into a float "
                              f"(given: {linker_minimum}).")
                else:
                    print(f"pLDDT cutoff value should be in [0, 100] (given: {plddt_cutoff}).")
            except:
                print(f"Value given for pLDDT cutoff should be possible to turn into a float (given: {plddt_cutoff}).")
        else:
            print(f"Homo-signal reevaluation outputfile was either not correctly given or has wrong extension "
                  f"(given: {input_filepath2}).")
    else:
        print(f"Distance reevaluation outputfile was either not correctly given or has wrong extension "
              f"(given: {input_filepath}).")
    return False
