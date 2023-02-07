import sys
import requests as r
import pandas as pd
import os
import time
import ast


def retrieve_oligomeric_states(data, intra_only):
    # retrieve oligomeric states known for each uniprot entry and add them to the dataset as a string
    #
    # input data: pd.DataFrame, intra_only: bool
    # return data: pd.DataFrame

    # container for already searched oligo-states
    known_ostates = {}

    data["swiss_model_homology"] = data.apply(lambda x: get_oligo_state_from_swiss(x, intra_only, known_ostates),
                                              axis=1)
    return data


def get_oligo_state_from_swiss(data, intra_only, known_ostates):
    # access SWISS-MODEL for given datapoint's uniprot id, if not previously encountered, else retrieve known result
    # from known_unips
    #
    # input data: pd.Series, intra_only: bool, known_ostates: dict(str: list(str)))
    # return oligo_states: str

    NUMBER_OF_CALL_REPEATS = 5
    DOWNLOAD_RATE_LIMITER_IN_SECONDS = .05

    base_url = "https://swissmodel.expasy.org/repository/uniprot/"

    unip_ids = [data['unip_id']] if intra_only else [data['unip_id_a'], data['unip_id_b']]

    # if uniprot id not in already searched entries, do search in SWISS-MODEL
    for unip_id in unip_ids:
        if unip_id not in known_ostates.keys():
            # set up json url with uniprot id
            url = f"{base_url}{unip_id}.json"
            # repeat SWISS-MODEL calls for consistency (SWISS-MODEL has shown to inconsistently return empty or only
            # partial API call results)
            ostates = []
            for _ in range(NUMBER_OF_CALL_REPEATS):
                list_of_states = {}
                # use get-request, retrieve all known multimer complexes, isolate homomeric oligomer-states into list
                try:
                    try:
                        list_of_states = {structure["template"].split('.')[0]: structure["oligo-state"]
                                          for structure in ast.literal_eval(r.get(url).text)["result"]["structures"]}
                        # add time skip to limit the rate of calls per second
                        # (see: https://swissmodel.expasy.org/docs/help, Modelling API)
                        time.sleep(DOWNLOAD_RATE_LIMITER_IN_SECONDS)
                    except KeyError as e:
                        print("Warning! Received 'Exceeding rate limit'-error from SWISS API. "
                              "Download speed will be reduced.")
                        DOWNLOAD_RATE_LIMITER_IN_SECONDS += .1
                except r.exceptions.Timeout:
                    pass
                except ConnectionError as e:
                    print("No connection to SWISS-MODEL API possible. Please try again later.")
                    print(e)
                    sys.exit()
                ostates.append(list_of_states)

            # Add resulting oligomeric states to list of known, if results not empty
            if ostates:
                ostates = {key: pd.unique([repeat[key] for repeat in ostates]).tolist() for key in ostates[0].keys()}

                # add result to known oligomeric states
                known_ostates[unip_id] = ostates

    # return homo-oligomer states if intra crosslink
    if intra_only or data['unip_id_a'] == data['unip_id_b']:
        unique_ostates = sorted(pd.unique([state for _, states in
                                           known_ostates[data['unip_id'] if intra_only else data['unip_id_a']].items()
                                           for state in states]).tolist())
        unique_ostates = [state.replace('-', '') for state in unique_ostates if state not in ["heteromer", "monomer"]]
        return '_'.join(unique_ostates)
    # else, compute intersecting set of structures and return their oligomeric states
    else:
        intsect_oligo_states = {key: value for key, value in known_ostates[data['unip_id_a']].items()
                                if key in known_ostates[data['unip_id_b']].keys()}
        unique_intersect_ostates = \
            sorted(pd.unique([state for _, states in intsect_oligo_states.items() for state in states]).tolist())
        unique_intersect_ostates = [state.replace('-', '') for state in unique_intersect_ostates
                                    if (state != "monomer") and not state.startswith("homo")]
        return '_'.join(unique_intersect_ostates)
