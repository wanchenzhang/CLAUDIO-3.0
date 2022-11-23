import sys
import requests as r
import pandas as pd

NUMBER_OF_CALL_REPEATS = 5


def retrieve_oligomeric_states(data):
    # retrieve oligomeric states known for each uniprot entry and add them to the dataset as a string
    #
    # input data: pd.DataFrame, filename: str, output_directory: str
    # return data: pd.DataFrame

    # container for already searched oligo-states
    known_oligo_states = []

    data["oligo_states"] = data.apply(lambda x: get_oligo_state_from_swiss(x, known_oligo_states), axis=1)
    return data


def get_oligo_state_from_swiss(data, known_oligo_states):
    # access SWISS-MODEL for given datapoint's uniprot id, if not previously encountered, else retrieve known result
    # from known_unips
    #
    # input data: pd.Series, known_oligo_states: list(tuple(str, str))
    # return oligo_states: str

    url = "https://swissmodel.expasy.org/repository/uniprot/"

    known_unips = [x[0] for x in known_oligo_states]

    # if uniprot id not in already searched entries, do search in SWISS-MODEL
    if not data['unip_id_a'] in known_unips:
        # set json url with uniprot id up
        url = f"{url}{data['unip_id_a']}.json"

        # repeat SWISS-MODEL calls for consistency (SWISS-MODEL has shown to inconsistently return empty API call
        # results)
        oligo_states = []
        for _ in range(NUMBER_OF_CALL_REPEATS):
            # use get-request, retrieve all known multimer complexes, isolate homomeric oligomer-states into list
            try:
                list_of_states = [l.strip().replace("oligo-state\": ", '')
                                  for l in r.get(url).text.split('\n') if "oligo" in l]
                list_of_states = sorted([state.replace('\"', '').replace(',', '').replace('-', '').replace("homo", '')
                                         for state in list_of_states if state not in ['""heteromer",', '""monomer",']])
            except ConnectionError as e:
                print("No connection to SWISS-MODEL API possible. Please try again later.")
                print(e)
                sys.exit()
            oligo_states.append(list_of_states)

        # join known homomeric states into one string
        oligo_states = pd.unique([states for repeat in oligo_states for states in repeat]).tolist()
        oligo_states = '_'.join(oligo_states)

        # add result to known oligomeric states
        known_oligo_states.append((data['unip_id_a'], oligo_states))
        return oligo_states
    # Else retrieve already known oligomeric state result from list
    else:
        return known_oligo_states[known_unips.index(data['unip_id_a'])][1]
