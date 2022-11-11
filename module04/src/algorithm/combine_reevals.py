import sys
import pandas as pd


def combine_inter_reevaluations(data, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness,
                                distance_maximum, cutoff):
    # combine distance and homo signal reevaluation, create score that represents inter interaction affinity (higher
    # value, higher affinity)
    #
    # input data: pd.DataFrame, plddt_cutoff: float, linker_minimum: float, linker_maximum: float,
    # euclidean_strictness: float, distance_maximum: float, cutoff: float
    # return data: pd.DataFrame

    data["inter_score"] = data.apply(lambda x: score_inter_potential(x, plddt_cutoff, linker_minimum, linker_maximum,
                                                                     euclidean_strictness, distance_maximum), axis=1)
    data["final_XL_type"] = data.apply(lambda x: "inter" if x.inter_score > cutoff else "intra", axis=1)
    return data


def score_inter_potential(datapoint, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness,
                          distance_maximum):
    # combine distance and homo signal reevaluation, create score that represents inter interaction affinity (higher
    # value, higher affinity)
    #
    # input datapoint: pd.Series, plddt_cutoff: float, linker_minimum: float, linker_maximum: float,
    # euclidean_strictness: float, distance_maximum: float
    # return score: float

    score = 0.0
    euclidean_linker_minimum = linker_minimum - euclidean_strictness
    euclidean_linker_minimum = euclidean_linker_minimum if euclidean_linker_minimum > 0 else 0
    euclidean_linker_maximum = linker_maximum - euclidean_strictness
    euclidean_linker_maximum = euclidean_linker_maximum if euclidean_linker_maximum > 0 else 0

    # calculate distance argument of inter protein crosslink confidence score, if plddt included, which is True if entry
    # is not retrieved from AlphaFold or if both plddts surpass or are equal to the given cutoff
    dist_argument = (not pd.isna(datapoint.eucl_dist_tplk)) and \
                    (not pd.isna(datapoint.topo_dist_tplk)) and \
                    (datapoint.pLDDT_a == '-' or float(datapoint.pLDDT_a) >= plddt_cutoff) and \
                    (datapoint.pLDDT_b == '-' or float(datapoint.pLDDT_b) >= plddt_cutoff)

    if dist_argument:
        # include euclidean distance score only if euclidean strictness is not None
        if euclidean_strictness is not None:
            # calculate euclidean distance score
            eucl_dist = datapoint.eucl_dist_tplk if datapoint.eucl_dist_tplk < distance_maximum else distance_maximum
            raw_eucl_score = (eucl_dist - euclidean_linker_maximum) / (distance_maximum - euclidean_linker_maximum)
            if datapoint.eucl_dist_tplk <= euclidean_linker_minimum:
                score += .25
            elif datapoint.eucl_dist_tplk >= euclidean_linker_maximum:
                score += raw_eucl_score * .25

        # calculate topological distance score
        topo_dist = datapoint.topo_dist_tplk if datapoint.topo_dist_tplk < distance_maximum else distance_maximum
        raw_topo_score = (topo_dist - linker_maximum) / (distance_maximum - linker_maximum)
        if datapoint.topo_dist_tplk <= linker_minimum:
            # If the euclidean strictness is set to None the max score for topological is 0.5, else 0.25
            score += .25 if euclidean_strictness is not None else .5
        elif datapoint.topo_dist_tplk >= linker_maximum:
            # If the euclidean strictness is set to None the max score for topological is 0.5, else 0.25
            score += raw_topo_score * .25 if euclidean_strictness is not None else raw_topo_score * .5

    # calculate OPS argument of inter protein crosslink confidence score
    if datapoint.homo_pep_overl:
        score += datapoint.homo_int_overl / 2

    return score
