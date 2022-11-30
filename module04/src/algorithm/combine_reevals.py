import sys
import pandas as pd


def combine_inter_reevaluations(data, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness,
                                distance_maximum, cutoff):
    # combine distance and homo signal reevaluation, create score that represents inter interaction affinity (higher
    # value, higher affinity), and evidence for this type evaluation. Base new crosslink types on each individually.
    #
    # input data: pd.DataFrame, plddt_cutoff: float, linker_minimum: float, linker_maximum: float,
    # euclidean_strictness: float, distance_maximum: float, cutoff: float
    # return data: pd.DataFrame

    # new crosslink type based on score
    data["inter_score"] = data.apply(lambda x: score_inter_potential(x, plddt_cutoff, linker_minimum, linker_maximum,
                                                                     euclidean_strictness, distance_maximum), axis=1)
    data["score_XL_type"] = data.apply(lambda x: "inter" if x.inter_score > cutoff else "intra", axis=1)

    # new crosslink type based on evidence
    data["evidence"] = data.apply(lambda x: write_evidence(x, plddt_cutoff, linker_minimum, linker_maximum,
                                                           euclidean_strictness), axis=1)
    data["final_XL_type"] = data.apply(lambda x: "inter" if x.evidence else "intra", axis=1)
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

    # calculate distance argument of inter protein crosslink confidence score, if plddt included, which is True if entry
    # is not retrieved from AlphaFold or if both plddts surpass or are equal to the given cutoff
    dist_argument = (not pd.isna(datapoint.eucl_dist_tplk)) and \
                    (not pd.isna(datapoint.topo_dist_tplk)) and \
                    (datapoint.pLDDT_a == '-' or float(datapoint.pLDDT_a) >= plddt_cutoff) and \
                    (datapoint.pLDDT_b == '-' or float(datapoint.pLDDT_b) >= plddt_cutoff)

    if dist_argument:
        # include euclidean distance score only if euclidean strictness is not None
        if euclidean_strictness is not None:
            # set euclidean linker minimum and maximum with euclidean strictness
            euclidean_linker_minimum = linker_minimum - euclidean_strictness
            euclidean_linker_minimum = euclidean_linker_minimum if euclidean_linker_minimum > 0 else 0
            euclidean_linker_maximum = linker_maximum - euclidean_strictness
            euclidean_linker_maximum = euclidean_linker_maximum if euclidean_linker_maximum > 0 else 0

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

    return round_self(score, 3)


def write_evidence(datapoint, plddt_cutoff, linker_minimum, linker_maximum, euclidean_strictness):
    # combine distance and homo signal reevaluation, create evidence string
    #
    # input datapoint: pd.Series, plddt_cutoff: float, linker_minimum: float, linker_maximum: float,
    # euclidean_strictness: float
    # return evidence: str

    evidence = '\''

    # calculate distance argument of inter protein crosslink confidence score, if plddt included, which is True if
    # entry is not retrieved from AlphaFold or if both plddts surpass or are equal to the given cutoff
    dist_arg = (not pd.isna(datapoint.eucl_dist_tplk)) and \
               (not pd.isna(datapoint.topo_dist_tplk)) and \
               (datapoint.pLDDT_a == '-' or float(datapoint.pLDDT_a) >= plddt_cutoff) and \
               (datapoint.pLDDT_b == '-' or float(datapoint.pLDDT_b) >= plddt_cutoff)

    if dist_arg:
        if euclidean_strictness is not None:
            # set euclidean linker minimum and maximum with euclidean strictness
            euclidean_linker_minimum = linker_minimum - euclidean_strictness
            euclidean_linker_minimum = euclidean_linker_minimum if euclidean_linker_minimum > 0 else 0
            euclidean_linker_maximum = linker_maximum - euclidean_strictness
            euclidean_linker_maximum = euclidean_linker_maximum if euclidean_linker_maximum > 0 else 0

        # set boolean arguments whether euclidean distance is not in linker range, whether topological distance is
        # not in linker range, both distances are lower than minimum, and both distances are higher than maximum
        e_dist_arg = (datapoint.eucl_dist_tplk <= euclidean_linker_minimum) or \
                     (datapoint.eucl_dist_tplk >= euclidean_linker_maximum) \
            if euclidean_strictness is not None else False
        t_dist_arg = (datapoint.topo_dist_tplk <= linker_minimum) or \
                     (datapoint.topo_dist_tplk >= linker_maximum)
        min_dist_arg = (datapoint.eucl_dist_tplk <= euclidean_linker_minimum) and \
                       (datapoint.topo_dist_tplk <= linker_minimum) \
            if euclidean_strictness is not None else datapoint.topo_dist_tplk <= linker_minimum
        max_dist_arg = (datapoint.eucl_dist_tplk >= euclidean_linker_maximum) and \
                       (datapoint.topo_dist_tplk >= linker_maximum) \
            if euclidean_strictness is not None else datapoint.topo_dist_tplk >= linker_maximum

    # write evidence for same peptide if residue positions are equal
    if datapoint.pos_a == datapoint.pos_b:
        evidence += "same peptide"
    # else write split ops and dist evidence
    else:
        ops_evidence = ""
        if datapoint.homo_pep_overl:
            ops_evidence += "peptides overlap" if datapoint.homo_pep_overl else ""

        dist_evidence = ""
        if dist_arg:
            if euclidean_strictness is not None:
                if e_dist_arg and t_dist_arg:
                    if min_dist_arg:
                        dist_evidence += "distances below range"
                    elif max_dist_arg:
                        dist_evidence += "distances above range"
                elif e_dist_arg:
                    if datapoint.eucl_dist_tplk <= euclidean_linker_minimum:
                        dist_evidence += "only euclidean distance below range"
                    elif datapoint.eucl_dist_tplk >= euclidean_linker_maximum:
                        dist_evidence += "only euclidean distance above range"
                elif t_dist_arg:
                    if datapoint.topo_dist_tplk <= linker_minimum:
                        dist_evidence += "only topological distance below range"
                    elif datapoint.topo_dist_tplk >= linker_maximum:
                        dist_evidence += "only topological distance above range"
            else:
                if t_dist_arg:
                    if min_dist_arg:
                        dist_evidence += "distance below range"
                    elif max_dist_arg:
                        dist_evidence += "distance above range"
        if ops_evidence and dist_evidence:
            evidence += ops_evidence + ';' + dist_evidence
        elif ops_evidence:
            evidence += ops_evidence
        elif dist_evidence:
            evidence += dist_evidence
    evidence += '\''
    return evidence if evidence != "\'\'" else ""


def round_self(value, decimals):
    # simple decimal rounding function (python by itself has a tendency to round fragmented with the buit-in function)
    #
    # input value: float, decimals: int
    # return rounded_value: float/int

    # If decimal less than 1, the resulting value will be an integer
    if pd.isna(value):
        return float("Nan")
    if decimals < 1:
        rounded_value = int(int((value * (10 ** decimals)) + .5) / (10 ** decimals))
        return rounded_value
    # Else, the resulting value will be a float
    else:
        rounded_value = int((value * (10 ** decimals)) + .5) / (10 ** decimals)
        return rounded_value
