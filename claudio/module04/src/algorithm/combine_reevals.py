"""Combines the results from Module02's distance analysis and
Module03's OPS analysis, and reevaluating assigned cross-link type."""
import pandas as pd

from claudio.utils.utils import round_self
from claudio.module04.src.algorithm.relationship_classification import (
    classify_model_candidates,
)

def combine_inter_reevaluations(data: pd.DataFrame, plddt_cutoff: float,
                                linker_minimum: float, linker_maximum: float,
                                euclidean_strictness: float,
                                distance_maximum: float, cutoff: float,
                                compute_scoring: bool):
    """Retain legacy heuristics as provenance and classify model candidates."""
    result = data.copy(deep=True)
    for column in (
        "XL_type", "XL_confirmed", "evidence", "inter_score", "score_XL_type",
    ):
        if column in result:
            if "XL_type_scope" not in result:
                result[f"input_{column}"] = result[column]
            result = result.drop(columns=column)

    columns = [
        "legacy_evidence", "legacy_XL_type", "legacy_XL_confirmed",
        "legacy_inference_status", "legacy_inference_reason",
    ]
    if compute_scoring:
        columns.extend(["legacy_inter_score", "legacy_score_XL_type"])
    else:
        result = result.drop(
            columns=["legacy_inter_score", "legacy_score_XL_type"],
            errors="ignore",
        )

    records = []
    for _, row in result.iterrows():
        record = dict.fromkeys(columns)
        record["legacy_inference_status"] = "computed"
        try:
            evidence = write_evidence(
                row, plddt_cutoff, linker_minimum, linker_maximum,
                euclidean_strictness,
            )
            record["legacy_evidence"] = evidence
            record["legacy_XL_type"] = (
                "intra" if row.unip_id_a == row.unip_id_b
                and row.chain_a == row.chain_b and not evidence else "inter"
            )
            record["legacy_XL_confirmed"] = (
                pd.notna(row.topo_dist) and "distance" not in evidence
            )
            if compute_scoring:
                score = score_inter_potential(
                    row, plddt_cutoff, linker_minimum, linker_maximum,
                    euclidean_strictness, distance_maximum,
                )
                record["legacy_inter_score"] = score
                record["legacy_score_XL_type"] = (
                    "intra" if row.unip_id_a == row.unip_id_b
                    and score <= cutoff else "inter"
                )
        except (
            AttributeError, TypeError, ValueError,
            ZeroDivisionError, OverflowError,
        ) as error:
            record["legacy_inference_status"] = "unavailable"
            record["legacy_inference_reason"] = type(error).__name__
        records.append(record)

    for column in columns:
        result[column] = pd.Series(
            [record[column] for record in records],
            index=result.index,
            dtype="object",
        )

    result["legacy_inference_rule_version"] = "legacy-m4-heuristic-v1"
    return classify_model_candidates(result)


def score_inter_potential(datapoint: pd.Series, plddt_cutoff: float,
                          linker_minimum: float, linker_maximum: float,
                          euclidean_strictness: float, distance_maximum: float):
    """
    combine distance and homo signal reevaluation, create score that 
    represents inter interaction affinity (higher value, higher affinity)
    
    Parameters
    ----------
    input datapoint : pd.Series,
    plddt_cutoff : float,
    linker_minimum : float,
    linker_maximum : float,
    euclidean_strictness : float,
    distance_maximum : float

    Returns
    -------
    score : (float | int)
    """

    score = 0.0

    # calculate distance argument of inter protein crosslink confidence score,
    # if plddt included, which is True if entry is not retrieved from
    # AlphaFold or if both plddts surpass or are equal to the given cutoff
    dist_argument = (not pd.isna(datapoint.eucl_dist)) and \
                    (not pd.isna(datapoint.topo_dist)) and \
                    (datapoint.pLDDT_a == '-' or float(datapoint.pLDDT_a) >= plddt_cutoff) and \
                    (datapoint.pLDDT_b == '-' or float(datapoint.pLDDT_b) >= plddt_cutoff)

    if dist_argument:
        # include euclidean distance score only if euclidean strictness
        # is not None
        if euclidean_strictness is not None:
            # set euclidean linker minimum and maximum with euclidean strictness
            euclidean_linker_minimum = linker_minimum - euclidean_strictness
            euclidean_linker_minimum = (euclidean_linker_minimum
                                        if euclidean_linker_minimum > 0 else 0)
            euclidean_linker_maximum = linker_maximum - euclidean_strictness
            euclidean_linker_maximum = (euclidean_linker_maximum
                                        if euclidean_linker_maximum > 0 else 0)

            # calculate euclidean distance score
            eucl_dist = (datapoint.eucl_dist
                         if datapoint.eucl_dist < distance_maximum else distance_maximum)
            raw_eucl_score = (eucl_dist - euclidean_linker_maximum) / (distance_maximum - euclidean_linker_maximum)
            if datapoint.eucl_dist <= euclidean_linker_minimum:
                score += .25
            elif datapoint.eucl_dist >= euclidean_linker_maximum:
                score += raw_eucl_score * .25

        # calculate topological distance score
        topo_dist = (datapoint.topo_dist
                     if datapoint.topo_dist < distance_maximum else distance_maximum)
        raw_topo_score = (topo_dist - linker_maximum) / (distance_maximum - linker_maximum)
        if datapoint.topo_dist <= linker_minimum:
            # If the euclidean strictness is set to None the max score for
            # topological is 0.5, else 0.25
            score += .25 if euclidean_strictness is not None else .5
        elif datapoint.topo_dist >= linker_maximum:
            # If the euclidean strictness is set to None the max score
            # for topological is 0.5, else 0.25
            score += (raw_topo_score * .25
                      if euclidean_strictness is not None else raw_topo_score * .5)

    # calculate OPS argument of inter protein crosslink confidence score
    overlap = datapoint.homo_pep_overl
    if pd.notna(overlap) and overlap:
        score += datapoint.homo_int_overl / 2

    return round_self(score, 3)


def write_evidence(datapoint: pd.Series, plddt_cutoff: float,
                   linker_minimum: float, linker_maximum: float,
                   euclidean_strictness: float):
    """
    combine distance and homo signal reevaluation, create evidence string

    Parameters
    ----------
    datapoint : pd.Series,
    plddt_cutoff : float,
    linker_minimum : float,
    linker_maximum : float,
    euclidean_strictness : float

    Returns
    -------
    evidence : str
    """

    evidence = '\''

    # calculate distance argument of inter protein crosslink confidence score,
    # if plddt included, which is True if entry is not retrieved from AlphaFold
    # or if both plddts surpass or are equal to the given cutoff
    dist_arg = (not pd.isna(datapoint.eucl_dist)) and \
               (not pd.isna(datapoint.topo_dist)) and \
               (datapoint.pLDDT_a == '-' or float(datapoint.pLDDT_a) >= plddt_cutoff) and \
               (datapoint.pLDDT_b == '-' or float(datapoint.pLDDT_b) >= plddt_cutoff)

    if dist_arg:
        if euclidean_strictness is not None:
            # set euclidean linker minimum and maximum with euclidean strictness
            euclidean_linker_minimum = linker_minimum - euclidean_strictness
            euclidean_linker_minimum = (euclidean_linker_minimum
                                        if euclidean_linker_minimum > 0 else 0)
            euclidean_linker_maximum = linker_maximum - euclidean_strictness
            euclidean_linker_maximum = (euclidean_linker_maximum
                                        if euclidean_linker_maximum > 0 else 0)

        # set boolean arguments whether euclidean distance is not in
        # linker range, whether topological distance is not in linker range,
        # both distances are lower than minimum, and both distances are higher
        # than maximum
        e_dist_arg = (datapoint.eucl_dist <= euclidean_linker_minimum) or \
                     (datapoint.eucl_dist >= euclidean_linker_maximum) \
            if euclidean_strictness is not None else False
        t_dist_arg = (datapoint.topo_dist <= linker_minimum) or \
                     (datapoint.topo_dist >= linker_maximum)
        min_dist_arg = (datapoint.eucl_dist <= euclidean_linker_minimum) and \
                       (datapoint.topo_dist <= linker_minimum) \
            if euclidean_strictness is not None else datapoint.topo_dist <= linker_minimum
        max_dist_arg = (datapoint.eucl_dist >= euclidean_linker_maximum) and \
                       (datapoint.topo_dist >= linker_maximum) \
            if euclidean_strictness is not None else datapoint.topo_dist >= linker_maximum

    # write evidence for same peptide if residue positions are equal
    if (datapoint.pos_a == datapoint.pos_b) and \
            (datapoint.unip_id_a == datapoint.unip_id_b) and \
            (datapoint.chain_a == datapoint.chain_b):
        evidence += "same peptide"
    # else write split ops and dist evidence
    else:
        ops_evidence = ""
        overlap = datapoint.homo_pep_overl
        if (
            pd.notna(overlap)
            and overlap
            and datapoint.unip_id_a == datapoint.unip_id_b
            and datapoint.chain_a == datapoint.chain_b
        ):
            ops_evidence += "peptides overlap"

        dist_evidence = ""
        if dist_arg:
            if euclidean_strictness is not None:
                if e_dist_arg and t_dist_arg:
                    if min_dist_arg:
                        dist_evidence += "distances below range"
                    elif max_dist_arg:
                        dist_evidence += "distances above range"
                elif e_dist_arg:
                    if datapoint.eucl_dist <= euclidean_linker_minimum:
                        dist_evidence += "only euclidean distance below range"
                    elif datapoint.eucl_dist >= euclidean_linker_maximum:
                        dist_evidence += "only euclidean distance above range"
                elif t_dist_arg:
                    if datapoint.topo_dist <= linker_minimum:
                        dist_evidence += "only topological distance below range"
                    elif datapoint.topo_dist >= linker_maximum:
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


def def_xl_type(datapoint: pd.Series):
    """Return the relationship of this resolved model candidate."""
    return classify_model_candidates(
        pd.DataFrame([datapoint])
    ).iloc[0]["XL_type"]