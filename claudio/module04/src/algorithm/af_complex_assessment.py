"""Compatibility interface for AF complex chain annotations."""
from claudio.module04.src.algorithm.structure_assessment import (
    annotate_structure_relationships,
)


def annotate_complex_relationships(data):
    result = annotate_structure_relationships(data, structure_source="af_complex")
    result["chain_relation_rule_version"] = "af-chain-v1"
    return result