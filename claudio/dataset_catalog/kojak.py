"""Read Kojak native results without assigning confidence or protein identities."""

import json

import csv
import hashlib
import io
import re
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path

from .xifdr import FDR_COLUMNS, _fdr_value, read_search_fasta


SOURCE_COLUMNS = [
    "Scan Number", "Ret Time", "Obs Mass", "Charge",
    "PSM Mass", "PPM Error", "Score", "dScore", "Pep. Diff.",
    "Peptide #1", "Link #1", "Protein #1",
    "Peptide #2", "Link #2", "Protein #2", "Linker Mass",
]
NUMBER = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _integer(raw):
    if re.fullmatch(r"-?\d+", raw or "") is None:
        return None
    return int(raw)


def _peptide_tokens(raw):
    """Separate numeric mass annotations while retaining their original text."""
    if raw == "-":
        return "", []
    sequence, modifications = [], []
    offset = 0
    while offset < len(raw):
        character = raw[offset]
        if "A" <= character <= "Z":
            sequence.append(character)
            offset += 1
            continue
        if character != "[":
            raise ValueError("Unsupported peptide annotation")
        closing = raw.find("]", offset + 1)
        if closing < 0:
            raise ValueError("Unclosed modification token")
        value = raw[offset + 1:closing]
        if NUMBER.fullmatch(value) is None:
            raise ValueError("Unsupported modification token")
        try:
            mass = Decimal(value)
        except InvalidOperation as error:
            raise ValueError("Invalid modification mass") from error
        if not mass.is_finite():
            raise ValueError("Nonfinite modification mass")
        modifications.append({
            "raw_token": raw[offset:closing + 1],
            "position_after_residue": len(sequence),
            "mass_delta_raw": value,
            "chemical_identity_assigned": False,
        })
        offset = closing + 1
    if not sequence:
        raise ValueError("Missing peptide sequence")
    return "".join(sequence), modifications


def _protein_tokens(raw):
    """Retain every reported protein association and its native coordinate."""
    if raw == "-":
        return []
    parts = raw.split(";")
    if parts[-1] == "":
        parts.pop()
    if not parts or any(not part for part in parts):
        raise ValueError("Malformed protein annotation")
    associations = []
    for part in parts:
        match = re.fullmatch(r"(.+)\((-?\d+)\)", part)
        if match:
            identifier, position = match.groups()
            associations.append({
                "raw": part,
                "protein_id": identifier,
                "protein_position": int(position),
            })
        else:
            associations.append({
                "raw": part,
                "protein_id": part,
                "protein_position": None,
            })
    return associations


def _endpoint(fields, number, proteins, decoy_prefix):
    result = {
        "raw_peptide": fields["Peptide #" + number],
        "peptide": None,
        "link_position": _integer(fields["Link #" + number]),
        "modifications": [],
        "protein_candidates": [],
        "status": "not_evaluable",
        "reason": "",
        "kind": "residue",
        "chemical_atom_type": "unspecified",
    }
    try:
        peptide, modifications = _peptide_tokens(result["raw_peptide"])
        associations = _protein_tokens(fields["Protein #" + number])
    except ValueError as error:
        result["reason"] = str(error)
        return result

    result.update(peptide=peptide, modifications=modifications)
    link_position = result["link_position"]
    candidates = []
    for association in associations:
        candidate = dict(association)
        identifier = association["protein_id"]
        position = association["protein_position"]
        candidate.update(
            decoy_label_present=identifier.startswith(decoy_prefix),
            status="not_evaluable",
            reason="",
        )
        if identifier not in proteins:
            candidate["reason"] = "protein_not_in_search_fasta"
        elif position is None or position <= 0:
            candidate["reason"] = "protein_position_unresolved"
        elif (
            link_position is None
            or link_position <= 0
            or link_position > len(peptide)
        ):
            candidate["reason"] = "peptide_link_position_unresolved"
        else:
            sequence = proteins[identifier]["sequence"]
            start = position - link_position + 1
            end = start + len(peptide) - 1
            candidate.update(peptide_start=start, peptide_end=end)
            if start < 1 or end > len(sequence):
                candidate["reason"] = "peptide_outside_search_sequence"
            elif sequence[start - 1:end] != peptide:
                candidate["reason"] = "peptide_sequence_or_coordinate_mismatch"
            elif sequence[position - 1] != peptide[link_position - 1]:
                candidate["reason"] = "link_residue_mismatch"
            else:
                candidate.update(
                    status="valid",
                    residue=sequence[position - 1],
                )
        candidates.append(candidate)

    result["protein_candidates"] = candidates
    distinct = {
        (item["protein_id"], item["protein_position"])
        for item in candidates if item["status"] == "valid"
    }
    all_valid = bool(candidates) and all(
        item["status"] == "valid" for item in candidates
    )
    if all_valid and len(distinct) == 1:
        result.update(status="valid", reason="")
    elif all_valid and len(distinct) > 1:
        result.update(
            status="ambiguous",
            reason="multiple_source_protein_associations",
        )
    else:
        result["reason"] = "source_protein_association_unresolved"
    return result


def _decoy_status(fields, proteins, decoy_prefix):
    associations = []
    for number in ("1", "2"):
        associations.extend(_protein_tokens(fields["Protein #" + number]))
    if not associations:
        return "unknown"
    labels = [
        item["protein_id"].startswith(decoy_prefix)
        for item in associations
    ]
    if all(labels):
        return "decoy_only"
    if any(labels):
        return "contains_target_and_decoy_labels"
    if all(item["protein_id"] in proteins for item in associations):
        return "target_labels_only"
    return "unknown"


def _record_kind(fields):
    first_present = fields["Peptide #1"] != "-"
    second_present = fields["Peptide #2"] != "-"
    first_link = _integer(fields["Link #1"])
    second_link = _integer(fields["Link #2"])
    both_positive = (
        first_link is not None and first_link > 0
        and second_link is not None and second_link > 0
    )
    if first_present and second_present and both_positive:
        return "paired_crosslink_candidate"
    if first_present and not second_present and both_positive:
        return "loop_link_candidate"
    if first_present and not second_present:
        return "single_peptide_record"
    if not first_present and not second_present:
        return "no_peptide_reported"
    return "source_unresolved"


def read_kojak_csms(
    results_path, fasta_path, *, dataset_accession, decoy_prefix
):
    """Import validated Kojak 1.5.1 and 1.5.3 rows without selection."""
    if (
        not isinstance(dataset_accession, str)
        or not dataset_accession.strip()
    ):
        raise ValueError("A nonempty dataset accession is required")
    if (
        not isinstance(decoy_prefix, str)
        or not decoy_prefix.strip()
        or decoy_prefix != decoy_prefix.strip()
    ):
        raise ValueError("An explicit nonempty decoy prefix is required")

    results_path = Path(results_path)
    fasta_path = Path(fasta_path)
    raw = results_path.read_bytes()
    psm_sha256 = hashlib.sha256(raw).hexdigest()
    fasta_sha256 = _sha256(fasta_path)
    proteins = read_search_fasta(fasta_path)

    reader = csv.reader(
        io.StringIO(raw.decode("utf-8-sig")), delimiter="\t"
    )
    rows = [
        (reader.line_num, values)
        for values in reader
        if any(value.strip() for value in values)
    ]
    if len(rows) < 2:
        raise ValueError("Kojak result lacks its version or header")
    version_values = rows[0][1]
    if len(version_values) != 1:
        raise ValueError("Missing Kojak native version marker")
    version = re.fullmatch(
        r"Kojak version ([0-9]+(?:\.[0-9]+)+)", version_values[0]
    )
    if version is None:
        raise ValueError("Missing Kojak native version marker")
    native_version = version.group(1)
    if native_version not in {"1.5.1", "1.5.3"}:
        raise ValueError(
            "Kojak native version is not yet validated: " + native_version
        )
    columns = rows[1][1]
    if columns != SOURCE_COLUMNS:
        raise ValueError("Unsupported Kojak native column layout")

    records = []
    for number, (line_number, values) in enumerate(rows[2:], start=1):
        record = {
            "source_record_id": (
                dataset_accession + ":" + psm_sha256 + ":" + str(number)
            ),
            "source_record_number": number,
            "source_line_number": line_number,
            "source_values": list(values),
            "source_fields": dict(zip(columns, values)),
            "record_kind": "source_unresolved",
            "record_status": "not_evaluable",
            "record_reason": "",
            "source_errors": [],
            "endpoints": {},
            "decoy_status": "unknown",
            "fdr_status": "not_provided",
            "fdr": {
                level: _fdr_value(None) for level in FDR_COLUMNS
            },
            "spectrum_file": None,
            "scan_number": None,
            "source_result_file_name": results_path.name,
            "spectrum_identity_scope": "source_result_file_and_scan",
        }
        records.append(record)
        if len(values) != len(columns):
            record["record_reason"] = "source_field_count_mismatch"
            record["source_errors"].append(record["record_reason"])
            continue

        fields = record["source_fields"]
        record["record_kind"] = _record_kind(fields)
        record["scan_number"] = _integer(fields["Scan Number"])
        record["linker_mass_raw"] = fields["Linker Mass"]

        try:
            record["decoy_status"] = _decoy_status(
                fields, proteins, decoy_prefix
            )
        except ValueError as error:
            record["source_errors"].append(str(error))

        if record["record_kind"] == "paired_crosslink_candidate":
            record["endpoints"] = {
                "a": _endpoint(fields, "1", proteins, decoy_prefix),
                "b": _endpoint(fields, "2", proteins, decoy_prefix),
            }
            if not all(
                endpoint["status"] == "valid"
                for endpoint in record["endpoints"].values()
            ):
                record["source_errors"].append(
                    "source_endpoint_unresolved_or_ambiguous"
                )
            if (
                record["scan_number"] is None
                or record["scan_number"] <= 0
            ):
                record["source_errors"].append(
                    "source_scan_number_unresolved"
                )
            if record["source_errors"]:
                record["record_reason"] = record["source_errors"][0]
            else:
                record["record_status"] = "valid"
        elif record["record_kind"] == "source_unresolved":
            record["record_reason"] = "source_record_kind_unresolved"
            record["source_errors"].append(record["record_reason"])
        else:
            record["record_status"] = "retained_not_pairwise"
            record["record_reason"] = {
                "single_peptide_record": "single_peptide_not_pairwise",
                "loop_link_candidate": "loop_link_requires_single_peptide_support",
                "no_peptide_reported": "source_no_peptide_reported",
            }[record["record_kind"]]

    if _sha256(results_path) != psm_sha256:
        raise ValueError("Kojak result source changed during import")
    if _sha256(fasta_path) != fasta_sha256:
        raise ValueError("Search FASTA source changed during import")

    paired = [
        record for record in records
        if record["record_kind"] == "paired_crosslink_candidate"
    ]
    valid = [
        record for record in paired
        if record["record_status"] == "valid"
    ]
    return {
        "schema_version": 1,
        "format": "kojak_csm",
        "source_columns": list(columns),
        "source_preamble": list(version_values),
        "native_version": native_version,
        "provenance": {
            "dataset_accession": dataset_accession,
            "psm_path": str(results_path.resolve()),
            "psm_sha256": psm_sha256,
            "fasta_path": str(fasta_path.resolve()),
            "fasta_sha256": fasta_sha256,
            "source_file_name": results_path.name,
        },
        "proteins": proteins,
        "records": records,
        "parser": {
            "backend": "native",
            "library": None,
            "library_version": None,
            "native_format_version": native_version,
            "source_records": len(records),
            "paired_records": len(paired),
            "validated_paired_records": len(valid),
            "modification_masses_assessed": False,
            "endpoint_order": "source_order",
        },
        "decoy_policy": {
            "kind": "caller_declared_prefix",
            "prefix": decoy_prefix,
            "target_decoy_validation_performed": False,
        },
        "fdr_policy": {
            "kind": "kojak_native_score_table",
            "verification_status": "unverified",
            "per_record_fdr_applied": False,
            "native_scores_used_as_fdr": False,
        },
        "record_kind_counts": dict(Counter(
            record["record_kind"] for record in records
        )),
        "identity_mapping_applied": False,
        "input_exported": False,
        "assessment_performed": False,
    }


from copy import deepcopy
from .xifdr import _reference_evidence, write_xifdr_outputs


def _verify_kojak_sources(bundle):
    for kind in ("psm", "fasta"):
        path = Path(bundle["provenance"][kind + "_path"])
        if _sha256(path) != bundle["provenance"][kind + "_sha256"]:
            raise ValueError("Kojak input source checksum changed")


def _verify_kojak_receipt(receipt):
    if not isinstance(receipt, dict) or not receipt.get("absolute_path"):
        raise ValueError("Evidence source path is unresolved")
    path = Path(receipt["absolute_path"])
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != receipt.get("sha256"):
        raise ValueError("Kojak reference evidence checksum changed")
    return raw


_KOJAK_INTACT_PEPTIDE_RULE = "intact-peptide-projection-v1"


def _kojak_reviewed_projection_context(
    bundle, identifier, profile, max_optimal_alignments
):
    """Verify reviewed reference coordinates without accepting identity."""
    from .xifdr import _projection_coordinate_options

    source = bundle["proteins"][identifier]
    source_sequence = source["sequence"]
    reference_sequence = profile["reference_sequence"]

    if not isinstance(profile.get("profile_id"), str) or not (
        profile["profile_id"].strip()
    ):
        raise ValueError("Reference profile identifier is missing")

    options, evidence = _projection_coordinate_options(
        source_sequence, reference_sequence, max_optimal_alignments
    )
    if not options or evidence.get("reason"):
        raise ValueError(
            "Reference projection coordinates are unresolved: "
            + str(evidence.get("reason", "no_mapping_options"))
        )
    if profile.get("coordinate_evidence") != evidence:
        raise ValueError("Reference projection coordinates changed")

    mismatches = []
    variable = []
    for position, source_residue in enumerate(source_sequence, start=1):
        positions = [option.get(position) for option in options]
        if any(value is None for value in positions) or len(set(positions)) != 1:
            variable.append(position)
            continue
        reference_position = positions[0]
        reference_residue = reference_sequence[reference_position - 1]
        if source_residue != reference_residue:
            mismatches.append({
                "source_position": position,
                "reference_position": reference_position,
                "source_residue": source_residue,
                "reference_residue": reference_residue,
            })

    review = profile.get("sequence_projection_review")
    if not isinstance(review, dict):
        raise ValueError("A construct sequence projection review is required")
    declaration = {
        key: value for key, value in review.items() if key != "sources"
    }
    required = {
        "schema_version": 1,
        "kind": "sequence_projection_review",
        "review_status": "reviewed",
        "projection_rule_version": _KOJAK_INTACT_PEPTIDE_RULE,
        "source_fasta_id": identifier,
        "source_fasta_sha256": bundle["provenance"]["fasta_sha256"],
        "source_sequence_sha256": hashlib.sha256(
            source_sequence.encode("ascii")
        ).hexdigest(),
        "declared_accession": source.get("uniprot_accession"),
        "reference_accession": profile["accession"],
        "reference_sequence_sha256": hashlib.sha256(
            reference_sequence.encode("ascii")
        ).hexdigest(),
        "taxon_id": profile["taxon_id"],
        "mapping_method": evidence["method"],
        "max_optimal_alignments": max_optimal_alignments,
        "coordinate_evidence": evidence,
        "invariant_mismatches": mismatches,
        "variable_or_unmapped_source_positions": variable,
        "require_all_optimal_alignments": True,
        "intact_peptide_match_required": True,
        "identity_accepted": False,
    }
    if any(declaration.get(key) != value for key, value in required.items()):
        raise ValueError("Construct review differs from recomputed projection")
    if not isinstance(declaration.get("reason"), str) or not (
        declaration["reason"].strip()
    ):
        raise ValueError("Construct review requires a reason")

    receipts = review.get("sources")
    if not isinstance(receipts, list) or not receipts:
        raise ValueError("Construct review evidence is missing")
    for receipt in receipts:
        if receipt.get("kind") != "sequence_projection_review":
            raise ValueError("Unsupported construct review source")
        raw = _verify_kojak_receipt(receipt)
        if json.loads(raw) != declaration:
            raise ValueError("Construct review document changed")

    return {
        "profile": profile,
        "coordinate_options": options,
        "projection_rule_version": _KOJAK_INTACT_PEPTIDE_RULE,
    }


def _kojak_project_reviewed_endpoint(endpoint, candidate, context):
    """Require the same intact peptide projection under every mapping."""
    from .xifdr import _project_peptide

    profile = context["profile"]
    options = context["coordinate_options"]
    source_endpoint = {
        "status": endpoint["status"],
        "peptide": endpoint["peptide"],
        "link_position": endpoint["link_position"],
        "peptide_start": candidate["peptide_start"],
        "protein_position": candidate["protein_position"],
    }
    result = _project_peptide(source_endpoint, profile, options)
    result["source_fasta_id"] = candidate["protein_id"]

    supported = (
        bool(options)
        and result["status"] == "projected"
        and result["compatible_alignment_count"] == len(options)
        and result["distinct_projection_count"] == 1
    )
    result["all_mapping_options_support_projection"] = supported

    if not supported:
        result["status"] = "not_evaluable"
        if not result.get("reason"):
            result["reason"] = "not_supported_by_all_optimal_alignments"
        return result

    position = result["reference_protein_position"]
    residue = profile["reference_sequence"][position - 1]
    if residue != candidate["residue"]:
        result.update(
            status="not_evaluable",
            reason="reference_residue_mismatch",
        )
        return result

    result["residue"] = residue
    return result

def _kojak_reference_profiles(bundle, mapping):
    if (
        not isinstance(mapping, dict)
        or mapping.get("schema_version") != 1
        or mapping.get("kind") != "kojak_reference_mapping"
    ):
        raise ValueError("Unsupported Kojak reference mapping")
    projection_rule = mapping.get("projection_rule_version")
    if projection_rule not in (None, _KOJAK_INTACT_PEPTIDE_RULE):
        raise ValueError("Unsupported Kojak reference projection rule")
    max_optimal_alignments = mapping.get("max_optimal_alignments")
    if projection_rule == _KOJAK_INTACT_PEPTIDE_RULE:
        if (
            type(max_optimal_alignments) is not int
            or not 1 <= max_optimal_alignments <= 64
        ):
            raise ValueError("Choose a reviewed alignment limit from 1 to 64")

    for kind in ("psm", "fasta"):
        if mapping.get("source_" + kind + "_sha256") != (
            bundle["provenance"][kind + "_sha256"]
        ):
            raise ValueError("Reference mapping source checksum differs")
    if (
        mapping.get("configuration_application_verified") is not False
        or mapping.get("crosslinker_identity_verified") is not False
    ):
        raise ValueError("Unsupported Kojak confidence or chemistry declaration")
    for receipt in mapping.get("context_sources", []):
        _verify_kojak_receipt(receipt)

    references = mapping.get("references")
    if not isinstance(references, dict) or not references:
        raise ValueError("Reference mapping has no profiles")

    checked = {}
    identities = {}
    for identifier, profile in references.items():
        if identifier not in bundle["proteins"]:
            raise ValueError("Reference source protein is missing")
        if identifier.startswith(bundle["decoy_policy"]["prefix"]):
            raise ValueError("Decoy proteins cannot receive target references")

        source = bundle["proteins"][identifier]
        source_sequence = source["sequence"]
        reference_sequence = profile.get("reference_sequence")
        if (
            profile.get("source_fasta_id") != identifier
            or profile.get("source_sequence") != source_sequence
        ):
            raise ValueError("Reference source sequence changed")
        if (
            not isinstance(reference_sequence, str)
            or re.fullmatch(r"[A-Z]+", reference_sequence) is None
        ):
            raise ValueError("Invalid reference sequence")
        reference_hash = hashlib.sha256(
            reference_sequence.encode("ascii")
        ).hexdigest()
        if profile.get("reference_sequence_sha256") != reference_hash:
            raise ValueError("Reference sequence checksum changed")
        if (
            profile.get("identity_scope") != "candidate_hypothesis"
            or profile.get("identity_accepted") is not False
            or profile.get("require_all_optimal_alignments") is not True
        ):
            raise ValueError("Unsupported reference identity declaration")

        receipts = profile.get("sources")
        if not isinstance(receipts, list) or not receipts:
            raise ValueError("Reference evidence is missing")
        for receipt in receipts:
            _verify_kojak_receipt(receipt)
        _reference_evidence(
            profile, profile, source,
            require_declared_accession=False,
        )

        if projection_rule == _KOJAK_INTACT_PEPTIDE_RULE:
            checked[identifier] = _kojak_reviewed_projection_context(
                bundle, identifier, profile, max_optimal_alignments
            )
        else:
            starts = [
                position for position in range(len(source_sequence))
                if source_sequence.startswith(reference_sequence, position)
            ]
            if len(starts) != 1:
                raise ValueError("Reference containment is unresolved or ambiguous")
            offset = starts[0]
            end = offset + len(reference_sequence)
            coordinate_evidence = {
                "method": "full_reference_exact_subsequence",
                "source_start": offset + 1,
                "source_end": end,
                "reference_start": 1,
                "reference_end": len(reference_sequence),
                "mapping_options": 1,
            }
            if profile.get("coordinate_evidence") != coordinate_evidence:
                raise ValueError("Reference projection coordinates changed")

            review = profile.get("sequence_projection_review")
            if not isinstance(review, dict):
                raise ValueError("A construct sequence projection review is required")
            declaration = {
                key: value for key, value in review.items()
                if key != "sources"
            }
            unmapped = [
                position for position in range(1, len(source_sequence) + 1)
                if not offset + 1 <= position <= end
            ]
            required = {
                "schema_version": 1,
                "kind": "sequence_projection_review",
                "review_status": "reviewed",
                "source_fasta_id": identifier,
                "source_fasta_sha256": bundle["provenance"]["fasta_sha256"],
                "source_sequence_sha256": hashlib.sha256(
                    source_sequence.encode("ascii")
                ).hexdigest(),
                "declared_accession": source.get("uniprot_accession"),
                "reference_accession": profile["accession"],
                "reference_sequence_sha256": reference_hash,
                "taxon_id": profile["taxon_id"],
                "mapping_method": "full_reference_exact_subsequence",
                "source_interval": [offset + 1, end],
                "reference_interval": [1, len(reference_sequence)],
                "allowed_source_unmapped_positions": unmapped,
                "require_all_optimal_alignments": True,
                "identity_accepted": False,
            }
            if any(declaration.get(key) != value for key, value in required.items()):
                raise ValueError("Construct review differs from recomputed projection")
            if not isinstance(declaration.get("reason"), str) or not (
                declaration["reason"].strip()
            ):
                raise ValueError("Construct review requires a reason")
            review_sources = review.get("sources")
            if not isinstance(review_sources, list) or not review_sources:
                raise ValueError("Construct review evidence is missing")
            for receipt in review_sources:
                if receipt.get("kind") != "sequence_projection_review":
                    raise ValueError("Unsupported construct review source")
                raw = _verify_kojak_receipt(receipt)
                if json.loads(raw) != declaration:
                    raise ValueError("Construct review document changed")

            checked[identifier] = {
                "profile": profile,
                "source_start": offset + 1,
                "source_end": end,
                "offset": offset,
            }
        identities[identifier] = {
            "status": "candidate_hypothesis",
            "accession": profile["accession"],
            "sequence": reference_sequence,
            "taxon_id": profile["taxon_id"],
            "identity_scope": "candidate_hypothesis",
            "identity_accepted": False,
            "sources": deepcopy(profile["sources"]),
        }
    return checked, identities


def write_kojak_claudio_input(
    bundle, reference_mapping, output_dir, *,
    allowed_linker_masses, allowed_residues,
    export_unverified_candidates=False,
):
    """Export explicitly selected native candidates with all source rows retained."""
    if export_unverified_candidates is not True:
        raise ValueError("Explicit unverified candidate export mode is required")
    if Path(output_dir).exists():
        raise FileExistsError("Output directory already exists")
    if bundle.get("format") != "kojak_csm":
        raise ValueError("Expected a native Kojak import")

    selections = {}
    for name, values in (
        ("linker_masses", allowed_linker_masses),
        ("residues", allowed_residues),
    ):
        if (
            not isinstance(values, (set, frozenset, list, tuple))
            or not values
            or any(
                not isinstance(value, str) or not value.strip()
                for value in values
            )
        ):
            raise ValueError("Explicit nonempty " + name + " are required")
        selections[name] = set(values)
    if any(
        re.fullmatch(r"[A-Z]", residue) is None
        for residue in selections["residues"]
    ):
        raise ValueError("Residues must be single uppercase letters")
    masses = set()
    for raw in selections["linker_masses"]:
        if NUMBER.fullmatch(raw) is None:
            raise ValueError("Invalid selected linker mass")
        value = Decimal(raw)
        if not value.is_finite():
            raise ValueError("Selected linker mass must be finite")
        masses.add(value)

    _verify_kojak_sources(bundle)
    provenance = bundle["provenance"]
    native = read_kojak_csms(
        provenance["psm_path"], provenance["fasta_path"],
        dataset_accession=provenance["dataset_accession"],
        decoy_prefix=bundle["decoy_policy"]["prefix"],
    )
    for key in (
        "provenance", "source_columns", "source_preamble", "native_version",
        "proteins", "records", "parser", "decoy_policy", "fdr_policy",
    ):
        if bundle.get(key) != native.get(key):
            raise ValueError("Native Kojak source evidence changed after import")

    profiles, identities = _kojak_reference_profiles(
        bundle, reference_mapping
    )
    rows, audit = [], []
    for record in bundle["records"]:
        entry = deepcopy(record)
        reasons = []
        paired = record["record_kind"] == "paired_crosslink_candidate"
        if not paired:
            reasons.append(
                record["record_reason"] or "not_pairwise_crosslink"
            )
        if record["decoy_status"] in (
            "decoy_only", "contains_target_and_decoy_labels",
        ):
            reasons.append("decoy_record")
        elif paired and record["decoy_status"] != "target_labels_only":
            reasons.append("source_decoy_status_unresolved")
        if paired and record["record_status"] != "valid":
            reasons.append(
                record["record_reason"] or "source_endpoint_unresolved"
            )

        projected = {}
        if not reasons:
            raw_mass = record["linker_mass_raw"]
            if NUMBER.fullmatch(raw_mass) is None:
                reasons.append("source_linker_mass_unresolved")
            elif Decimal(raw_mass) not in masses:
                reasons.append("linker_mass_not_selected")

        if not reasons:
            for side in ("a", "b"):
                endpoint = record["endpoints"][side]
                candidates = endpoint["protein_candidates"]
                if len(candidates) != 1:
                    reasons.append("source_protein_association_unresolved_" + side)
                    continue
                candidate = candidates[0]
                identifier = candidate["protein_id"]
                context = profiles.get(identifier)
                if context is None:
                    reasons.append("protein_identity_unresolved_" + side)
                    continue
                if (
                    context.get("projection_rule_version")
                    == _KOJAK_INTACT_PEPTIDE_RULE
                ):
                    projection = _kojak_project_reviewed_endpoint(
                        endpoint, candidate, context
                    )
                    projected[side] = projection
                    if projection["status"] != "projected":
                        reasons.append("reference_projection_unresolved_" + side)
                        continue
                    if projection["residue"] not in selections["residues"]:
                        reasons.append("residue_not_selected_" + side)
                    continue

                profile = context["profile"]
                start = candidate["peptide_start"]
                end = candidate["peptide_end"]
                if (
                    start < context["source_start"]
                    or end > context["source_end"]
                ):
                    reasons.append("peptide_overlaps_unmapped_source_" + side)
                    continue
                reference_start = start - context["offset"]
                reference_end = end - context["offset"]
                position = (
                    candidate["protein_position"] - context["offset"]
                )
                sequence = profile["reference_sequence"]
                if sequence[reference_start - 1:reference_end] != (
                    endpoint["peptide"]
                ):
                    reasons.append("reference_peptide_mismatch_" + side)
                    continue
                residue = sequence[position - 1]
                if residue != candidate["residue"]:
                    reasons.append("reference_residue_mismatch_" + side)
                    continue
                projected[side] = {
                    "source_fasta_id": identifier,
                    "accession": profile["accession"],
                    "identity_scope": "candidate_hypothesis",
                    "source_protein_position": candidate["protein_position"],
                    "reference_protein_position": position,
                    "reference_peptide_start": reference_start,
                    "reference_peptide_end": reference_end,
                    "residue": residue,
                    "status": "projected",
                }
                if residue not in selections["residues"]:
                    reasons.append("residue_not_selected_" + side)

        entry.update(
            input_status="not_evaluable" if reasons else "candidate_ready",
            input_reasons=list(dict.fromkeys(reasons)),
            reference_endpoints=deepcopy(projected),
        )
        audit.append(entry)
        if reasons:
            continue

        fields = record["source_fields"]
        row = {
            "source_record_id": record["source_record_id"],
            "source_record_number": record["source_record_number"],
            "source_line_number": record["source_line_number"],
            "dataset_accession": provenance["dataset_accession"],
            "source_psm_sha256": provenance["psm_sha256"],
            "source_fasta_sha256": provenance["fasta_sha256"],
            "source_result_file_name": record["source_result_file_name"],
            "source_scan_number": record["scan_number"],
            "source_spectrum_file": record["spectrum_file"],
            "source_linker_mass": record["linker_mass_raw"],
            "source_decoy_status": record["decoy_status"],
            "source_score": fields["Score"],
            "source_delta_score": fields["dScore"],
            "source_peptide_difference": fields["Pep. Diff."],
            "source_retention_time_raw": fields["Ret Time"],
            "source_fields_json": json.dumps(
                fields, ensure_ascii=False, sort_keys=True
            ),
            "input_scope": "unverified_candidate",
            "identity_scope": "candidate_hypothesis",
            "protein_identity_accepted": False,
            "crosslinker_identity_verified": False,
            "configuration_application_verified": False,
            "fdr_verification_status": "unverified",
            "selected_fdr_level": None,
            "selected_fdr_value": None,
            "selected_fdr_cutoff": None,
        }
        for side in ("a", "b"):
            endpoint = record["endpoints"][side]
            candidate = endpoint["protein_candidates"][0]
            projection = projected[side]
            profile = profiles[candidate["protein_id"]]["profile"]
            row.update({
                "pep_" + side: endpoint["peptide"],
                "pos_" + side: projection["reference_protein_position"],
                "res_pos_" + side: endpoint["link_position"],
                "unip_id_" + side: profile["accession"],
                "seq_" + side: profile["reference_sequence"],
                "source_protein_id_" + side: candidate["protein_id"],
                "source_protein_position_" + side: candidate["protein_position"],
                "source_peptide_start_" + side: candidate["peptide_start"],
                "source_peptide_raw_" + side: endpoint["raw_peptide"],
                "source_modifications_json_" + side: json.dumps(
                    endpoint["modifications"], sort_keys=True
                ),
            })
        rows.append(row)

    bodies = {}
    for kind, suffix in (("psm", ".kojak.txt"), ("fasta", ".fasta")):
        path = Path(provenance[kind + "_path"])
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != provenance[kind + "_sha256"]:
            raise ValueError("Kojak input source checksum changed before output")
        bodies[
            "native_sources/" + provenance[kind + "_sha256"] + suffix
        ] = raw

    prepared = {
        "rows": rows,
        "audit": audit,
        "extra_source_bodies": bodies,
        "reference_projection": deepcopy(reference_mapping),
        "reference_projection_policy": {
            "kind": (
                "reviewed_intact_peptide_projection"
                if reference_mapping.get("projection_rule_version")
                == _KOJAK_INTACT_PEPTIDE_RULE
                else "reviewed_full_reference_exact_subsequence"
            ),
            "require_all_optimal_alignments": True,
            "identity_accepted": False,
        },
        "fdr_policy": {
            "kind": "unverified_candidates",
            "verification_status": "unverified",
            "level": None,
            "cutoff": None,
            "per_record_fdr_applied": False,
            "author_filter_applied": False,
            "native_scores_used_as_fdr": False,
            "target_decoy_validation_performed": False,
        },
    }
    return write_xifdr_outputs(
        bundle, identities, prepared, output_dir,
        pipeline_policy={
            "input_scope": "unverified_candidate",
            "linker_masses": sorted(str(value) for value in masses),
            "residues": sorted(selections["residues"]),
            "coordinate_basis": 1,
            "decoy_policy": "exclude_decoy_and_unresolved_labels",
            "configuration_application_verified": False,
            "crosslinker_identity_verified": False,
            "protein_association_policy": "single_reported_association",
        },
    )
