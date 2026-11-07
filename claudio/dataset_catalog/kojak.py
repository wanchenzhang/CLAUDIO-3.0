"""Read Kojak native results without assigning confidence or protein identities."""

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
    """Import the validated Kojak 1.5.3 layout without selecting source rows."""
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
    if native_version != "1.5.3":
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
