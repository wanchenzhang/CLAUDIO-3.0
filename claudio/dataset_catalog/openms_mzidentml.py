"""Read legacy OpenMS labeled mzIdentML without rewriting source evidence."""

import gzip
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
from pathlib import Path

from .xifdr import FDR_COLUMNS, _fdr_value, read_search_fasta


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tag(element):
    return element.tag.rsplit("}", 1)[-1]


def _children(element, name):
    return [child for child in element if _tag(child) == name]


def _parameters(element):
    return [
        {"parameter_type": _tag(child), **child.attrib}
        for child in element
        if _tag(child) in ("cvParam", "userParam")
    ]


def _item_snapshot(item):
    return {
        "attributes": dict(item.attrib),
        "parameters": _parameters(item),
        "peptide_evidence_refs": [
            dict(ref.attrib)
            for ref in _children(item, "PeptideEvidenceRef")
        ],
    }


def _sequence(element, name, *, optional=False):
    entries = _children(element, name)
    if not entries and optional:
        return None
    if len(entries) != 1:
        raise ValueError("Sequence element is missing or ambiguous")
    return "".join("".join(entries[0].itertext()).split())


def _integer(raw):
    if re.fullmatch(r"-?\d+", raw or "") is None:
        raise ValueError("Invalid native integer")
    return int(raw)


def _read_xml(path, maximum):
    with Path(path).open("rb") as handle:
        signature = handle.read(2)
        handle.seek(0)
        if signature == b"\x1f\x8b":
            with gzip.GzipFile(fileobj=handle) as decoded:
                raw = decoded.read(maximum + 1)
        else:
            raw = handle.read(maximum + 1)
    if len(raw) > maximum:
        raise ValueError("Native XML exceeds the configured size limit")
    return raw


def _index_elements(root, name):
    result = {}
    for element in root.iter():
        if _tag(element) != name:
            continue
        identifier = element.attrib.get("id")
        if not identifier or identifier in result:
            raise ValueError("Missing or duplicate native " + name + " ID")
        result[identifier] = element
    return result


def _evidence_check(identifier, index, proteins, decoy_prefix):
    result = {
        "peptide_evidence_id": identifier,
        "status": "not_evaluable",
        "reason": "",
        "native_is_decoy": None,
    }
    try:
        evidence = index["PeptideEvidence"][identifier]
        result["native_attributes"] = dict(evidence.attrib)
        raw_decoy = evidence.attrib.get("isDecoy")
        decoy = {
            "0": False, "1": True, "false": False, "true": True,
        }.get(raw_decoy)
        result.update(
            native_is_decoy_raw=raw_decoy,
            native_is_decoy=decoy,
        )
        if raw_decoy is not None and decoy is None:
            raise ValueError("invalid_native_decoy_annotation")

        peptide = index["Peptide"][evidence.attrib["peptide_ref"]]
        database = index["DBSequence"][evidence.attrib["dBSequence_ref"]]
        peptide_sequence = _sequence(peptide, "PeptideSequence")
        protein_id = database.attrib["accession"]
        native_start = _integer(evidence.attrib.get("start"))
        native_end = _integer(evidence.attrib.get("end"))
        start, end = native_start + 1, native_end + 1
        result.update(
            source_protein_id=protein_id,
            peptide=peptide_sequence,
            native_start=native_start,
            native_end=native_end,
            peptide_start=start,
            peptide_end=end,
        )
        if protein_id not in proteins:
            raise ValueError("protein_not_in_saved_fasta")

        source_sequence = proteins[protein_id]["sequence"]
        embedded = _sequence(database, "Seq", optional=True)
        result["embedded_sequence_status"] = (
            "not_present" if embedded is None else "present"
        )
        if embedded is not None and embedded != source_sequence:
            raise ValueError("embedded_sequence_differs_from_fasta")
        if (
            start < 1 or end > len(source_sequence)
            or end - start + 1 != len(peptide_sequence)
            or source_sequence[start - 1:end] != peptide_sequence
        ):
            raise ValueError("intact_peptide_or_coordinate_mismatch")
        result["fasta_decoy_label"] = protein_id.startswith(decoy_prefix)
        if decoy is not None and decoy != result["fasta_decoy_label"]:
            raise ValueError("native_decoy_annotation_conflict")
        result["status"] = "valid"
    except (KeyError, ValueError) as error:
        result["reason"] = str(error)
    return result


def _role_sites(peptide_id, peptide):
    sites = []
    for modification in _children(peptide, "Modification"):
        for parameter in _parameters(modification):
            if parameter.get("accession") in ("MS:1002509", "MS:1002510"):
                sites.append({
                    "peptide_ref": peptide_id,
                    "role": parameter["accession"],
                    "link_token": parameter.get("value"),
                    "location_raw": modification.attrib.get("location"),
                    "modification_attributes": dict(modification.attrib),
                    "modification_parameters": _parameters(modification),
                })
    return sites


def _normalize_site(site, items, evidence, proteins):
    relevant = [
        item for item in items
        if item.attrib.get("peptide_ref") == site["peptide_ref"]
    ]
    evidence_ids = list(dict.fromkeys(
        ref.attrib["peptideEvidence_ref"]
        for item in relevant
        for ref in _children(item, "PeptideEvidenceRef")
    ))
    checked = [evidence[value] for value in evidence_ids]
    signatures = {
        (
            value.get("source_protein_id"),
            value.get("peptide_start"),
            value.get("peptide_end"),
            value.get("peptide"),
        )
        for value in checked if value["status"] == "valid"
    }
    if (
        not checked
        or any(value["status"] != "valid" for value in checked)
        or len(signatures) != 1
    ):
        raise ValueError("source_association_unresolved")
    protein_id, start, end, peptide_sequence = next(iter(signatures))
    location = _integer(site["location_raw"])
    if not 1 <= location <= len(peptide_sequence):
        raise ValueError("terminal_or_invalid_site_requires_review")
    position = start + location - 1
    residue = proteins[protein_id]["sequence"][position - 1]
    if residue != peptide_sequence[location - 1]:
        raise ValueError("source_residue_mismatch")
    return {
        **site,
        "source_protein_id": protein_id,
        "peptide": peptide_sequence,
        "peptide_start": start,
        "peptide_end": end,
        "link_position": location,
        "protein_position": position,
        "residue": residue,
        "evidence_ids": evidence_ids,
    }


def _native_fdr(items):
    annotations = [
        {"item_id": item.attrib["id"], "parameter": parameter}
        for item in items
        for parameter in _parameters(item)
        if parameter.get("name") == "OpenXQuest:fdr"
    ]
    raw_values = sorted({
        item["parameter"].get("value", "") for item in annotations
    })
    parsed = [_fdr_value(value) for value in raw_values]
    result = {"raw_values": raw_values, "value": None, "status": "missing"}
    if parsed:
        invalid = [value for value in parsed if value["status"] != "valid"]
        if invalid:
            result["status"] = invalid[0]["status"]
        elif len({value["value"] for value in parsed}) != 1:
            result["status"] = "conflicting"
        else:
            result.update(status="valid", value=parsed[0]["value"])
    return result, annotations


def read_openms_labeled_mzidentml(
    results_path, fasta_path, *, dataset_accession,
    decoy_prefix="decoy_", max_xml_bytes=16 * 1024 * 1024,
):
    """Read the legacy labeled profile using zero-based PeptideEvidence positions."""
    if not isinstance(dataset_accession, str) or not dataset_accession.strip():
        raise ValueError("A nonempty dataset accession is required")
    if not isinstance(decoy_prefix, str) or not decoy_prefix.strip():
        raise ValueError("A nonempty decoy prefix is required")
    if type(max_xml_bytes) is not int or max_xml_bytes <= 0:
        raise ValueError("The XML size limit must be a positive integer")

    results_path, fasta_path = Path(results_path), Path(fasta_path)
    psm_sha256, fasta_sha256 = _sha256(results_path), _sha256(fasta_path)
    raw = _read_xml(results_path, max_xml_bytes)
    root = ET.fromstring(raw)
    if _tag(root) != "MzIdentML":
        raise ValueError("Expected an mzIdentML document")
    proteins = read_search_fasta(fasta_path)
    names = (
        "DBSequence", "Peptide", "PeptideEvidence", "SpectraData",
        "SpectrumIdentificationResult", "SpectrumIdentificationItem",
    )
    index = {name: _index_elements(root, name) for name in names}
    sirs = list(index["SpectrumIdentificationResult"].values())

    groups, unresolved_items = [], []
    source_identifications = []
    for sir in sirs:
        local = defaultdict(list)
        source_items = _children(sir, "SpectrumIdentificationItem")
        source_identifications.append({
            "attributes": dict(sir.attrib),
            "parameters": _parameters(sir),
            "items": [_item_snapshot(item) for item in source_items],
        })
        for item in source_items:
            markers = [
                value for value in _parameters(item)
                if value.get("accession") == "MS:1002511"
            ]
            if len(markers) != 1 or not markers[0].get("value"):
                unresolved_items.append({
                    "sir_attributes": dict(sir.attrib),
                    "item": _item_snapshot(item),
                    "reason": "missing_or_ambiguous_source_group_marker",
                })
                continue
            local[markers[0]["value"]].append(item)
        groups.extend(
            (sir, token, items) for token, items in local.items()
        )

    used_evidence = {
        ref.attrib["peptideEvidence_ref"]
        for _, _, items in groups for item in items
        for ref in _children(item, "PeptideEvidenceRef")
    }
    evidence = {
        identifier: _evidence_check(
            identifier, index, proteins, decoy_prefix
        )
        for identifier in sorted(used_evidence)
    }

    records = []
    for ordinal, (sir, token, items) in enumerate(groups, start=1):
        identifier_text = json.dumps(
            [sir.attrib["id"], token],
            ensure_ascii=False, separators=(",", ":"),
        )
        record_id = (
            dataset_accession + ":" + psm_sha256 + ":"
            + hashlib.sha256(identifier_text.encode("utf-8")).hexdigest()
        )
        peptide_ids = list(dict.fromkeys(
            item.attrib.get("peptide_ref") for item in items
        ))
        errors, roles = [], []
        for peptide_id in peptide_ids:
            peptide = index["Peptide"].get(peptide_id)
            if peptide is None:
                errors.append("source_peptide_reference_unresolved")
            else:
                roles.extend(_role_sites(peptide_id, peptide))
        donors = [
            site for site in roles if site["role"] == "MS:1002509"
        ]
        acceptors = [
            site for site in roles if site["role"] == "MS:1002510"
        ]
        if len(donors) == len(acceptors) == 1:
            kind = (
                "loop_link_group"
                if donors[0]["peptide_ref"] == acceptors[0]["peptide_ref"]
                else "crosslink_group"
            )
            expected_peptides = 1 if kind == "loop_link_group" else 2
            if len(peptide_ids) != expected_peptides:
                errors.append("source_chain_count_unresolved")
        elif not roles and len(peptide_ids) == 1:
            kind = "single_peptide_group"
        else:
            kind = "unresolved_group"
            errors.append("donor_acceptor_layout_unresolved")

        native_ids = sir.attrib.get("spectrumID", "").split(",")
        if len(native_ids) != 2 or any(
            re.fullmatch(
                r"controllerType=\d+\s+controllerNumber=\d+\s+scan=\d+",
                value,
            ) is None for value in native_ids
        ):
            errors.append("paired_native_id_layout_unresolved")
        if sir.attrib.get("spectraData_ref") not in index["SpectraData"]:
            errors.append("source_spectra_reference_unresolved")

        variants = defaultdict(list)
        for item in items:
            key = tuple(
                item.attrib.get(name) for name in (
                    "calculatedMassToCharge",
                    "experimentalMassToCharge", "chargeState",
                )
            )
            variants[key].append(item)
        if len(variants) != len(native_ids):
            errors.append("native_id_and_precursor_variant_counts_differ")

        variant_views = []
        for number, (key, members) in enumerate(variants.items(), start=1):
            member_refs = [
                member.attrib.get("peptide_ref") for member in members
            ]
            if (
                len(member_refs) != len(peptide_ids)
                or set(member_refs) != set(peptide_ids)
            ):
                errors.append("precursor_variant_chain_membership_unresolved")
            variant_views.append({
                "source_channel_ordinal": number,
                "ordered_native_id_candidate": (
                    native_ids[number - 1]
                    if number <= len(native_ids) else None
                ),
                "calculated_mass_to_charge_raw": key[0],
                "experimental_mass_to_charge_raw": key[1],
                "charge_raw": key[2],
                "item_ids": [member.attrib["id"] for member in members],
                "peptide_refs": member_refs,
            })

        sites = {}
        if kind in ("crosslink_group", "loop_link_group"):
            if not (
                donors[0]["link_token"] == acceptors[0]["link_token"] == token
            ):
                errors.append("source_link_tokens_disagree")
            for name, site in (
                ("donor", donors[0]), ("acceptor", acceptors[0])
            ):
                try:
                    sites[name] = _normalize_site(
                        site, items, evidence, proteins
                    )
                except (KeyError, ValueError) as error:
                    errors.append(name + ":" + str(error))

        record_evidence = list(dict.fromkeys(
            ref.attrib["peptideEvidence_ref"]
            for item in items
            for ref in _children(item, "PeptideEvidenceRef")
        ))
        decoy_flags = [
            evidence[value].get("native_is_decoy")
            for value in record_evidence
        ]
        decoy_status = (
            "unknown" if not decoy_flags or any(
                value is None for value in decoy_flags
            )
            else "source_decoy" if any(decoy_flags)
            else "target_target"
        )
        neutral_shift = None
        if len(variant_views) == 2:
            try:
                first, second = variant_views
                if first["charge_raw"] == second["charge_raw"]:
                    neutral_shift = str(
                        (
                            Decimal(second["calculated_mass_to_charge_raw"])
                            - Decimal(first["calculated_mass_to_charge_raw"])
                        ) * Decimal(first["charge_raw"])
                    )
            except (InvalidOperation, TypeError):
                errors.append("isotope_mass_shift_unresolved")

        native_fdr, fdr_annotations = _native_fdr(items)
        snapshots = [_item_snapshot(item) for item in items]
        errors = list(dict.fromkeys(errors))
        records.append({
            "source_record_id": record_id,
            "source_record_number": ordinal,
            "source_sir_attributes": dict(sir.attrib),
            "source_group_value": token,
            "source_items": snapshots,
            "source_fields": {
                "sir_attributes": dict(sir.attrib),
                "group_value": token,
                "items": snapshots,
            },
            "source_values": [],
            "source_peptide_refs": peptide_ids,
            "source_rank_values": [
                item.attrib.get("rank") for item in items
            ],
            "source_native_ids": native_ids,
            "scan_number": None,
            "record_kind": kind,
            "peptide_multiplicity": (
                1 if kind in ("loop_link_group", "single_peptide_group")
                else 2 if kind == "crosslink_group" else None
            ),
            "record_status": (
                "not_evaluable" if errors
                else "retained_nonpairwise" if kind == "single_peptide_group"
                else "valid"
            ),
            "record_reason": errors[0] if errors else "",
            "source_errors": errors,
            "source_precursor_variants": variant_views,
            "calculated_neutral_isotope_shift": neutral_shift,
            "normalized_sites": sites,
            "source_crosslinker_terms": [
                parameter
                for site in roles
                for parameter in site["modification_parameters"]
                if parameter.get("accession", "").startswith("XLMOD:")
            ],
            "crosslinker_identity_verified": False,
            "native_decoy_status": decoy_status,
            "native_fdr": native_fdr,
            "native_fdr_annotations": fdr_annotations,
            "native_fdr_scope_verified": False,
            "fdr": {
                level: _fdr_value(None) for level in FDR_COLUMNS
            },
        })

    sequence_collection = {
        "db_sequences": [
            {
                "attributes": dict(element.attrib),
                "parameters": _parameters(element),
                "sequence": _sequence(element, "Seq", optional=True),
            }
            for element in index["DBSequence"].values()
        ],
        "peptides": [
            {
                "attributes": dict(element.attrib),
                "parameters": _parameters(element),
                "sequence": _sequence(element, "PeptideSequence"),
                "modifications": [
                    {
                        "attributes": dict(child.attrib),
                        "parameters": _parameters(child),
                    }
                    for child in _children(element, "Modification")
                ],
            }
            for element in index["Peptide"].values()
        ],
        "peptide_evidence": [
            {
                "attributes": dict(element.attrib),
                "parameters": _parameters(element),
            }
            for element in index["PeptideEvidence"].values()
        ],
    }
    if _sha256(results_path) != psm_sha256:
        raise ValueError("Native mzIdentML source changed during import")
    if _sha256(fasta_path) != fasta_sha256:
        raise ValueError("Native FASTA source changed during import")

    return {
        "schema_version": 1,
        "format": "mzidentml_csm",
        "source_profile": "openms_labeled",
        "source_record_level": "crosslink_spectrum_group",
        "source_columns": [
            "source_sir_attributes", "source_group_value", "source_items",
        ],
        "provenance": {
            "dataset_accession": dataset_accession,
            "psm_path": str(results_path.resolve()),
            "psm_sha256": psm_sha256,
            "fasta_path": str(fasta_path.resolve()),
            "fasta_sha256": fasta_sha256,
            "source_file_name": results_path.name,
            "xml_sha256": hashlib.sha256(raw).hexdigest(),
            "source_profile": "openms_labeled",
            "max_xml_bytes": max_xml_bytes,
        },
        "proteins": proteins,
        "records": records,
        "source_identifications": source_identifications,
        "source_sequence_collection": sequence_collection,
        "unresolved_source_items": unresolved_items,
        "coordinate_policy": {
            "source_peptide_evidence_base": 0,
            "source_modification_location_base": 1,
            "normalized_position_base": 1,
            "source_rank_rewritten": False,
        },
        "evidence_validation": evidence,
        "native_metadata": {
            "root_attributes": dict(root.attrib),
            "spectra_data": [
                {
                    "attributes": dict(element.attrib),
                    "parameters": _parameters(element),
                }
                for element in index["SpectraData"].values()
            ],
            "protocol_xml": [
                ET.tostring(element, encoding="unicode")
                for element in root.iter()
                if _tag(element) == "SpectrumIdentificationProtocol"
            ],
        },
        "parser": {
            "backend": "native",
            "name": "openms_labeled",
            "library": None,
            "library_validation_performed": False,
            "parsed_records": 0,
            "source_groups": len(records),
            "source_identification_items": len(
                index["SpectrumIdentificationItem"]
            ),
            "validated_groups": sum(
                record["record_status"] == "valid" for record in records
            ),
            "modification_masses_assessed": False,
        },
        "record_kind_counts": dict(Counter(
            record["record_kind"] for record in records
        )),
        "decoy_policy": {
            "prefix": decoy_prefix,
            "native_is_decoy_preserved": True,
            "target_decoy_validation_performed": False,
        },
        "fdr_policy": {
            "kind": "native_openms_annotations",
            "parameter_name": "OpenXQuest:fdr",
            "level": None,
            "cutoff": None,
            "verification_status": "unverified",
            "per_record_fdr_applied": False,
            "negative_values_clamped": False,
        },
        "fdr_filter_applied": False,
        "identity_mapping_applied": False,
        "input_exported": False,
        "assessment_performed": False,
    }
