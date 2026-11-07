"""Import paired mzIdentML annotations without inferring missing confidence."""

import csv
import gzip
import hashlib
import io
import json
import math
import re
import warnings
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from copy import deepcopy
from importlib.metadata import version
from pathlib import Path

from .xifdr import (
    FDR_COLUMNS,
    _fdr_value,
    _reference_evidence,
    _projection_coordinate_options,
    _project_peptide,
    read_search_fasta,
    write_xifdr_outputs,
)


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


def _params(element):
    return [
        {"parameter_type": _tag(child), **child.attrib}
        for child in element
        if _tag(child) in ("cvParam", "userParam")
    ]


def _boolean(raw):
    return {"true": True, "false": False, "1": True, "0": False}.get(raw)


def _integer(raw):
    if not isinstance(raw, str) or re.fullmatch(r"[1-9]\d*", raw) is None:
        raise ValueError(f"Expected positive source coordinate: {raw!r}")
    return int(raw)


def _scan(raw):
    matches = re.findall(
        r"(?:^|[\s,;])scan=(\d+)(?=$|[\s,;])", str(raw)
    )
    if len(matches) != 1 or int(matches[0]) < 1:
        raise ValueError("No unique explicit scan number")
    return int(matches[0])


def _xml_bytes(path, maximum):
    path = Path(path)
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    try:
        with (gzip.open(path, "rb") if compressed else path.open("rb")) as handle:
            raw = handle.read(maximum + 1)
    except (OSError, EOFError) as error:
        raise ValueError("Cannot read mzIdentML payload") from error
    if len(raw) > maximum:
        raise ValueError("mzIdentML XML exceeds the configured byte limit")
    return raw


def looks_like_mzidentml(path):
    """Recognize the XML root from a bounded plain or gzip prefix."""
    path = Path(path)
    with path.open("rb") as handle:
        signature = handle.read(4)
    if not (
        signature.startswith(b"\x1f\x8b")
        or signature.lstrip(b"\xef\xbb\xbf \r\n\t").startswith(b"<")
        or signature[:1] in (b" ", b"\r", b"\n", b"\t")
    ):
        return False
    try:
        opener = gzip.open if signature.startswith(b"\x1f\x8b") else open
        with opener(path, "rb") as handle:
            prefix = handle.read(65536)
        pull = ET.XMLPullParser(events=("start",))
        pull.feed(prefix)
        event = next(iter(pull.read_events()), None)
        return event is not None and _tag(event[1]) == "MzIdentML"
    except (ET.ParseError, OSError, EOFError):
        return False


def verify_mzidentml_sources(bundle):
    if bundle.get("format") != "mzidentml_csm":
        raise ValueError("Expected an mzIdentML bundle")
    for kind in ("psm", "fasta"):
        provenance = bundle["provenance"]
        if _sha256(provenance[kind + "_path"]) != provenance[kind + "_sha256"]:
            raise ValueError(f"mzIdentML source checksum changed: {kind}")


def _item_snapshot(item):
    return {
        "attributes": dict(item.attrib),
        "parameters": _params(item),
        "peptide_evidence_refs": [
            dict(child.attrib) for child in _children(item, "PeptideEvidenceRef")
        ],
    }


def _native_endpoint(item, peptides, evidence, databases, proteins, aliases):
    peptide_id = item.get("peptide_ref")
    if peptide_id not in peptides:
        raise ValueError("Unresolved peptide reference")
    peptide = peptides[peptide_id]
    sequence = peptide["sequence"]
    sites = []
    for modification in peptide["modifications"]:
        for parameter in modification["parameters"]:
            if parameter.get("accession") in ("MS:1002509", "MS:1002510"):
                sites.append({
                    "role": parameter["accession"],
                    "link_token": parameter.get("value"),
                    "location": modification["attributes"].get("location"),
                    "modification_attributes": modification["attributes"],
                    "modification_parameters": modification["parameters"],
                })
    if len(sites) != 1:
        raise ValueError("Expected one donor or receiver site per peptide")
    site = sites[0]
    position = _integer(site["location"])
    if position > len(sequence) or site["link_token"] in (None, ""):
        raise ValueError("Invalid or missing donor/receiver site")
    associations = []
    references = _children(item, "PeptideEvidenceRef")
    if not references:
        raise ValueError("Missing peptide evidence")
    for reference in references:
        evidence_id = reference.get("peptideEvidence_ref")
        if evidence_id not in evidence:
            raise ValueError("Unresolved peptide evidence")
        entry = evidence[evidence_id]
        if entry.get("peptide_ref") != peptide_id:
            raise ValueError("Peptide evidence refers to another peptide")
        db_id = entry.get("dBSequence_ref")
        if db_id not in databases:
            raise ValueError("Unresolved database sequence")
        database = databases[db_id]
        accession = database["attributes"].get("accession", "")
        start, end = _integer(entry.get("start")), _integer(entry.get("end"))
        if (
            end != start + len(sequence) - 1
            or database["sequence"][start - 1:end] != sequence
        ):
            raise ValueError("Peptide span differs from embedded protein")
        candidates = [
            identifier for identifier in sorted(aliases.get(accession, set()))
            if proteins[identifier]["sequence"] == database["sequence"]
        ]
        associations.append({
            "peptide_evidence_id": evidence_id,
            "peptide_evidence_attributes": deepcopy(entry),
            "db_sequence_id": db_id,
            "protein_accession": accession,
            "peptide_start": start,
            "peptide_end": end,
            "protein_position": start + position - 1,
            "is_decoy_raw": entry.get("isDecoy"),
            "fasta_candidates": candidates,
        })
    if len(associations) != 1 or len(associations[0]["fasta_candidates"]) != 1:
        raise ValueError("Protein evidence does not resolve uniquely in search FASTA")
    association = associations[0]
    identifier = association["fasta_candidates"][0]
    decoy = _boolean(association["is_decoy_raw"])
    return {
        "peptide_id": peptide_id,
        "source_protein_id": association["protein_accession"],
        "source_fasta_id": identifier,
        "peptide": sequence,
        "peptide_start": association["peptide_start"],
        "peptide_end": association["peptide_end"],
        "protein_position": association["protein_position"],
        "link_position": position,
        "residue": sequence[position - 1],
        "kind": "residue",
        "status": "valid",
        "reason": "",
        "source_link_site": deepcopy(site),
        "associations": associations,
        "decoy": decoy,
    }


def _minimal_key(filename, scan, endpoints):
    return (
        filename, scan,
        tuple(sorted(
            (endpoint["peptide"], endpoint["link_position"])
            for endpoint in endpoints
        )),
    )


def _parsed_endpoint_matches(csm, prefix, endpoint):
    if (
        csm[prefix + "_peptide"] != endpoint["peptide"]
        or csm[prefix + "_peptide_crosslink_position"] != endpoint["link_position"]
    ):
        return False
    expected = {
        "_proteins": [endpoint["source_protein_id"]],
        "_proteins_crosslink_positions": [endpoint["protein_position"]],
        "_proteins_peptide_positions": [endpoint["peptide_start"]],
    }
    for suffix, values in expected.items():
        actual = csm[prefix + suffix]
        if actual is not None and sorted(actual) != sorted(values):
            return False
    actual = csm[prefix + "_decoy"]
    if endpoint["decoy"] is None:
        return actual is None
    return actual is None or actual == endpoint["decoy"]


def _read_with_pyxlms(bundle, raw):
    from pyXLMS import parser

    records = [
        record for record in bundle["records"]
        if record["record_status"] == "valid"
        and record["rank_supported_by_parser"]
    ]
    if not records:
        return
    allowed_items = {
        item["attributes"]["id"]
        for record in records for item in record["source_items"]
    }
    # Filter only the parser view; the complete native XML ledger stays intact.
    view = ET.fromstring(raw)
    for parent in list(view.iter()):
        for sir in list(parent):
            if _tag(sir) != "SpectrumIdentificationResult":
                continue
            for item in list(sir):
                if (
                    _tag(item) == "SpectrumIdentificationItem"
                    and item.get("id") not in allowed_items
                ):
                    sir.remove(item)
            if not _children(sir, "SpectrumIdentificationItem"):
                parent.remove(sir)

    info = bundle["parser"]
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always")
            result = parser.read_mzid(
                io.BytesIO(ET.tostring(view, encoding="utf-8")),
                scan_nr_parser=_scan, decoy=None, verbose=1,
            )
            parsed = result.model_dump(mode="json")["crosslink_spectrum_matches"]
        info["warnings"] = [str(item.message) for item in captured]
        info["library_csm_records"] = len(parsed)
        index = defaultdict(list)
        for record in records:
            index[_minimal_key(
                record["spectrum_file"], record["scan_number"],
                list(record["endpoints"].values()),
            )].append(record)
        matched = {}
        for csm in parsed:
            key = _minimal_key(
                csm["spectrum_file"], csm["scan_nr"],
                [
                    {
                        "peptide": csm[prefix + "_peptide"],
                        "link_position": csm[prefix + "_peptide_crosslink_position"],
                    }
                    for prefix in ("alpha", "beta")
                ],
            )
            matches = []
            for record in index.get(key, []):
                first, second = record["endpoints"]["a"], record["endpoints"]["b"]
                direct = (
                    _parsed_endpoint_matches(csm, "alpha", first)
                    and _parsed_endpoint_matches(csm, "beta", second)
                )
                swapped = (
                    _parsed_endpoint_matches(csm, "alpha", second)
                    and _parsed_endpoint_matches(csm, "beta", first)
                )
                if direct or swapped:
                    matches.append((
                        record,
                        "indistinguishable" if direct and swapped
                        else "same" if direct else "swapped",
                    ))
            if len(matches) == 1:
                record, orientation = matches[0]
                identifier = record["source_record_id"]
                matched[identifier] = (
                    orientation if identifier not in matched else None
                )
        for record in records:
            orientation = matched.get(record["source_record_id"])
            if orientation is None:
                record.update(
                    record_status="not_evaluable",
                    record_reason="pyxlms_source_association_unresolved",
                )
            else:
                record["parser_orientation"] = orientation
        info["parsed_records"] = sum(
            record["record_status"] == "valid" for record in records
        )
    except Exception as error:
        for record in records:
            record.update(
                record_status="not_evaluable",
                record_reason="pyxlms_parse_failed",
                parser_error=f"{type(error).__name__}: {error}",
            )


def read_mzidentml_csms(
    psm_path, fasta_path, *, dataset_accession, backend="pyxlms",
    max_xml_bytes=64 * 1024 * 1024,
):
    """Retain all identifications and validate explicitly paired crosslink groups."""
    if backend not in ("native", "pyxlms"):
        raise ValueError("Backend must be native or pyxlms")
    if not isinstance(dataset_accession, str) or not dataset_accession.strip():
        raise ValueError("Dataset accession is required")
    if type(max_xml_bytes) is not int or max_xml_bytes < 1:
        raise ValueError("XML byte limit must be positive")
    psm_path, fasta_path = Path(psm_path).resolve(), Path(fasta_path).resolve()
    raw = _xml_bytes(psm_path, max_xml_bytes)
    try:
        root = ET.fromstring(raw)
    except ET.ParseError as error:
        raise ValueError("Cannot parse mzIdentML XML") from error
    if _tag(root) != "MzIdentML":
        raise ValueError("Input XML is not mzIdentML")
    elements = list(root.iter())

    def unique_elements(name):
        result = {}
        for element in elements:
            if _tag(element) != name:
                continue
            identifier = element.get("id")
            if not identifier or identifier in result:
                raise ValueError(f"Missing or repeated {name} identifier")
            result[identifier] = element
        return result

    databases = {}
    for identifier, element in unique_elements("DBSequence").items():
        text = "".join(
            "".join(child.itertext()) for child in _children(element, "Seq")
        )
        databases[identifier] = {
            "attributes": dict(element.attrib),
            "sequence": "".join(text.split()),
        }
    peptides = {}
    for identifier, element in unique_elements("Peptide").items():
        sequences = _children(element, "PeptideSequence")
        peptides[identifier] = {
            "attributes": dict(element.attrib),
            "sequence": (
                (sequences[0].text or "").strip() if len(sequences) == 1 else ""
            ),
            "modifications": [
                {"attributes": dict(child.attrib), "parameters": _params(child)}
                for child in _children(element, "Modification")
            ],
        }
    evidence = {
        identifier: dict(element.attrib)
        for identifier, element in unique_elements("PeptideEvidence").items()
    }
    spectra = unique_elements("SpectraData")
    proteins = read_search_fasta(fasta_path)
    aliases = defaultdict(set)
    for identifier, protein in proteins.items():
        aliases[identifier].add(identifier)
        if protein["uniprot_accession"]:
            aliases[protein["uniprot_accession"]].add(identifier)

    provenance = {
        "dataset_accession": dataset_accession,
        "psm_path": str(psm_path),
        "psm_sha256": _sha256(psm_path),
        "fasta_path": str(fasta_path),
        "fasta_sha256": _sha256(fasta_path),
        "xml_sha256": hashlib.sha256(raw).hexdigest(),
        "source_file_name": psm_path.name,
        "max_xml_bytes": max_xml_bytes,
    }
    identifications, records = [], []
    for sir in elements:
        if _tag(sir) != "SpectrumIdentificationResult":
            continue
        items = _children(sir, "SpectrumIdentificationItem")
        identifications.append({
            "attributes": dict(sir.attrib),
            "parameters": _params(sir),
            "items": [_item_snapshot(item) for item in items],
        })
        grouped = defaultdict(list)
        for item in items:
            markers = [
                value for value in _params(item)
                if value.get("accession") == "MS:1002511"
            ]
            if not markers:
                continue
            if len(markers) != 1 or markers[0].get("value") in (None, ""):
                grouped["unresolved-item:" + str(item.get("id"))].append(item)
            else:
                grouped[markers[0]["value"]].append(item)
        for group_value, members in grouped.items():
            identity = json.dumps(
                [sir.get("id"), group_value],
                ensure_ascii=False, separators=(",", ":"),
            )
            source_items = [_item_snapshot(item) for item in members]
            record = {
                "source_record_id": (
                    dataset_accession + ":" + provenance["psm_sha256"] + ":"
                    + hashlib.sha256(identity.encode("utf-8")).hexdigest()
                ),
                "source_record_number": len(records) + 1,
                "source_sir_attributes": dict(sir.attrib),
                "source_group_value": group_value,
                "source_items": source_items,
                "source_fields": {
                    "sir_attributes": dict(sir.attrib),
                    "group_value": group_value,
                    "items": deepcopy(source_items),
                },
                "source_values": [],
                "record_kind": "paired_crosslink_candidate",
                "record_status": "not_evaluable",
                "record_reason": "",
                "decoy_status": "unknown",
                "fdr": {key: _fdr_value(None) for key in FDR_COLUMNS},
                "endpoints": {},
                "source_crosslinker": None,
                "source_crosslinker_terms": [],
                "native_pass_thresholds": [
                    item.get("passThreshold") for item in members
                ],
            }
            records.append(record)
            try:
                if group_value.startswith("unresolved-item:"):
                    raise ValueError("Ambiguous crosslink group annotation")
                if len(members) != 2:
                    raise ValueError("Crosslink group does not have two items")
                endpoints = [
                    _native_endpoint(
                        item, peptides, evidence, databases, proteins, aliases
                    )
                    for item in members
                ]
                first, second = [
                    endpoint["source_link_site"] for endpoint in endpoints
                ]
                if (
                    {first["role"], second["role"]}
                    != {"MS:1002509", "MS:1002510"}
                    or first["link_token"] != second["link_token"]
                ):
                    raise ValueError("Donor and receiver linkage tokens do not pair")
                record["endpoints"] = dict(zip(("a", "b"), endpoints))
                flags = [
                    _boolean(association["is_decoy_raw"])
                    for endpoint in endpoints
                    for association in endpoint["associations"]
                ]
                record["decoy_status"] = (
                    "source_decoy" if any(value is True for value in flags)
                    else "target_target" if flags and all(
                        value is False for value in flags
                    ) else "unknown"
                )
                terms = sorted({
                    (value.get("accession"), value.get("name", ""))
                    for endpoint in endpoints
                    for value in endpoint["source_link_site"]["modification_parameters"]
                    if value.get("accession", "").startswith("XLMOD:")
                })
                record["source_crosslinker_terms"] = [
                    {"accession": accession, "name": name}
                    for accession, name in terms
                ]
                if len(terms) == 1 and terms[0][1]:
                    record["source_crosslinker"] = terms[0][1]
                spectrum = spectra.get(sir.get("spectraData_ref"))
                if spectrum is None or not spectrum.get("location"):
                    raise ValueError("Missing source spectrum file location")
                record["spectrum_file"] = spectrum.get("location")
                record["scan_number"] = _scan(sir.get("spectrumID"))
                record["rank_supported_by_parser"] = all(
                    item.get("rank") is None or int(item.get("rank")) == 1
                    for item in members
                )
                if not record["rank_supported_by_parser"]:
                    record["record_reason"] = "rank_not_supported_by_parser"
                else:
                    record["record_status"] = "valid"
            except (ValueError, KeyError, TypeError) as error:
                record.update(
                    record_reason="source_group_unresolved",
                    source_error=str(error),
                )

    thresholds = [
        [
            dict(child.attrib) for child in element.iter()
            if _tag(child) in ("cvParam", "userParam")
        ]
        for element in elements if _tag(element) == "Threshold"
    ]
    bundle = {
        "format": "mzidentml_csm",
        "source_columns": [],
        "source_representation": "mzidentml_xml_groups",
        "records": records,
        "proteins": proteins,
        "source_identifications": identifications,
        "source_sequence_collection": {
            "database_sequences": databases,
            "peptides": peptides,
            "peptide_evidence": evidence,
        },
        "native_metadata": {
            "root_tag": root.tag,
            "root_attributes": dict(root.attrib),
            "threshold_parameters": thresholds,
            "spectra_data": {
                identifier: dict(element.attrib)
                for identifier, element in spectra.items()
            },
            "analysis_software": [
                {
                    "attributes": dict(element.attrib),
                    "parameters": [
                        dict(child.attrib) for child in element.iter()
                        if _tag(child) in ("cvParam", "userParam")
                    ],
                }
                for element in elements if _tag(element) == "AnalysisSoftware"
            ],
        },
        "provenance": provenance,
        "fdr_policy": {
            "kind": "mzidentml_source_annotations",
            "verification_status": "unverified",
            "per_record_fdr_applied": False,
            "pass_threshold_used_as_fdr": False,
            "native_pass_threshold_policy": "exclude_explicit_failure",
            "native_threshold_parameters": thresholds,
        },
        "parser": {
            "name": "pyXLMS" if backend == "pyxlms" else "native",
            "version": version("pyXLMS") if backend == "pyxlms" else None,
            "parsed_records": 0,
            "modification_masses_assessed": False,
            "input_policy": "validated_rank_one_source_groups",
        },
    }
    if backend == "pyxlms":
        _read_with_pyxlms(bundle, raw)
    verify_mzidentml_sources(bundle)
    return bundle


def _verify_review(review, source, profile, kind):
    if not isinstance(review, dict):
        raise ValueError("Candidate reference requires a bound review")
    expected = {
        "schema_version": 1,
        "review_status": "reviewed",
        "kind": kind,
        "source_fasta_id": source["id"],
        "source_fasta_sha256": source["fasta_sha256"],
        "source_sequence_sha256": hashlib.sha256(
            source["sequence"].encode("ascii")
        ).hexdigest(),
        "declared_accession": source["uniprot_accession"],
        "reference_accession": profile["accession"],
        "reference_sequence_sha256": profile["reference_sequence_sha256"],
        "taxon_id": profile["taxon_id"],
    }
    if any(review.get(key) != value for key, value in expected.items()):
        raise ValueError("Candidate reference review bindings differ")
    if not isinstance(review.get("reason"), str) or not review["reason"].strip():
        raise ValueError("Reference review requires a reason")
    receipts = review.get("sources")
    if not isinstance(receipts, list) or not receipts:
        raise ValueError("Reference review evidence is missing")
    declaration = {key: value for key, value in review.items() if key != "sources"}
    for receipt in receipts:
        path = Path(receipt["absolute_path"])
        if _sha256(path) != receipt["sha256"]:
            raise ValueError("Reference review evidence checksum changed")
        if json.loads(path.read_text(encoding="utf-8-sig")) != declaration:
            raise ValueError("Reference review document differs from declaration")


def _reference_profiles(bundle, reference_mapping):
    if (
        reference_mapping.get("kind") != "mzidentml_reference_mapping"
        or reference_mapping.get("schema_version") != 1
        or reference_mapping.get("source_fasta_sha256")
        != bundle["provenance"]["fasta_sha256"]
    ):
        raise ValueError("Reference mapping differs from source FASTA")
    references = reference_mapping.get("references")
    if not isinstance(references, dict) or not references:
        raise ValueError("Per-protein references are required")
    profiles, options, mappings = {}, {}, {}
    for identifier, supplied in references.items():
        source = bundle["proteins"].get(identifier)
        if source is None:
            raise ValueError("Unknown source protein in reference mapping")
        profile = deepcopy(supplied)
        reference = profile.get("reference_sequence")
        if (
            profile.get("source_fasta_id") != identifier
            or profile.get("source_sequence") != source["sequence"]
            or not isinstance(reference, str)
            or re.fullmatch(r"[A-Z]+", reference) is None
            or profile.get("reference_sequence_sha256")
            != hashlib.sha256(reference.encode("ascii")).hexdigest()
        ):
            raise ValueError("Reference profile differs from source or sequence")
        if profile["accession"] != source["uniprot_accession"]:
            _verify_review(
                profile.get("reference_accession_review"),
                source, profile, "reference_accession_review",
            )
        if reference != source["sequence"]:
            review = profile.get("sequence_projection_review")
            _verify_review(review, source, profile, "sequence_projection_review")
            if (
                review.get("require_all_optimal_alignments") is not True
                or type(review.get("max_optimal_alignments")) is not int
                or not 1 <= review["max_optimal_alignments"] <= 64
            ):
                raise ValueError("Invalid sequence projection review policy")
            maximum = review["max_optimal_alignments"]
        else:
            maximum = 64
        _reference_evidence(profile, profile, source)
        coordinate_options, evidence = _projection_coordinate_options(
            source["sequence"], reference, maximum
        )
        if reference != source["sequence"]:
            if (
                evidence["method"] != "full_reference_exact_subsequence"
                or len(coordinate_options) != 1
            ):
                raise ValueError("Sequence review requires one intact contained reference")
            unmapped = sorted(
                set(range(1, len(source["sequence"]) + 1))
                - set(coordinate_options[0])
            )
            if review.get("allowed_source_unmapped_positions") != unmapped:
                raise ValueError("Sequence review unmapped positions differ")
        exact = (
            reference == source["sequence"]
            and profile["accession"] == source["uniprot_accession"]
        )
        scope = "exact_sequence" if exact else "candidate_hypothesis"
        if profile.get("identity_scope") != scope:
            raise ValueError("Reference identity scope differs from its evidence")
        profile["coordinate_evidence"] = evidence
        profiles[identifier] = profile
        options[identifier] = coordinate_options
        receipts = deepcopy(profile["sources"])
        for key in ("reference_accession_review", "sequence_projection_review"):
            if key in profile:
                receipts.extend(deepcopy(profile[key]["sources"]))
        mappings[identifier] = {
            "status": "mapped" if exact else "reference_projected",
            "accession": profile["accession"],
            "sequence": reference,
            "taxon_id": profile["taxon_id"],
            "identity_scope": scope,
            "sources": receipts,
        }
    return profiles, options, mappings


def write_mzidentml_claudio_input(
    bundle, reference_mapping, output_dir, *,
    allowed_crosslinkers, allowed_residues,
    export_unverified_candidates=False,
):
    """Export explicit candidates while retaining unknown decoys and native flags."""
    if export_unverified_candidates is not True:
        raise ValueError("Explicit unverified candidate export mode is required")
    if Path(output_dir).exists():
        raise FileExistsError("Output directory already exists")
    verify_mzidentml_sources(bundle)
    selections = {}
    for name, values in (
        ("crosslinkers", allowed_crosslinkers), ("residues", allowed_residues),
    ):
        if (
            not isinstance(values, (set, frozenset, tuple, list))
            or not values
            or any(not isinstance(value, str) or not value.strip() for value in values)
        ):
            raise ValueError("Explicit nonempty " + name + " are required")
        selections[name] = set(values)
    if any(re.fullmatch(r"[A-Z]", value) is None for value in selections["residues"]):
        raise ValueError("Residues must be single uppercase letters")
    provenance = bundle["provenance"]
    native = read_mzidentml_csms(
        provenance["psm_path"], provenance["fasta_path"],
        dataset_accession=provenance["dataset_accession"], backend="native",
        max_xml_bytes=provenance["max_xml_bytes"],
    )
    for key in (
        "proteins", "source_identifications", "source_sequence_collection",
        "native_metadata", "fdr_policy",
    ):
        if native[key] != bundle[key]:
            raise ValueError("mzIdentML bundle differs from native sources")
    if len(native["records"]) != len(bundle["records"]):
        raise ValueError("Native crosslink group count changed")
    for original, checked in zip(bundle["records"], native["records"]):
        for key in (
            "source_record_id", "source_record_number",
            "source_sir_attributes", "source_group_value", "source_items",
            "source_fields", "record_kind", "decoy_status", "fdr",
            "endpoints", "native_pass_thresholds", "source_crosslinker",
            "source_crosslinker_terms",
        ):
            if original[key] != checked[key]:
                raise ValueError("Native crosslink evidence changed after import")
        if original["record_status"] == "valid" and checked["record_status"] != "valid":
            raise ValueError("Native crosslink validation changed")

    profiles, options, mappings = _reference_profiles(bundle, reference_mapping)
    rows, audit = [], []
    for record in bundle["records"]:
        entry = deepcopy(record)
        reasons = []
        if record["record_status"] != "valid":
            reasons.append(record["record_reason"] or "source_record_unresolved")
        if record["decoy_status"] == "source_decoy":
            reasons.append("decoy_record")
        if any(_boolean(raw) is False for raw in record["native_pass_thresholds"]):
            reasons.append("native_pass_threshold_failed")
        crosslinker = record["source_crosslinker"]
        if not crosslinker:
            reasons.append("crosslinker_annotation_unresolved")
        elif crosslinker not in selections["crosslinkers"]:
            reasons.append("crosslinker_not_selected")

        projected = {}
        if not reasons:
            for site in ("a", "b"):
                endpoint = record["endpoints"][site]
                identifier = endpoint["source_fasta_id"]
                profile = profiles.get(identifier)
                if profile is None:
                    reasons.append("protein_identity_unresolved_" + site)
                    continue
                result = _project_peptide(endpoint, profile, options[identifier])
                projected[site] = result
                if result["status"] != "projected":
                    reasons.append("reference_projection_unresolved_" + site)
                    continue
                if result["compatible_alignment_count"] != len(options[identifier]):
                    reasons.append("reference_projection_not_invariant_" + site)
                    continue
                residue = profile["reference_sequence"][
                    result["reference_protein_position"] - 1
                ]
                if residue not in selections["residues"]:
                    reasons.append("residue_not_selected_" + site)

        entry.update(
            input_status="not_evaluable" if reasons else "candidate_ready",
            input_reasons=list(dict.fromkeys(reasons)),
            reference_endpoints=deepcopy(projected),
        )
        audit.append(entry)
        if reasons:
            continue

        scope = (
            "exact_sequence" if all(
                profiles[record["endpoints"][site]["source_fasta_id"]][
                    "identity_scope"
                ] == "exact_sequence"
                for site in ("a", "b")
            ) else "candidate_hypothesis"
        )
        row = {
            "source_record_id": record["source_record_id"],
            "source_record_number": record["source_record_number"],
            "dataset_accession": provenance["dataset_accession"],
            "source_psm_sha256": provenance["psm_sha256"],
            "source_fasta_sha256": provenance["fasta_sha256"],
            "source_reported_file_name": provenance["source_file_name"],
            "source_sir_id": record["source_sir_attributes"].get("id"),
            "source_group_value": record["source_group_value"],
            "source_crosslinker": crosslinker,
            "source_crosslinker_terms_json": json.dumps(
                record["source_crosslinker_terms"], sort_keys=True
            ),
            "source_decoy_status": record["decoy_status"],
            "source_pass_thresholds_json": json.dumps(
                record["native_pass_thresholds"]
            ),
            "source_spectrum_file": record["spectrum_file"],
            "source_scan_number": record["scan_number"],
            "source_items_json": json.dumps(record["source_items"], sort_keys=True),
            "input_scope": "unverified_candidate",
            "identity_scope": scope,
            "fdr_verification_status": "unverified",
            "selected_fdr_level": None,
            "selected_fdr_value": None,
            "selected_fdr_cutoff": None,
        }
        for site in ("a", "b"):
            endpoint = record["endpoints"][site]
            result = projected[site]
            profile = profiles[endpoint["source_fasta_id"]]
            row.update({
                "pep_" + site: endpoint["peptide"],
                "pos_" + site: result["reference_protein_position"],
                "res_pos_" + site: endpoint["link_position"],
                "unip_id_" + site: profile["accession"],
                "seq_" + site: profile["reference_sequence"],
                "source_protein_id_" + site: endpoint["source_protein_id"],
                "source_fasta_id_" + site: endpoint["source_fasta_id"],
                "source_protein_position_" + site: endpoint["protein_position"],
                "source_peptide_start_" + site: endpoint["peptide_start"],
            })
        rows.append(row)

    bodies = {}
    for kind in ("psm", "fasta"):
        path = Path(provenance[kind + "_path"])
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != provenance[kind + "_sha256"]:
            raise ValueError("Native source checksum changed before export")
        suffix = ".mzid.gz" if kind == "psm" and raw[:2] == b"\x1f\x8b" else path.suffix
        bodies["native_sources/" + provenance[kind + "_sha256"] + suffix] = raw
    for name, data in (
        ("identification_audit.json", bundle["source_identifications"]),
        ("sequence_collection.json", bundle["source_sequence_collection"]),
        ("native_metadata.json", bundle["native_metadata"]),
    ):
        bodies[name] = (
            json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        ).encode("utf-8")
    prepared = {
        "rows": rows,
        "audit": audit,
        "fdr_policy": deepcopy(bundle["fdr_policy"]),
        "extra_source_bodies": bodies,
        "reference_projection": deepcopy(reference_mapping),
        "reference_projection_policy": {
            "kind": "source_bound_mzidentml_reference_mapping",
            "require_all_optimal_alignments": True,
        },
    }
    return write_xifdr_outputs(
        bundle, mappings, prepared, output_dir,
        pipeline_policy={
            "input_scope": "unverified_candidate",
            "crosslinkers": sorted(selections["crosslinkers"]),
            "residues": sorted(selections["residues"]),
            "native_pass_threshold_policy": "exclude_explicit_failure",
            "missing_decoy_policy": "retain_unknown",
            "coordinate_basis": 1,
        },
    )
