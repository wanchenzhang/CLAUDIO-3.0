"""Import MaxLynx semicolon-localized CSM tables and export candidates."""

import csv
import hashlib
import importlib
import json
import math
import re
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from copy import deepcopy
from importlib.metadata import version
from pathlib import Path
from unittest.mock import patch

from .xifdr import (
    FDR_COLUMNS,
    _fdr_value,
    _positive_integer,
    _merox_multi_reference_profiles,
    _project_peptide,
    read_search_fasta,
    write_xifdr_outputs,
)


REQUIRED_COLUMNS = {
    "Raw file", "Scan number", "Charge", "Sequence1", "Sequence2",
    "Modified sequence1", "Modified sequence2", "Proteins1", "Proteins2",
    "Decoy", "Score", "Partial score 1", "Partial score 2", "Delta score",
    "Crosslink Product Type", "PEP", "InterLinks",
    "Pro_InterLink1", "Pep_InterLink1", "AA_InterLink1",
    "Pro_InterLink2", "Pep_InterLink2", "AA_InterLink2",
}
LOCALIZATION_FIELDS = (
    "Pro_InterLink1", "Pep_InterLink1", "AA_InterLink1",
    "Pro_InterLink2", "Pep_InterLink2", "AA_InterLink2", "InterLinks",
)
NONPAIR_REASONS = {
    "SingleProtein": "single_peptide_identification_not_pairwise",
    "Monolinked": "monolink_not_pairwise_crosslink",
    "Looplinked": "loop_link_requires_single_peptide_support",
}
KINDS = {
    "SingleProtein": "linear_peptide",
    "Monolinked": "dead_end",
    "Looplinked": "loop_link",
    "IntraProXL": "paired_crosslink_candidate",
    "InterProXL": "paired_crosslink_candidate",
}


NATIVE_57_COLUMNS = frozenset(('Raw file', 'Scan number', 'Scan index', 'Sequence1', 'Sequence2', 'Modified sequence1', 'Modified sequence2', 'Length1', 'Length2', '#Matches1', '#Matches2', 'Decoy', 'Score', 'Crosslink product type', 'PEP', 'Search mode', 'Proteins1', 'Proteins2', 'Pair of linked residues', 'Protein index of Crosslink 1', 'Peptide index of Crosslink 1', 'AA name of Crosslink 1', 'Protein index of Crosslink 2', 'Peptide index of Crosslink 2', 'AA name of Crosslink 2', 'Protein index of Looplink 1 (first residue)', 'Peptide index of Looplink 1 (first residue)', 'AA name of Looplink 1 (first residue)', 'Protein index of Looplink 1 (second residue)', 'Peptide index of Looplink 1 (second residue)', 'AA name of Looplink 1 (second residue)', 'Protein index of Looplink 2 (first residue)', 'Peptide index of Looplink 2 (first residue)', 'AA name of Looplink 2 (first residue)', 'Protein index of Looplink 2 (second residue)', 'Peptide index of Looplink 2 (second residue)', 'AA name of Looplink 2 (second residue)', 'Protein index of Monolink 1', 'Peptide index of Monolink 1', 'AA name of of Monolink 1', 'Protein index of Monolink 2', 'Peptide index of Monolink 2', 'AA name of of Monolink 2', 'All possible linking sites on peptide 1', 'All possible linking sites on peptide 2', 'Delta score', 'Partial score 1', 'Partial score 2', 'Charge', 'm/z', 'Mass', 'Mass error [ppm]', 'Annotation ions', 'Annotation m/z', 'Query Type', '#Overlapped fragments 1', '#Overlapped fragments 2'))

NATIVE_57_COLUMN_MAP = {
    "Crosslink product type": "Crosslink Product Type",
    "Pair of linked residues": "InterLinks",
    "Protein index of Crosslink 1": "Pro_InterLink1",
    "Peptide index of Crosslink 1": "Pep_InterLink1",
    "AA name of Crosslink 1": "AA_InterLink1",
    "Protein index of Crosslink 2": "Pro_InterLink2",
    "Peptide index of Crosslink 2": "Pep_InterLink2",
    "AA name of Crosslink 2": "AA_InterLink2",
}
NATIVE_57_PRODUCT_MAP = {
    "Linear": "SingleProtein",
    "Mono-linked": "Monolinked",
    "Loop-linked": "Looplinked",
    "Intra-protein link": "IntraProXL",
    "Inter-protein link": "InterProXL",
    "Mono- and loop-linked": "MonoAndLooplinked",
}
KINDS["MonoAndLooplinked"] = "mono_and_loop_link"
NONPAIR_REASONS["MonoAndLooplinked"] = "mono_and_loop_linked_not_pairwise"


def _maxlynx_header_profile(columns):
    """Recognize reviewed layouts without modifying source columns."""
    names = set(columns)
    if len(columns) != len(names):
        raise ValueError("Unsupported or ambiguous MaxLynx table header")
    if {
        "Crosslink Product Type", "Crosslink product type"
    } <= names:
        raise ValueError("Ambiguous MaxLynx product type columns")
    if REQUIRED_COLUMNS <= names:
        return "semicolon_localized"
    if names == NATIVE_57_COLUMNS:
        return "native_57_column"
    raise ValueError("Unsupported or ambiguous MaxLynx table header")


def _maxlynx_validation_fields(fields, profile):
    """Create a validation view while retaining native fields separately."""
    view = deepcopy(fields)
    if profile == "native_57_column":
        for native, canonical in NATIVE_57_COLUMN_MAP.items():
            view[canonical] = fields.get(native, "")
        product = fields.get("Crosslink product type", "")
        view["Crosslink Product Type"] = NATIVE_57_PRODUCT_MAP.get(
            product, product
        )
    return view


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_maxlynx_sources(bundle):
    """Reject changed native inputs before creating output files."""
    if bundle.get("format") != "maxlynx_csm":
        raise ValueError("Expected a MaxLynx CSM bundle")
    for kind in ("psm", "fasta", "parameters"):
        provenance = bundle["provenance"]
        path = Path(provenance[kind + "_path"])
        if _sha256(path) != provenance[kind + "_sha256"]:
            raise ValueError(f"MaxLynx source checksum changed: {kind}")


def _settings(path):
    try:
        document = ET.fromstring(Path(path).read_bytes())
    except ET.ParseError as error:
        raise ValueError("Cannot parse MaxQuant parameters XML") from error

    def values(tag):
        return [
            (element.text or "").strip()
            for element in document.iter()
            if element.tag.rsplit("}", 1)[-1] == tag
            and (element.text or "").strip()
        ]

    crosslinkers = sorted(set(values("crossLinker")))
    rules = sorted(set(values("identifierParseRule")))
    if len(crosslinkers) != 1 or len(rules) != 1:
        raise ValueError(
            "This layout requires one crosslinker and one identifier parse rule"
        )
    try:
        rule = re.compile(rules[0])
    except re.error as error:
        raise ValueError("Invalid FASTA identifier parse rule") from error
    if rule.groups < 1:
        raise ValueError("FASTA identifier rule requires a capture group")

    raw_cutoffs = values("psmFdrCrosslink")
    cutoff = None
    if len(set(raw_cutoffs)) == 1:
        try:
            candidate = float(raw_cutoffs[0])
            if math.isfinite(candidate) and 0 <= candidate <= 1:
                cutoff = candidate
        except ValueError:
            pass

    settings = {
        "crosslinker": crosslinkers[0],
        "identifier_parse_rule": rules[0],
        "maxquant_versions": values("maxQuantVersion"),
        "decoy_mode": values("decoyMode"),
        "fasta_paths": values("fastaFilePath"),
        "raw_psm_fdr_crosslink": raw_cutoffs,
    }
    policy = {
        "kind": "maxquant_saved_configuration",
        "configured_level": "psm",
        "configured_cutoff": cutoff,
        "verification_status": (
            "configured_only" if cutoff is not None else "unverified"
        ),
        "per_record_fdr_applied": False,
        "application_verified": False,
        "target_decoy_validation_performed": False,
        "source_pep_interpretation": "retained_uninterpreted",
    }
    return settings, policy, rule


def _tokens(raw):
    values = raw.strip().split(";")
    if values and values[-1] == "":
        values.pop()
    values = [value.strip() for value in values]
    if not values or any(not value for value in values):
        raise ValueError("Missing or internally empty localization token")
    return values


def _localizations(fields, proteins, aliases):
    lists = {name: _tokens(fields[name]) for name in LOCALIZATION_FIELDS}
    if len({len(values) for values in lists.values()}) != 1:
        raise ValueError("Localization lists have unequal lengths")
    candidates = []
    for index in range(len(lists["InterLinks"])):
        endpoints = {}
        for side, site in (("1", "a"), ("2", "b")):
            annotation = fields[f"Proteins{side}"].strip()
            matches = sorted(aliases.get(annotation, set()))
            if len(matches) != 1:
                raise ValueError("Protein annotation does not resolve uniquely")
            peptide = fields[f"Sequence{side}"].strip()
            if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", peptide) is None:
                raise ValueError("Unsupported peptide sequence")
            protein_position = _positive_integer(
                lists[f"Pro_InterLink{side}"][index]
            )
            link_position = _positive_integer(
                lists[f"Pep_InterLink{side}"][index]
            )
            residue = lists[f"AA_InterLink{side}"][index]
            start = protein_position - link_position + 1
            end = start + len(peptide) - 1
            sequence = proteins[matches[0]]["sequence"]
            if not (
                1 <= link_position <= len(peptide)
                and 1 <= start <= end <= len(sequence)
                and sequence[start - 1:end] == peptide
                and peptide[link_position - 1] == residue
            ):
                raise ValueError("Source peptide or site differs from search FASTA")
            endpoints[site] = {
                "source_protein_id": annotation,
                "source_fasta_id": matches[0],
                "source_peptide": fields[f"Sequence{side}"],
                "peptide": peptide,
                "peptide_start": start,
                "peptide_end": end,
                "protein_position": protein_position,
                "link_position": link_position,
                "residue": residue,
                "kind": "residue",
                "status": "valid",
                "reason": "",
            }
        match = re.fullmatch(
            r"([A-Z])([1-9]\d*)-([A-Z])([1-9]\d*)",
            lists["InterLinks"][index],
        )
        expected = (
            endpoints["a"]["residue"], endpoints["a"]["protein_position"],
            endpoints["b"]["residue"], endpoints["b"]["protein_position"],
        )
        if match is None or (
            match[1], int(match[2]), match[3], int(match[4])
        ) != expected:
            raise ValueError("InterLinks differs from aligned endpoint lists")
        candidates.append({
            "candidate_number": index + 1,
            "endpoints": endpoints,
        })
    return candidates


def _read_with_pyxlms(bundle):
    from pyXLMS import parser
    import pandas as pd

    eligible = [
        record for record in bundle["records"]
        if record["record_status"] == "valid"
    ]
    info = bundle["parser"]
    if not eligible:
        return

    marker = "__claudio_maxlynx_source_id"
    while marker in bundle["source_columns"]:
        marker += "_"
    views = []
    missing_numeric = 0
    for record in eligible:
        view = _maxlynx_validation_fields(
            record["source_fields"],
            _maxlynx_header_profile(bundle["source_columns"]),
        )
        for side, site in (("1", "a"), ("2", "b")):
            endpoint = record["endpoints"][site]
            view[f"Protein index of Crosslink {side}"] = endpoint["protein_position"]
            view[f"Peptide index of Crosslink {side}"] = endpoint["link_position"]
        for name in ("Score", "Partial score 1", "Partial score 2", "Delta score"):
            if not str(view[name]).strip():
                view[name] = float("nan")
                missing_numeric += 1
        view[marker] = record["source_record_id"]
        views.append(view)
    info["missing_numeric_adapter_values"] = missing_numeric
    frame = pd.DataFrame(views)
    module = importlib.import_module(parser.read_maxquant.__module__)

    def source_endpoint(record, site):
        endpoint = record["endpoints"][site]
        return (
            endpoint["peptide"], endpoint["link_position"],
            (endpoint["source_protein_id"],),
            (endpoint["protein_position"],),
            (endpoint["peptide_start"],), False,
        )

    def parsed_endpoint(csm, prefix):
        return (
            csm[prefix + "_peptide"],
            csm[prefix + "_peptide_crosslink_position"],
            tuple(csm[prefix + "_proteins"] or []),
            tuple(csm[prefix + "_proteins_crosslink_positions"] or []),
            tuple(csm[prefix + "_proteins_peptide_positions"] or []),
            csm[prefix + "_decoy"],
        )

    try:
        with patch.object(module.pd, "read_csv", return_value=frame):
            parsed = parser.read_maxlynx(
                bundle["provenance"]["psm_path"],
                crosslinker=bundle["native_settings"]["crosslinker"],
                parse_modifications=False,
            ).model_dump(mode="json")["crosslink_spectrum_matches"]
        expected = {record["source_record_id"]: record for record in eligible}
        seen = []
        orientations = {}
        for csm in parsed:
            source_id = csm["additional_information"]["source"][marker]
            if source_id not in expected or source_id in seen:
                raise ValueError("Unknown or repeated parser source association")
            seen.append(source_id)
            record = expected[source_id]
            native = (
                source_endpoint(record, "a"), source_endpoint(record, "b")
            )
            actual = (
                parsed_endpoint(csm, "alpha"), parsed_endpoint(csm, "beta")
            )
            direct, swapped = actual == native, actual[::-1] == native
            if not (direct or swapped):
                raise ValueError("Parsed endpoints differ from source endpoints")
            fields = record["source_fields"]
            if (
                csm["spectrum_file"] != fields["Raw file"].strip()
                or csm["scan_nr"] != _positive_integer(fields["Scan number"])
                or csm["charge"] != _positive_integer(fields["Charge"])
            ):
                raise ValueError("Parsed spectrum identity differs from source")
            orientations[source_id] = (
                "indistinguishable" if direct and swapped
                else "same" if direct else "swapped"
            )
        if set(seen) != set(expected):
            raise ValueError("Parser omitted source records")
        for record in eligible:
            record["parser_orientation"] = orientations[record["source_record_id"]]
        info["parsed_records"] = len(parsed)
        info["return_order_matches_source"] = seen == list(expected)
    except Exception as error:
        for record in eligible:
            record.update(
                record_status="not_evaluable",
                record_reason="pyxlms_parse_failed",
                parser_error=f"{type(error).__name__}: {error}",
            )


def read_maxlynx_csms(
    psm_path, fasta_path, *, dataset_accession, parameters_path,
    backend="pyxlms",
):
    """Retain every native row and validate single paired localizations."""
    if backend not in ("native", "pyxlms"):
        raise ValueError("Backend must be native or pyxlms")
    if not isinstance(dataset_accession, str) or not dataset_accession.strip():
        raise ValueError("Dataset accession is required")
    paths = {
        "psm": Path(psm_path).resolve(),
        "fasta": Path(fasta_path).resolve(),
        "parameters": Path(parameters_path).resolve(),
    }
    provenance = {"dataset_accession": dataset_accession}
    for kind, path in paths.items():
        provenance[kind + "_path"] = str(path)
        provenance[kind + "_sha256"] = _sha256(path)

    settings, policy, rule = _settings(paths["parameters"])
    proteins = read_search_fasta(paths["fasta"])
    aliases = defaultdict(set)
    for identifier, protein in proteins.items():
        aliases[identifier].add(identifier)
        match = rule.match(">" + protein["header"])
        if match is not None and match.group(1):
            aliases[match.group(1)].add(identifier)

    with paths["psm"].open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        columns = next(reader, [])
        profile = _maxlynx_header_profile(columns)
        records = []
        for values in reader:
            if not any(value.strip() for value in values):
                continue
            fields = dict(zip(columns, values))
            validation_fields = _maxlynx_validation_fields(fields, profile)
            number = len(records) + 1
            label = fields.get("Decoy", "").strip().casefold()
            status = (
                "target_target" if label == "forward"
                else "source_decoy" if label == "decoy" else "unknown"
            )
            product = validation_fields.get("Crosslink Product Type", "")
            record = {
                "source_record_id": (
                    f"{dataset_accession}:{provenance['psm_sha256']}:{number:06d}"
                ),
                "source_record_number": number,
                "source_values": values,
                "source_fields": fields,
                "record_kind": KINDS.get(product, "unsupported_product_type"),
                "decoy_status": status,
                "fdr": {key: _fdr_value(None) for key in FDR_COLUMNS},
                "source_pep": {
                    "raw": fields.get("PEP"),
                    "interpretation": "unverified",
                },
                "endpoints": {},
                "localization": {"status": "not_applicable", "candidates": []},
                "record_status": "not_evaluable",
                "record_reason": "",
            }
            records.append(record)
            if len(values) != len(columns):
                record["record_reason"] = "row_width_mismatch"
                continue
            if product in NONPAIR_REASONS:
                record["record_reason"] = NONPAIR_REASONS[product]
                continue
            if record["record_kind"] != "paired_crosslink_candidate":
                record["record_reason"] = "unsupported_product_type"
                continue
            if status != "target_target":
                record["record_reason"] = (
                    "decoy_record" if status == "source_decoy"
                    else "unknown_native_decoy_label"
                )
                continue
            if any("REV__" in fields[f"Proteins{side}"] for side in ("1", "2")):
                record["record_reason"] = "conflicting_target_decoy_annotations"
                continue
            try:
                candidates = _localizations(
                    validation_fields, proteins, aliases
                )
                record["localization"] = {
                    "status": (
                        "single_candidate" if len(candidates) == 1
                        else "multiple_candidates"
                    ),
                    "candidates": candidates,
                }
                if len(candidates) != 1:
                    record["record_reason"] = "multiple_source_localizations"
                    continue
                record["endpoints"] = deepcopy(candidates[0]["endpoints"])
                _positive_integer(fields["Scan number"])
                _positive_integer(fields["Charge"])
                if not fields["Raw file"].strip():
                    raise ValueError("Missing spectrum file label")
                record["record_status"] = "valid"
            except (ValueError, KeyError) as error:
                record.update(
                    record_reason="source_endpoint_unresolved",
                    source_error=str(error),
                )

    bundle = {
        "format": "maxlynx_csm",
        "format_variant": "semicolon_interlinks",
        "source_columns": columns,
        "records": records,
        "proteins": proteins,
        "provenance": provenance,
        "native_settings": settings,
        "fdr_policy": policy,
        "parser": {
            "name": "pyXLMS" if backend == "pyxlms" else "native",
            "version": version("pyXLMS") if backend == "pyxlms" else None,
            "entrypoint": "read_maxlynx",
            "parsed_records": 0,
            "modification_masses_assessed": False,
        },
    }
    if backend == "pyxlms":
        _read_with_pyxlms(bundle)
    verify_maxlynx_sources(bundle)
    return bundle


def write_maxlynx_claudio_input(
    bundle, reference_mapping, output_dir, *,
    allowed_crosslinkers, allowed_residues,
    export_unverified_candidates=False,
):
    """Export explicit candidates without assigning per-record FDR values."""
    if export_unverified_candidates is not True:
        raise ValueError("Explicit unverified candidate export mode is required")
    if Path(output_dir).exists():
        raise FileExistsError("Output directory already exists")
    verify_maxlynx_sources(bundle)
    selections = {}
    for name, values in (
        ("crosslinkers", allowed_crosslinkers),
        ("residues", allowed_residues),
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
    native = read_maxlynx_csms(
        provenance["psm_path"], provenance["fasta_path"],
        dataset_accession=provenance["dataset_accession"],
        parameters_path=provenance["parameters_path"], backend="native",
    )
    for key in ("source_columns", "proteins", "native_settings", "fdr_policy"):
        if native[key] != bundle[key]:
            raise ValueError("MaxLynx bundle differs from its bound native sources")
    if len(native["records"]) != len(bundle["records"]):
        raise ValueError("MaxLynx source record count changed")
    for original, checked in zip(bundle["records"], native["records"]):
        for key in (
            "source_record_id", "source_record_number", "source_values",
            "source_fields", "record_kind", "decoy_status", "fdr",
            "source_pep", "endpoints", "localization",
        ):
            if original[key] != checked[key]:
                raise ValueError("MaxLynx source record changed after import")
        if original["record_status"] == "valid" and checked["record_status"] != "valid":
            raise ValueError("MaxLynx source validation changed")

    profiles, options, mappings = _merox_multi_reference_profiles(
        bundle, reference_mapping
    )
    crosslinker = bundle["native_settings"]["crosslinker"]
    rows, audit = [], []
    for record in bundle["records"]:
        entry = deepcopy(record)
        reasons = []
        if record["record_status"] != "valid":
            reasons.append(record["record_reason"] or "source_record_unresolved")
        if record["decoy_status"] != "target_target":
            reasons.append(
                "decoy_record" if record["decoy_status"] == "source_decoy"
                else "source_decoy_status_unresolved"
            )
        if crosslinker not in selections["crosslinkers"]:
            reasons.append("crosslinker_not_selected")
        projected = {}
        if not reasons:
            for site in ("a", "b"):
                endpoint = record["endpoints"][site]
                identifier = endpoint["source_fasta_id"]
                profile = profiles.get(identifier)
                if profile is None:
                    reasons.append("reference_mapping_unavailable_" + site)
                    continue
                result = _project_peptide(endpoint, profile, options[identifier])
                projected[site] = result
                if result["status"] != "projected":
                    reasons.append("reference_projection_unresolved_" + site)
                    continue
                if result["compatible_alignment_count"] != len(options[identifier]):
                    reasons.append("reference_projection_not_invariant_" + site)
                    continue
                interval = profile.get("reference_interval")
                if interval is not None and not (
                    interval[0] <= result["reference_peptide_start"]
                    and result["reference_peptide_end"] <= interval[1]
                ):
                    reasons.append("reference_peptide_outside_reviewed_core_" + site)
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

        fields = record["source_fields"]
        scope = (
            "exact_sequence"
            if all(
                profiles[record["endpoints"][site]["source_fasta_id"]][
                    "identity_scope"
                ] == "exact_sequence"
                for site in ("a", "b")
            )
            else "candidate_hypothesis"
        )
        row = {
            "source_record_id": record["source_record_id"],
            "dataset_accession": provenance["dataset_accession"],
            "source_record_number": record["source_record_number"],
            "source_psm_sha256": provenance["psm_sha256"],
            "source_fasta_sha256": provenance["fasta_sha256"],
            "source_parameters_sha256": provenance["parameters_sha256"],
            "source_crosslinker": crosslinker,
            "source_decoy_status": record["decoy_status"],
            "source_native_decoy_label": fields["Decoy"],
            "source_native_product_type": fields.get(
                "Crosslink product type",
                fields.get("Crosslink Product Type", ""),
            ),
            "source_pep_raw": fields["PEP"],
            "source_score_raw": fields["Score"],
            "source_spectrum_file": fields["Raw file"],
            "source_scan_number": fields["Scan number"],
            "input_scope": "unverified_candidate",
            "identity_scope": scope,
            "fdr_verification_status": bundle["fdr_policy"]["verification_status"],
            "native_configured_fdr_level": "psm",
            "native_configured_fdr_cutoff": bundle["fdr_policy"]["configured_cutoff"],
            "selected_fdr_level": None,
            "selected_fdr_value": None,
            "selected_fdr_cutoff": None,
            "source_fields_json": json.dumps(
                fields, ensure_ascii=False, sort_keys=True
            ),
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
                "source_modified_sequence_" + site: fields[
                    "Modified sequence1" if site == "a" else "Modified sequence2"
                ],
            })
        rows.append(row)

    bodies = {}
    for kind in ("psm", "fasta", "parameters"):
        path = Path(provenance[kind + "_path"])
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != provenance[kind + "_sha256"]:
            raise ValueError("Native source checksum changed before export")
        name = "native_sources/" + provenance[kind + "_sha256"] + path.suffix
        bodies[name] = raw
    bodies["native_settings.json"] = (
        json.dumps(bundle["native_settings"], indent=2) + "\n"
    ).encode("utf-8")
    prepared = {
        "rows": rows,
        "audit": audit,
        "fdr_policy": deepcopy(bundle["fdr_policy"]),
        "extra_source_bodies": bodies,
        "reference_projection": deepcopy(reference_mapping),
        "reference_projection_policy": {
            "kind": "source_bound_per_protein_reference_evidence",
            "require_all_optimal_alignments": True,
            "source_fasta_sha256": provenance["fasta_sha256"],
        },
    }
    return write_xifdr_outputs(
        bundle, mappings, prepared, output_dir,
        pipeline_policy={
            "input_scope": "unverified_candidate",
            "crosslinkers": sorted(selections["crosslinkers"]),
            "residues": sorted(selections["residues"]),
            "localization": "single_source_candidate",
            "coordinate_basis": 1,
        },
    )
