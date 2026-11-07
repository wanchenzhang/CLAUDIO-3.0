"""Import XL-MS PSMs, verify source coordinates, and prepare CLAUDIO inputs.

Coordinates are 1-based. Identity resolution requires full sequence equality.
FDR levels and source classifications are preserved rather than inferred.
"""

import configparser
import csv
import hashlib
import io
import json
import math
import re
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
from numbers import Real
from pathlib import Path
from urllib.parse import urlsplit
import zipfile
import requests


UNIPROT_REST = "https://rest.uniprot.org"
ACCESSION = re.compile(r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[1-9][0-9]*)?")
MODIFICATIONS = frozenset({"cm", "ox", "bs3", "bs3oh", "bs3nh2"})
FDR_COLUMNS = {"psm": "fdr", "peptide_pair": "PeptidePairFDR", "protein_a": "Protein1FDR",
               "protein_b": "Protein2FDR", "link": "LinkFDR", "protein_pair": "PPIFDR"}
REQUIRED_COLUMNS = {
    "PSMID", "Protein1", "Protein2", "Decoy1", "Decoy2", "PepSeq1", "PepSeq2",
    "PepPos1", "PepPos2", "PeptideLength1", "PeptideLength2", "LinkPos1", "LinkPos2",
    "ProteinLinkPos1", "ProteinLinkPos2", "isDecoy", "isTT", "isTD", "isDD", "fdr",
}


def _write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _fasta_records(text):
    records, header, parts = {}, None, []

    def finish():
        if header is None:
            return
        identifier = header.split()[0]
        if identifier in records:
            raise ValueError(f"Duplicate FASTA identifier: {identifier}")
        sequence = "".join(parts).upper()
        terminal_stop = sequence.endswith("*")
        if terminal_stop:
            sequence = sequence[:-1]
        if not sequence or re.fullmatch(r"[A-Z]+", sequence) is None:
            raise ValueError(f"Invalid FASTA sequence: {identifier}")
        fields = identifier.split("|")
        accession = fields[1] if len(fields) == 3 and fields[0] in ("sp", "tr") else identifier
        taxon = re.search(r"(?:^|\s)OX=(\d+)(?:\s|$)", header)
        records[identifier] = {
            "id": identifier, "header": header, "sequence": sequence,
            "uniprot_accession": accession if ACCESSION.fullmatch(accession) else None,
            "taxon_id": int(taxon[1]) if taxon else None,
            "terminal_stop_removed": terminal_stop,
        }

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish()
            header, parts = line[1:].strip(), []
            if not header:
                raise ValueError("Empty FASTA identifier")
        elif header is None:
            raise ValueError("Sequence appears before a FASTA header")
        else:
            parts.append("".join(line.split()))
    finish()
    if not records:
        raise ValueError("FASTA file has no sequences")
    return records


def read_search_fasta(path):
    """Read exact source identifiers, retaining UniProt accessions where explicit."""
    path = Path(path)
    raw = path.read_bytes()
    records = _fasta_records(raw.decode("utf-8-sig"))
    for record in records.values():
        record.update(fasta_path=str(path.resolve()), fasta_sha256=hashlib.sha256(raw).hexdigest())
    return records


def normalize_xifdr_peptide(raw, *, allowed_modifications=None):
    """Read supported xi inline modifications without changing residue positions."""
    allowed = MODIFICATIONS if allowed_modifications is None else frozenset(allowed_modifications)
    if not isinstance(raw, str) or not raw:
        raise ValueError("Empty peptide")
    sequence, modifications, position = [], [], 0
    while position < len(raw):
        residue = raw[position]
        if not "A" <= residue <= "Z":
            raise ValueError(f"Unsupported peptide annotation: {raw}")
        sequence.append(residue)
        position += 1
        match = re.match(r"[a-z][a-z0-9]*(?:-[a-z][a-z0-9]*)*", raw[position:])
        if match:
            token = match[0]
            if token not in allowed:
                raise ValueError(f"Unsupported xi modification: {token}")
            modifications.append({"position": len(sequence), "residue": residue, "token": token})
            position += len(token)
    return "".join(sequence), modifications


def _positive_integer(raw):
    if not isinstance(raw, str) or re.fullmatch(r"[0-9]+", raw.strip()) is None or int(raw) <= 0:
        raise ValueError("Expected a positive integer")
    return int(raw)


def _boolean(raw):
    return {"true": True, "false": False, "1": True, "0": False}.get(str(raw).strip().lower())


def _decoy_status(fields):
    values = [_boolean(fields.get(key)) for key in ("Decoy1", "Decoy2", "isDecoy", "isTT", "isTD", "isDD")]
    if any(value is None for value in values):
        return "unknown"
    a, b, decoy, tt, td, dd = values
    if (decoy, tt, td, dd) != (a or b, not a and not b, a != b, a and b):
        return "unknown"
    return "decoy_decoy" if dd else "target_decoy" if td else "target_target"


def _fdr_value(raw):
    result = {"raw": raw, "value": None, "status": "missing"}
    if raw is None or not str(raw).strip():
        return result
    try:
        value = float(raw)
    except (ValueError, TypeError):
        result["status"] = "invalid"
    else:
        if not math.isfinite(value):
            result["status"] = "non_finite"
        elif not 0 <= value <= 1:
            result["status"] = "out_of_range"
        else:
            result.update(value=value, status="valid")
    return result


def _endpoint(fields, site, proteins, allowed_modifications):
    index = "1" if site == "a" else "2"
    protein_id = fields.get("Protein" + index, "")
    result = {"source_protein_id": protein_id, "source_peptide": fields.get("PepSeq" + index),
              "source_fasta_id": None, "peptide": None, "modifications": [],
              "peptide_start": None, "link_position": None, "protein_position": None,
              "status": "not_evaluable", "reason": ""}
    try:
        result.update(peptide_start=_positive_integer(fields.get("PepPos" + index)),
                      link_position=_positive_integer(fields.get("LinkPos" + index)),
                      protein_position=_positive_integer(fields.get("ProteinLinkPos" + index)))
        length = _positive_integer(fields.get("PeptideLength" + index))
    except ValueError:
        result["reason"] = "invalid_position"
        return result
    try:
        peptide, modifications = normalize_xifdr_peptide(result["source_peptide"], allowed_modifications=allowed_modifications)
    except ValueError:
        result["reason"] = "unsupported_peptide_annotation"
        return result
    result.update(peptide=peptide, modifications=modifications)
    if len(peptide) != length:
        result["reason"] = "peptide_length_mismatch"
        return result
    matches = ([proteins[protein_id]] if protein_id in proteins else
               [p for p in proteins.values() if p["uniprot_accession"] == protein_id])
    if len(matches) != 1:
        result["reason"] = "protein_identity_ambiguous" if len(matches) > 1 or re.search(r"[;,]", protein_id) else "protein_not_in_search_fasta"
        return result
    protein = matches[0]
    result["source_fasta_id"] = protein["id"]
    start, link, position = result["peptide_start"], result["link_position"], result["protein_position"]
    sequence = protein["sequence"]
    if position > len(sequence):
        result["reason"] = "protein_position_out_of_range"
    elif link > len(peptide):
        result["reason"] = "peptide_link_position_out_of_range"
    elif position != start + link - 1:
        result["reason"] = "coordinate_inconsistent"
    elif sequence[start - 1:start - 1 + len(peptide)] != peptide:
        result["reason"] = "peptide_sequence_mismatch"
    else:
        result.update(status="valid", reason="")
    return result


def read_xifdr_psms(psm_path, fasta_path, *, dataset_accession,
                    allowed_modifications=None, backend="native"):
    """Keep all nonblank PSM records, including decoys and unresolved source rows.

    Missing or duplicated required column labels reject the file. Ragged data
    records stay in the audit without assigning possibly shifted field values.
    The pyxlms backend checks parsed endpoints against the source FASTA and
    restores the original endpoint order. Modification masses are not assessed.
    """
    if backend not in ("native", "pyxlms"):
        raise ValueError("backend must be 'native' or 'pyxlms'")
    if not isinstance(dataset_accession, str) or re.fullmatch(r"PXD\d{6,}", dataset_accession) is None:
        raise ValueError("A ProteomeXchange PXD accession is required")
    psm_path, fasta_path = Path(psm_path), Path(fasta_path)
    raw = psm_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    proteins = read_search_fasta(fasta_path)
    reader = csv.reader(io.StringIO(raw.decode("utf-8-sig")))
    columns = next(reader, [])
    named = [column for column in columns if column]
    if len(named) != len(set(named)):
        raise ValueError("xiFDR PSM table has duplicate named columns")
    missing = REQUIRED_COLUMNS - set(columns)
    if missing:
        raise ValueError("Not a supported xiFDR PSM table; missing columns: " + ", ".join(sorted(missing)))
    records, previous_line = [], reader.line_num
    for number, values in enumerate(reader, 1):
        line_start, line_end = previous_line + 1, reader.line_num
        previous_line = line_end
        if not any(value.strip() for value in values):
            continue
        fields = {key: value for key, value in zip(columns, values) if key} if len(values) == len(columns) else {}
        record = {"source_record_id": f"{dataset_accession}:{digest}:{number:06d}",
                  "source_record_number": number, "source_line_start": line_start, "source_line_end": line_end,
                  "source_values": values, "source_fields": fields,
                  "decoy_status": _decoy_status(fields),
                  "fdr": {key: _fdr_value(fields.get(column)) for key, column in FDR_COLUMNS.items()},
                  "endpoints": {}, "record_status": "not_evaluable", "record_reason": ""}
        if len(values) != len(columns):
            record["record_reason"] = "row_width_mismatch"
        else:
            record["endpoints"] = {site: _endpoint(fields, site, proteins, allowed_modifications) for site in ("a", "b")}
            if record["decoy_status"] == "unknown":
                record["record_reason"] = "invalid_decoy_flags"
            elif any(e["status"] != "valid" for e in record["endpoints"].values()):
                record["record_reason"] = "source_endpoint_unresolved"
            else:
                record["record_status"] = "valid"
        records.append(record)
    bundle = {"format": "xifdr_psm", "format_version": 1, "source_columns": columns,
              "records": records, "proteins": proteins,
              "provenance": {"dataset_accession": dataset_accession, "psm_path": str(psm_path.resolve()),
                             "psm_sha256": digest, "fasta_path": str(fasta_path.resolve()),
                             "fasta_sha256": hashlib.sha256(fasta_path.read_bytes()).hexdigest()}}
    if backend == "pyxlms":
        _read_with_pyxlms(bundle)
    return bundle


def _read_with_pyxlms(bundle):
    """Read validated rows with pyXLMS, retaining the original source ledger."""
    try:
        from importlib.metadata import version
        from pyXLMS import parser
    except ImportError as error:
        raise RuntimeError("The pyxlms backend requires pyXLMS==2.0.7") from error

    bundle["parser"] = {"name": "pyXLMS", "version": version("pyXLMS"),
                        "modification_masses_assessed": False}
    records = [record for record in bundle["records"] if record["record_status"] == "valid"]
    if not records:
        return

    # The marker exists only in the in-memory copy; original fields stay intact.
    marker = "__claudio_source_record_id"
    while marker in bundle["source_columns"]:
        marker += "_"
    stream = io.StringIO()
    writer = csv.writer(stream)
    writer.writerow(bundle["source_columns"] + [marker])
    for record in records:
        writer.writerow(record["source_values"] + [record["source_record_id"]])
    stream.seek(0)

    def reject(record, reason, message):
        record.update(record_status="not_evaluable", record_reason=reason,
                      parser_error=message)

    try:
        result = parser.read_xi(stream, decoy_prefix="", parse_modifications=False,
                                dtype=str, keep_default_na=False)
        parsed = result["crosslink-spectrum-matches"] or []
        by_source = {}
        expected_ids = {record["source_record_id"] for record in records}
        for csm in parsed:
            source = (csm["additional_information"] or {}).get("source", {})
            source_id = source.get(marker)
            if source_id not in expected_ids or source_id in by_source:
                raise ValueError("pyXLMS returned an unknown or repeated source record")
            by_source[source_id] = csm
    except (ValueError, TypeError, KeyError, RuntimeError, IndexError, AttributeError, OverflowError) as error:
        for record in records:
            reject(record, "pyxlms_parse_failed", f"{type(error).__name__}: {error}")
        return

    def original_endpoint(endpoint, decoy):
        return (endpoint["peptide"], endpoint["link_position"],
                (endpoint["source_protein_id"],), (endpoint["protein_position"],),
                (endpoint["peptide_start"],), decoy)

    def parsed_endpoint(csm, prefix):
        return (csm[prefix + "_peptide"], csm[prefix + "_peptide_crosslink_position"],
                tuple(csm[prefix + "_proteins"] or []),
                tuple(csm[prefix + "_proteins_crosslink_positions"] or []),
                tuple(csm[prefix + "_proteins_peptide_positions"] or []),
                csm[prefix + "_decoy"])

    for record in records:
        csm = by_source.get(record["source_record_id"])
        if csm is None:
            reject(record, "pyxlms_parse_failed", "pyXLMS did not return this source record")
            continue
        fields = record["source_fields"]
        expected = tuple(original_endpoint(record["endpoints"][site], _boolean(fields["Decoy" + index]))
                         for site, index in (("a", "1"), ("b", "2")))
        try:
            actual = tuple(parsed_endpoint(csm, prefix) for prefix in ("alpha", "beta"))
            direct, swapped = expected == actual, expected == actual[::-1]
            if not (direct or swapped):
                reject(record, "pyxlms_endpoint_mismatch", "Parsed endpoints differ from the validated source")
                continue
            record["parser_orientation"] = "indistinguishable" if direct and swapped else "same" if direct else "swapped"
            prefixes = ("alpha", "beta") if direct else ("beta", "alpha")
            updates = {}
            for site, prefix in zip(("a", "b"), prefixes):
                endpoint = record["endpoints"][site]
                tokens = parser.parse_modifications_from_xi_sequence(endpoint["source_peptide"])
                if tokens != {mod["position"]: mod["token"] for mod in endpoint["modifications"]}:
                    raise ValueError("Parsed modification tokens differ from the source")
                peptide = csm[prefix + "_peptide"]
                updates[site] = {
                    "peptide": peptide,
                    "link_position": csm[prefix + "_peptide_crosslink_position"],
                    "protein_position": csm[prefix + "_proteins_crosslink_positions"][0],
                    "peptide_start": csm[prefix + "_proteins_peptide_positions"][0],
                    "modifications": [{"position": position, "residue": peptide[position - 1], "token": token}
                                      for position, token in tokens.items()],
                }
            for site, update in updates.items():
                record["endpoints"][site].update(update)
        except (ValueError, TypeError, KeyError, RuntimeError, IndexError, AttributeError, OverflowError) as error:
            reject(record, "pyxlms_parse_failed", f"{type(error).__name__}: {error}")


def _plink_endpoint(site, proteins, modification_status):
    """Validate one reported alternative with the existing FASTA checks."""
    if len(site["proteins"]) != 1 or len(site["protein_positions"]) != 1:
        raise ValueError("A pLink protein alternative must contain one site per endpoint")
    peptide, link = site["peptide"], site["peptide_position"]
    position = site["protein_positions"][0]
    fields = {"Protein1": site["proteins"][0], "PepSeq1": peptide,
              "PepPos1": str(position - link + 1), "PeptideLength1": str(len(peptide)),
              "LinkPos1": str(link), "ProteinLinkPos1": str(position)}
    endpoint = _endpoint(fields, "a", proteins, ())
    endpoint["modifications"] = [] if modification_status == "none_reported" else None
    return endpoint


PLINK_JOIN_FIELDS = ("Title", "Charge", "Peptide", "Proteins", "FileID", "Modifications")


def _plink_source_fields(columns, values):
    """Recognize pLink's unquoted File/NativeID title without merging arbitrary cells."""
    if len(values) == len(columns):
        return {key: value for key, value in zip(columns, values) if key}, None
    standard = (
        "Order", "Title", "Charge", "Precursor_Mass", "Peptide", "Peptide_Type", "Linker",
        "Peptide_Mass", "Modifications", "Evalue", "Score", "Precursor_Mass_Error(Da)",
        "Precursor_Mass_Error(ppm)", "Proteins", "Protein_Type", "FileID", "LabelID",
        "Alpha_Matched", "Beta_Matched", "Alpha_Evalue", "Beta_Evalue",
    )
    if len(values) <= len(columns) or tuple(columns[:len(standard)]) != standard:
        return {}, None
    title_end = 2 + len(values) - len(columns)
    title = ",".join(values[1:title_end])
    if re.fullmatch(r'[^,\r\n]+_File:"[^"\r\n]+",_NativeID:"[^"\r\n]+"', title) is None:
        return {}, None
    aligned = [values[0], title] + values[title_end:]
    fields = {key: value for key, value in zip(columns, aligned) if key}
    if (fields["Peptide_Type"].strip().lower() != "cross-linked"
            or fields["Protein_Type"].strip().lower() not in ("intra-protein", "inter-protein")
            or not fields["Linker"].strip() or not fields["Proteins"].strip()
            or re.fullmatch(r"[A-Z]+\([1-9][0-9]*\)-[A-Z]+\([1-9][0-9]*\)", fields["Peptide"]) is None):
        return {}, None
    try:
        for name in ("Order", "Charge", "FileID", "LabelID"):
            _positive_integer(fields[name])
        for name in ("Precursor_Mass", "Peptide_Mass", "Evalue", "Score",
                     "Precursor_Mass_Error(Da)", "Precursor_Mass_Error(ppm)",
                     "Alpha_Matched", "Beta_Matched", "Alpha_Evalue", "Beta_Evalue"):
            value = float(fields[name])
            signed = name in ("Precursor_Mass_Error(Da)", "Precursor_Mass_Error(ppm)")
            if not math.isfinite(value) or (not signed and value < 0):
                return {}, None
    except ValueError:
        return {}, None
    return fields, {
        "method": "plink_unquoted_title_metadata", "original_field_count": len(values),
        "interpreted_field_count": len(columns),
        "title_source_field_indices_zero_based": list(range(1, title_end)),
    }


def _plink_filter_policy(path):
    """Read the author's settings without treating them as independent validation."""
    path = Path(path)
    raw = path.read_bytes()
    settings = configparser.ConfigParser(interpolation=None)
    try:
        settings.read_string(raw.decode("utf-8-sig"))
        software_version = settings.get("version", "version")
        if not re.fullmatch(r"2(?:\.[0-9]+)+", software_version):
            raise ValueError("Native evidence currently supports pLink 2 parameter files")
        cutoff = _fdr_value(settings.get("filter", "FDR"))
        if cutoff["status"] != "valid":
            raise ValueError("The pLink parameter FDR must be a finite fraction between 0 and 1")
        count = _positive_integer(settings.get("linker", "linker_total"))
        linkers = [settings.get("linker", f"linker{number}") for number in range(1, count + 1)]
        grouping = settings.get("filter", "fdr_separate_or_global")
        database = settings.get("database", "db_name")
        if not grouping.strip() or not database.strip() or any(not linker.strip() for linker in linkers):
            raise ValueError("The pLink parameters have empty database, linker or FDR grouping settings")
    except configparser.Error as error:
        raise ValueError(f"Invalid pLink parameters: {error}") from error
    return {
        "level": "psm", "cutoff": cutoff["value"], "cutoff_raw": cutoff["raw"],
        "grouping_raw": grouping, "software_version": software_version,
        "database_name": database, "crosslinkers": linkers,
        "filter_tolerance_raw": settings.get("filter", "filter_tol", fallback=None),
        "filter_tolerance_unit_raw": settings.get("filter", "filter_tol_type", fallback=None),
        "evalue_setting_raw": settings.get("filter", "evalue", fallback=None),
        "independent_validation_performed": False,
    }, hashlib.sha256(raw).hexdigest()


def _add_plink_native_evidence(bundle, unfiltered_path, parameters_path):
    """Join native rows by source fields; never transfer confidence by row order."""
    unfiltered_path, parameters_path = Path(unfiltered_path), Path(parameters_path)
    policy, parameters_hash = _plink_filter_policy(parameters_path)
    raw = unfiltered_path.read_bytes()
    reader = csv.reader(io.StringIO(raw.decode("utf-8-sig")))
    columns = next(reader, [])
    named = [name for name in columns if name]
    required = set(PLINK_JOIN_FIELDS) | {"Peptide_Type", "Target_Decoy", "Q-value", "Protein_Type"}
    if len(named) != len(set(named)):
        raise ValueError("The native pLink unfiltered table has duplicate named columns")
    if not required <= set(columns):
        raise ValueError("The native pLink unfiltered table is missing columns: "
                         + ", ".join(sorted(required - set(columns))))

    wanted = {tuple(record["source_fields"][name] for name in PLINK_JOIN_FIELDS)
              for record in bundle["records"] if set(PLINK_JOIN_FIELDS) <= record["source_fields"].keys()}
    index, previous_line, ragged = {}, reader.line_num, []
    for number, values in enumerate(reader, 1):
        start, end = previous_line + 1, reader.line_num
        previous_line = end
        if not any(value.strip() for value in values):
            continue
        if len(values) != len(columns):
            ragged.append(number)
            continue
        fields = {name: value for name, value in zip(columns, values) if name}
        key = tuple(fields[name] for name in PLINK_JOIN_FIELDS)
        if key not in wanted:
            continue
        entry = {"source_record_number": number, "source_line_start": start,
                 "source_line_end": end, "source_values": values, "source_fields": fields}
        index.setdefault(key, []).append(entry)

    bundle["native_filter_policy"] = policy
    bundle["provenance"].update(
        plink_unfiltered_path=str(unfiltered_path.resolve()),
        plink_unfiltered_sha256=hashlib.sha256(raw).hexdigest(),
        plink_parameters_path=str(parameters_path.resolve()),
        plink_parameters_sha256=parameters_hash,
    )
    bundle["parser"]["native_evidence"] = {
        "join_fields": list(PLINK_JOIN_FIELDS), "q_value_level": "psm",
        "unfiltered_ragged_record_numbers": ragged,
    }
    decoy_codes = {"0": "decoy_decoy", "1": "target_decoy", "2": "target_target"}
    for record in bundle["records"]:
        fields = record["source_fields"]
        evidence = {"status": "not_evaluable", "reason": "", "candidate_record_numbers": [],
                    "within_reported_fdr_cutoff": None}
        record["native_evidence"] = evidence
        if (not set(PLINK_JOIN_FIELDS) <= fields.keys()
                or any(not fields[name].strip() for name in PLINK_JOIN_FIELDS[:-1])):
            evidence["reason"] = "missing_source_join_fields"
            continue
        if fields.get("Peptide_Type", "").strip().lower() != "cross-linked":
            evidence["reason"] = "unsupported_source_peptide_type"
            continue
        if fields.get("Linker") not in policy["crosslinkers"]:
            evidence["reason"] = "crosslinker_conflict"
            continue
        matches = index.get(tuple(fields[name] for name in PLINK_JOIN_FIELDS), [])
        evidence["candidate_record_numbers"] = [entry["source_record_number"] for entry in matches]
        if len(matches) != 1:
            evidence["reason"] = "ambiguous_native_psm" if matches else "native_psm_not_found"
            continue
        native = matches[0]
        reported = native["source_fields"]
        evidence.update(deepcopy(native))
        evidence["native_order"] = reported.get("Order")
        evidence["q_value"] = _fdr_value(reported["Q-value"])
        if reported["Peptide_Type"].strip() != "3":
            evidence["reason"] = "native_record_not_crosslinked"
            continue
        source_type = fields.get("Protein_Type", "").strip().lower()
        expected_type = {"intra-protein": "1", "inter-protein": "2"}.get(source_type)
        if expected_type is not None and reported["Protein_Type"].strip() != expected_type:
            evidence["reason"] = "native_protein_type_conflict"
            continue
        record["decoy_status"] = decoy_codes.get(reported["Target_Decoy"].strip(), "unknown")
        record["fdr"]["psm"] = deepcopy(evidence["q_value"])
        if record["decoy_status"] == "unknown":
            evidence["reason"] = "unknown_native_decoy_code"
            record["decoy_reason"] = evidence["reason"]
        elif evidence["q_value"]["status"] != "valid":
            evidence["reason"] = "invalid_native_q_value"
            record["decoy_reason"] = ""
        elif reported["Protein_Type"].strip() not in ("1", "2"):
            evidence["reason"] = "unknown_native_protein_type"
            record["decoy_reason"] = ""
        else:
            evidence.update(status="verified", reason="",
                            within_reported_fdr_cutoff=(record["decoy_status"] == "target_target"
                                                       and evidence["q_value"]["value"] <= policy["cutoff"]))
            record["decoy_reason"] = ""
        if record["record_reason"] == "decoy_policy_not_provided":
            if evidence["status"] == "verified":
                record.update(record_status="valid", record_reason="")
            else:
                record["record_reason"] = evidence["reason"]
    return bundle


def _verify_plink_native_sources(bundle):
    if (
        bundle.get("format") != "plink_psm"
        or "native_filter_policy" not in bundle
    ):
        return
    for kind in ("plink_unfiltered", "plink_parameters"):
        source = Path(bundle["provenance"][kind + "_path"])
        digest = hashlib.sha256()
        with source.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != bundle["provenance"][kind + "_sha256"]:
            raise ValueError(f"A native pLink source changed after parsing: {source}")


def _author_filter_source_bodies(bundle):
    """Verify the reviewed declaration and return portable copies of its evidence."""
    if "author_filter_policy" not in bundle:
        return {}
    policy = bundle["author_filter_policy"]
    path = Path(bundle["provenance"]["author_filter_path"])
    raw = path.read_bytes()
    if (hashlib.sha256(raw).hexdigest() != bundle["provenance"]["author_filter_sha256"]
            or json.loads(raw.decode("utf-8-sig")) != policy):
        raise ValueError("The author filter source declaration changed after parsing")
    for kind in ("psm", "fasta"):
        digest = hashlib.sha256()
        with Path(bundle["provenance"][kind + "_path"]).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != policy[kind + "_sha256"]:
            raise ValueError("An input source changed after author filter review")
    sources = policy.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("Author filter evidence requires supporting sources")
    bodies = {"author_filter_sources/" + path.name: raw}
    root = path.parent.resolve()
    for source in sources:
        if not isinstance(source, dict) or not isinstance(source.get("path"), str) or not source["path"]:
            raise ValueError("Invalid author filter source path")
        relative = Path(source["path"])
        resolved = (root / relative).resolve()
        if relative.is_absolute() or ".." in relative.parts or root not in resolved.parents:
            raise ValueError("Author filter source paths must stay within the evidence directory")
        name = "author_filter_sources/" + relative.as_posix()
        if name in bodies or not isinstance(source.get("role"), str) or not source["role"].strip():
            raise ValueError("Author filter sources must have unique paths and explicit roles")
        body = resolved.read_bytes()
        if hashlib.sha256(body).hexdigest() != source.get("sha256"):
            raise ValueError("An author filter source checksum does not match")
        bodies[name] = body
    return bodies


def _write_source_bodies(output_dir, bodies):
    for name, raw in bodies.items():
        target = output_dir / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)


def _add_author_filter_evidence(bundle, path):
    """Retain an author's table-level acceptance without assigning per-PSM FDR."""
    path = Path(path)
    raw = path.read_bytes()
    policy = json.loads(raw.decode("utf-8-sig"))
    json.dumps(policy, allow_nan=False)
    if (not isinstance(policy, dict) or policy.get("schema_version") != 1
            or policy.get("kind") != "xlms_author_filter_evidence"
            or policy.get("review_status") != "reviewed"
            or policy.get("result_scope") != "author_filtered_psm_table"
            or policy.get("target_only") is not True
            or policy.get("fasta_role") != "target_search_database"
            or policy.get("filter_level") != "psm"
            or policy.get("per_record_fdr_available") is not False):
        raise ValueError("Unsupported or unreviewed author filter declaration")
    _fdr_policy(policy["filter_level"], policy.get("reported_fdr_cutoff"))
    provenance = bundle["provenance"]
    if any(policy.get(key) != provenance[key] for key in ("dataset_accession", "psm_sha256", "fasta_sha256")):
        raise ValueError("The author filter dataset or source checksum does not match this import")
    bundle["author_filter_policy"] = policy
    provenance.update(author_filter_path=str(path.resolve()), author_filter_sha256=hashlib.sha256(raw).hexdigest())
    _author_filter_source_bodies(bundle)
    for record in bundle["records"]:
        pairs = record["protein_pairs"]
        target_member = bool(pairs) and all(
            pair[side]["source_protein_id"] in bundle["proteins"]
            for pair in pairs for side in ("a", "b")
        )
        record["source_filter_status"] = "reported_accepted" if target_member else "not_evaluable"
        record["source_filter_reason"] = "" if target_member else "target_database_membership_unresolved"
        if target_member:
            record.update(decoy_status="target_target", decoy_reason="",
                          decoy_evidence="author_filtered_table_and_declared_target_database")
            if record["record_reason"] == "decoy_policy_not_provided":
                record.update(record_status="valid", record_reason="")
    return bundle

def _read_plink_workbook(raw):
    """Read one pLink worksheet without evaluating or repairing source cells."""
    import posixpath
    import xml.etree.ElementTree as ET
    from zipfile import ZipFile
    from openpyxl import load_workbook
    from openpyxl.utils import get_column_letter

    required = {"Title", "Peptide", "Proteins", "Peptide_Type", "Linker"}
    with ZipFile(io.BytesIO(raw)) as archive:
        names = archive.namelist()
        required_parts = (
            "[Content_Types].xml",
            "xl/workbook.xml",
            "xl/_rels/workbook.xml.rels",
        )
        if any(names.count(name) != 1 for name in required_parts):
            raise ValueError("Not a supported Excel workbook")

        workbook_xml = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = {
            node.get("Id"): node
            for node in ET.fromstring(
                archive.read("xl/_rels/workbook.xml.rels")
            )
        }
        workbook = load_workbook(
            io.BytesIO(raw),
            read_only=True,
            data_only=False,
            keep_links=False,
        )
        try:
            matches = []
            for sheet in workbook.worksheets:
                for row_number, cells in enumerate(sheet.iter_rows(), 1):
                    values = [cell.value for cell in cells]
                    if not any(
                        value.strip() if isinstance(value, str)
                        else value is not None
                        for value in values
                    ):
                        continue
                    while values and values[-1] in (None, ""):
                        values.pop()
                    if required <= {
                        value for value in values if isinstance(value, str)
                    }:
                        if any(
                            value is not None and not isinstance(value, str)
                            for value in values
                        ):
                            raise ValueError("Non-text pLink worksheet headers")
                        matches.append((
                            sheet,
                            row_number,
                            ["" if value is None else value for value in values],
                        ))
                    break

            if len(matches) != 1:
                raise ValueError(
                    "Expected exactly one pLink PSM worksheet; "
                    f"found {len(matches)}"
                )
            sheet, header_row, columns = matches[0]
            named = [column for column in columns if column]
            if len(named) != len(set(named)):
                raise ValueError("pLink worksheet has duplicate named columns")

            sheet_nodes = workbook_xml.findall("./{*}sheets/{*}sheet")
            selected = [
                (index, node)
                for index, node in enumerate(sheet_nodes, 1)
                if node.get("name") == sheet.title
            ]
            if len(selected) != 1:
                raise ValueError("Ambiguous worksheet package identity")
            sheet_index, node = selected[0]
            relation_id = node.get(
                "{http://schemas.openxmlformats.org/officeDocument/"
                "2006/relationships}id"
            )
            relation = relationships.get(relation_id)
            if relation is None or relation.get("TargetMode") == "External":
                raise ValueError("Unsupported worksheet relationship")
            target = relation.get("Target", "")
            member = posixpath.normpath(
                target.lstrip("/") if target.startswith("/")
                else posixpath.join("xl", target)
            )
            if (
                not relation.get("Type", "").endswith("/worksheet")
                or not member.startswith("xl/worksheets/")
                or names.count(member) != 1
            ):
                raise ValueError("Invalid worksheet package location")

            raw_cells = {}
            with archive.open(member) as handle:
                for event, element in ET.iterparse(handle, events=("end",)):
                    if element.tag.rsplit("}", 1)[-1] != "c":
                        continue
                    coordinate = element.get("r")
                    if not coordinate or coordinate in raw_cells:
                        raise ValueError(
                            "Missing or duplicate worksheet cell address"
                        )
                    value = element.find("{*}v")
                    formula = element.find("{*}f")
                    raw_cells[coordinate] = {
                        "raw_xml_value": None if value is None else value.text,
                        "raw_xml_type": element.get("t"),
                        "raw_formula": None if formula is None else formula.text,
                    }
                    element.clear()

            rows = []
            for row_number, cells in enumerate(
                sheet.iter_rows(min_row=header_row + 1), header_row + 1
            ):
                if not any(
                    cell.value.strip() if isinstance(cell.value, str)
                    else cell.value is not None
                    for cell in cells
                ):
                    continue
                occupied = [
                    index for index, cell in enumerate(cells, 1)
                    if cell.value is not None and cell.value != ""
                ]
                width = max(len(columns), max(occupied, default=0))
                original, values, notes, issue = [], [], [], ""
                for index, cell in enumerate(cells[:width], 1):
                    coordinate = f"{get_column_letter(index)}{row_number}"
                    value = cell.value
                    data_type = cell.data_type
                    if data_type == "f":
                        issue = issue or "xlsx_formula_cell"
                    elif data_type == "e":
                        issue = issue or "xlsx_error_cell"
                    elif data_type == "d":
                        issue = issue or "xlsx_date_cell"
                    elif data_type not in ("s", "n", "inlineStr"):
                        issue = issue or "xlsx_unsupported_cell_type"

                    if isinstance(value, float) and not math.isfinite(value):
                        issue = issue or "xlsx_nonfinite_number"
                        value = repr(value)
                    elif value is not None and not isinstance(
                        value, (str, int, float, bool)
                    ):
                        value = (
                            value.isoformat()
                            if hasattr(value, "isoformat") else str(value)
                        )
                    original.append(value)
                    values.append("" if value is None else str(value))
                    notes.append({
                        "coordinate": coordinate,
                        "value": value,
                        "data_type": data_type,
                        "number_format": cell.number_format,
                        **raw_cells.get(coordinate, {
                            "raw_xml_value": None,
                            "raw_xml_type": None,
                            "raw_formula": None,
                        }),
                    })
                rows.append({
                    "source_record_number": len(rows) + 1,
                    "source_line_start": None,
                    "source_line_end": None,
                    "source_values": original,
                    "values": values,
                    "source_sheet": sheet.title,
                    "source_sheet_index": sheet_index,
                    "source_worksheet_row": row_number,
                    "source_cells": notes,
                    "source_cell_issue": issue,
                })
            return columns, rows
        finally:
            workbook.close()
def read_plink_psms(psm_path, fasta_path, *, dataset_accession,
                    decoy_prefix=None, fdr_columns=None, unfiltered_path=None, parameters_path=None,
                    author_filter_path=None):
    """Read pLink spectra tables while retaining unverified source information.

    Protein alternatives keep their original pairings. The caller must supply
    the decoy prefix and mappings from FDR levels to reported column names, or
    provide the paired native unfiltered results and pLink 2 parameters. A
    reviewed author_filter_path can instead document a target-only filtered
    PSM table while leaving all unreported per-record FDR values missing.
    Missing policies stay unresolved. Scores are not interpreted as FDR.
    Modification strings and spectrum titles are retained without inventing
    modification masses or numeric scan identifiers.
    """
    if not isinstance(dataset_accession, str) or re.fullmatch(r"PXD\d{6,}", dataset_accession) is None:
        raise ValueError("A ProteomeXchange PXD accession is required")
    if author_filter_path is not None and (decoy_prefix is not None or fdr_columns
                                           or unfiltered_path is not None or parameters_path is not None):
        raise ValueError("Do not mix author filter evidence with native or manual FDR policies")
    if (unfiltered_path is None) != (parameters_path is None):
        raise ValueError("Provide unfiltered_path and parameters_path together for native pLink evidence")
    if unfiltered_path is not None and (decoy_prefix is not None or fdr_columns):
        raise ValueError("Do not mix native pLink evidence with manual decoy or FDR column policies")
    if decoy_prefix is not None and (not isinstance(decoy_prefix, str)
                                    or not decoy_prefix or decoy_prefix != decoy_prefix.strip()):
        raise ValueError("decoy_prefix must be a nonempty prefix or None")
    fdr_columns = {} if fdr_columns is None else fdr_columns
    if (not isinstance(fdr_columns, dict)
            or any(level not in FDR_COLUMNS or not isinstance(column, str) or not column
                   for level, column in fdr_columns.items())):
        raise ValueError("FDR columns must map supported levels to reported column names")
    try:
        from importlib.metadata import version
        from pyXLMS import parser
    except ImportError as error:
        raise RuntimeError("Reading pLink PSMs requires pyXLMS==2.0.7") from error

    psm_path, fasta_path = Path(psm_path), Path(fasta_path)
    raw = psm_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    proteins = read_search_fasta(fasta_path)
    from zipfile import is_zipfile

    is_workbook = is_zipfile(io.BytesIO(raw))
    if is_workbook:
        columns, source_rows = _read_plink_workbook(raw)
    else:
        reader = csv.reader(io.StringIO(raw.decode("utf-8-sig")))
        columns = next(reader, [])

        def csv_source_rows():
            previous_line = reader.line_num
            for number, values in enumerate(reader, 1):
                line_start, line_end = previous_line + 1, reader.line_num
                previous_line = line_end
                yield {
                    "source_record_number": number,
                    "source_line_start": line_start,
                    "source_line_end": line_end,
                    "source_values": values,
                    "values": values,
                }

        source_rows = csv_source_rows()
    named = [column for column in columns if column]
    if len(named) != len(set(named)):
        raise ValueError("pLink PSM table has duplicate named columns")
    required = {"Title", "Peptide", "Proteins", "Peptide_Type", "Linker"}
    if not required <= set(columns):
        raise ValueError("Not a supported pLink PSM spectra table; missing columns: "
                         + ", ".join(sorted(required - set(columns))))
    if not set(fdr_columns.values()) <= set(columns):
        raise ValueError("Selected FDR columns are missing from the pLink PSM table")

    records = []
    for source_row in source_rows:
        number = source_row["source_record_number"]
        values = source_row["values"]
        if not any(value.strip() for value in values):
            continue
        if is_workbook:
            fields = {
                key: value for key, value in zip(columns, values) if key
            }
            interpretation = None
            source_key = (
                f"worksheet_{source_row['source_sheet_index']}:"
                f"{source_row['source_worksheet_row']:06d}"
            )
        else:
            fields, interpretation = _plink_source_fields(columns, values)
            source_key = f"{number:06d}"

        annotation = fields.get("Modifications", "").strip().lower()
        modification_status = (
            "none_reported" if annotation in ("null", "none", "nan")
            else "not_parsed" if annotation else "not_provided"
        )
        record = {
            "source_record_id": f"{dataset_accession}:{digest}:{source_key}",
            "source_record_number": number,
            "source_line_start": source_row["source_line_start"],
            "source_line_end": source_row["source_line_end"],
            "source_values": source_row["source_values"],
            "source_fields": fields,
            "source_scan": None,
            "decoy_status": "unknown",
            "decoy_reason": (
                "decoy_policy_not_provided" if decoy_prefix is None else ""
            ),
            "fdr": {
                level: _fdr_value(fields.get(fdr_columns.get(level)))
                for level in FDR_COLUMNS
            },
            "endpoints": {},
            "protein_pairs": [],
            "modification_status": modification_status,
            "record_status": "not_evaluable",
            "record_reason": "",
        }
        if is_workbook:
            record.update({
                key: source_row[key] for key in (
                    "source_sheet", "source_sheet_index",
                    "source_worksheet_row", "source_cells",
                    "source_cell_issue",
                )
            })
        records.append(record)
        if is_workbook and source_row["source_cell_issue"]:
            record["record_reason"] = source_row["source_cell_issue"]
            continue
        if interpretation is not None:
            record["source_csv_interpretation"] = interpretation
        if len(values) != len(columns) and interpretation is None:
            record["record_reason"] = "row_width_mismatch"
            continue
        if fields["Peptide_Type"].strip().lower() != "cross-linked":
            record["record_reason"] = "unsupported_peptide_type"
            continue
        peptide = fields["Peptide"].strip()
        expected = re.fullmatch(r"([A-Z]+)\(([0-9]+)\)-([A-Z]+)\(([0-9]+)\)", peptide)
        if expected is None:
            record["record_reason"] = "unsupported_peptide_annotation"
            continue
        try:
            # Parse each alternative separately to avoid introducing a Cartesian product.
            alternatives = fields["Proteins"].strip().rstrip("/").split("/")
            for alternative in alternatives:
                sites = parser.parse_sites_from_plink(peptide, alternative.strip())
                if ((sites[0]["peptide"], sites[0]["peptide_position"],
                     sites[1]["peptide"], sites[1]["peptide_position"])
                        != (expected[1], int(expected[2]), expected[3], int(expected[4]))):
                    raise ValueError("Parsed peptide annotations differ from the reported pLink sites")
                record["protein_pairs"].append({
                    "a": _plink_endpoint(sites[0], proteins, modification_status),
                    "b": _plink_endpoint(sites[1], proteins, modification_status),
                })
        except (ValueError, TypeError, KeyError, RuntimeError, IndexError, AttributeError, OverflowError) as error:
            record.update(record_reason="pyxlms_parse_failed", parser_error=f"{type(error).__name__}: {error}")
            continue

        for side in ("a", "b"):
            candidates = {}
            for pair in record["protein_pairs"]:
                endpoint = pair[side]
                key = (endpoint["source_protein_id"], endpoint["protein_position"],
                       endpoint["link_position"], endpoint["peptide_start"])
                candidates[key] = endpoint
            endpoint = deepcopy(next(iter(candidates.values())))
            if len(candidates) > 1:
                endpoint.update(source_protein_id=None, source_fasta_id=None,
                                protein_position=None, peptide_start=None, status="not_evaluable",
                                reason="protein_identity_ambiguous", candidates=list(candidates.values()))
            record["endpoints"][side] = endpoint

        if decoy_prefix is not None:
            categories = set()
            for pair in record["protein_pairs"]:
                a, b = (pair[side]["source_protein_id"].startswith(decoy_prefix) for side in ("a", "b"))
                categories.add("decoy_decoy" if a and b else "target_decoy" if a or b else "target_target")
            if len(categories) == 1:
                record["decoy_status"] = categories.pop()
            else:
                record["decoy_reason"] = "ambiguous_decoy_status"
        if any(endpoint["status"] != "valid" for endpoint in record["endpoints"].values()):
            record["record_reason"] = "source_endpoint_unresolved"
        elif record["decoy_status"] == "unknown":
            record["record_reason"] = record["decoy_reason"]
        else:
            record.update(record_status="valid", record_reason="")

    bundle = {"format": "plink_psm", "format_version": 1, "source_columns": columns,
            "records": records, "proteins": proteins,
            "parser": {"name": "pyXLMS", "version": version("pyXLMS"),
                       "entrypoint": "parse_sites_from_plink", "modification_masses_assessed": False,
                       "decoy_prefix": decoy_prefix, "fdr_columns": dict(fdr_columns)},
            "provenance": {"dataset_accession": dataset_accession, "psm_path": str(psm_path.resolve()),
                           "psm_sha256": digest, "fasta_path": str(fasta_path.resolve()),
                           "fasta_sha256": next(iter(proteins.values()))["fasta_sha256"]}}
    if is_workbook:
        bundle["parser"].update({
            "format_variant": "plink_psm_xlsx",
            "source_reader": "openpyxl",
            "numeric_source_values_repaired": False,
            "source_cells_retained": True,
        })
    if unfiltered_path is not None:
        return _add_plink_native_evidence(bundle, unfiltered_path, parameters_path)
    if author_filter_path is not None:
        return _add_author_filter_evidence(bundle, author_filter_path)
    return bundle


def write_import_report(bundle, output_dir):
    """Save the complete import without selecting records or applying identities."""
    output_dir = Path(output_dir)
    if output_dir.exists():
        raise FileExistsError("Import report output directory already exists")
    _verify_msannika_sources(bundle)
    _verify_plink_native_sources(bundle)
    author_sources = _author_filter_source_bodies(bundle)
    for kind in ("psm", "fasta"):
        source = Path(bundle["provenance"][kind + "_path"])
        digest = hashlib.sha256()
        with source.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        if digest.hexdigest() != bundle["provenance"][kind + "_sha256"]:
            raise ValueError(f"An import source changed after parsing: {source}")

    output_dir.mkdir(parents=True, exist_ok=False)
    _write_source_bodies(output_dir, author_sources)
    report = output_dir / "import_report.json"
    _write_json(report, bundle)
    records = bundle["records"]
    valid = sum(record["record_status"] == "valid" for record in records)
    errors = sum(record["record_reason"].startswith("pyxlms_") for record in records)
    manifest = {
        "schema_version": 1, "kind": "xlms_psm_import", "format": bundle["format"],
        "provenance": deepcopy(bundle["provenance"]), "parser": deepcopy(bundle.get("parser")),
        "counts": {"source_records": len(records), "validated_source_records": valid,
                   "unresolved_source_records": len(records) - valid, "parser_errors": errors},
        "record_reason_counts": dict(Counter(record["record_reason"] for record in records
                                             if record["record_reason"])),
        "fdr_filter_applied": False, "identity_mapping_applied": False,
        "input_exported": False, "assessment_performed": False,
        "output_sha256": {"import_report.json": hashlib.sha256(report.read_bytes()).hexdigest()},
    }
    if "native_filter_policy" in bundle:
        manifest["native_filter_policy"] = deepcopy(bundle["native_filter_policy"])
        manifest["native_evidence_counts"] = dict(Counter(record["native_evidence"]["status"] for record in records))
    if author_sources:
        manifest["author_filter_policy"] = deepcopy(bundle["author_filter_policy"])
        manifest["author_filter_evidence_snapshot"] = "author_filter_sources/" + Path(bundle["provenance"]["author_filter_path"]).name
        manifest["author_filter_record_counts"] = dict(Counter(record["source_filter_status"] for record in records))
        manifest["output_sha256"].update({name: hashlib.sha256(raw).hexdigest() for name, raw in author_sources.items()})
    _write_json(output_dir / "manifest.json", manifest)
    return output_dir


def target_proteins(bundle):
    """Select only proteins used by valid target PSMs, not the entire search FASTA."""
    used = {endpoint["source_fasta_id"] for row in bundle["records"]
            if row["record_status"] == "valid" and row["decoy_status"] == "target_target"
            for endpoint in row["endpoints"].values()}
    return {key: protein for key, protein in bundle["proteins"].items() if key in used}


def _taxon_id(value):
    if value is not None and (type(value) is not int or value <= 0):
        raise ValueError("taxon_id must be a positive integer or None")
    return value


def _identity(candidates, sources, *, reason="no_exact_sequence_match"):
    unique = {candidate["accession"]: candidate for candidate in candidates}
    candidates = [unique[key] for key in sorted(unique)]
    status = "mapped" if len(candidates) == 1 else "ambiguous" if candidates else "unmapped"
    match = candidates[0] if status == "mapped" else {}
    return {"status": status, "reason": "" if status == "mapped" else "multiple_accessions_for_sequence" if candidates else reason,
            "accession": match.get("accession"), "sequence": match.get("sequence"),
            "taxon_id": match.get("taxon_id"), "match_method": "full_sequence_exact",
            "candidates": candidates, "sources": deepcopy(sources)}


def map_reference_proteins(proteins, references, *, taxon_id=None):
    """Map against an explicitly supplied UniProt FASTA; never match by name alone."""
    taxon_id = _taxon_id(taxon_id)
    mappings = {}
    for identifier, source in proteins.items():
        candidates, receipts = [], []
        for reference in references.values():
            receipt = {"kind": "reference_fasta", "path": reference["fasta_path"],
                       "absolute_path": reference["fasta_path"], "sha256": reference["fasta_sha256"]}
            if receipt not in receipts:
                receipts.append(receipt)
            accession = reference["uniprot_accession"]
            if (not accession or reference["sequence"] != source["sequence"]
                    or taxon_id is not None and reference["taxon_id"] != taxon_id
                    or source["uniprot_accession"] and source["uniprot_accession"] != accession):
                continue
            candidates.append({"accession": accession, "sequence": reference["sequence"], "taxon_id": reference["taxon_id"]})
        mappings[identifier] = _identity(candidates, receipts)
    return mappings


def _lookup_get(url, sources_dir, cache, receipts, *, params=None, fasta=False):
    uri = urlsplit(url)
    if uri.scheme != "https" or uri.netloc != "rest.uniprot.org" or uri.username or uri.password:
        raise ValueError("Unexpected UniProt API address")
    key = json.dumps([url, params], sort_keys=True)
    if key in cache:
        data, response, receipt = cache[key]
    else:
        response = requests.get(url, params=params, headers={"Accept": "text/plain" if fasta else "application/json"}, timeout=(15, 120))
        response.raise_for_status()
        raw = response.content
        data = _fasta_records(raw.decode("utf-8-sig")) if fasta else json.loads(raw)
        if not isinstance(data, dict):
            raise ValueError("Unexpected UniProt response format")
        path = sources_dir / f"request_{len(cache) + 1:04d}{'.fasta' if fasta else '.json'}"
        with path.open("xb") as handle:
            handle.write(raw)
        receipt = {"kind": "uniprot_api", "request_url": url, "params": params or {},
                   "response_url": response.url, "response_headers": dict(response.headers),
                   "status": response.status_code, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                   "path": path.name, "absolute_path": str(path.resolve()), "sha256": hashlib.sha256(raw).hexdigest()}
        cache[key] = data, response, receipt
    if receipt not in receipts:
        receipts.append(receipt)
    return data, response


def resolve_uniprot_proteins(proteins, sources_dir, *, taxon_id=None):
    """Resolve exact full sequences via UniParc and active UniProtKB references.

    Both the UniParc sequence and the current UniProt FASTA must equal the
    source sequence. A failed lookup is retained as a technical failure.
    The caller explicitly chooses whether to constrain the organism.
    """
    taxon_id = _taxon_id(taxon_id)
    sources_dir = Path(sources_dir)
    sources_dir.mkdir(parents=True, exist_ok=False)
    cache, mappings = {}, {}
    for identifier, source in proteins.items():
        receipts, candidates = [], []
        try:
            sequence = source["sequence"]
            digest = hashlib.md5(sequence.encode("ascii")).hexdigest().upper()
            url, params, seen, entries = UNIPROT_REST + "/uniparc/search", {"query": "checksum:" + digest, "format": "json", "size": 100}, set(), []
            while url:
                if url in seen or len(seen) >= 10:
                    raise ValueError("UniParc pagination is repeated or incomplete")
                seen.add(url)
                data, response = _lookup_get(url, sources_dir, cache, receipts, params=params)
                if not isinstance(data.get("results"), list):
                    raise ValueError("UniParc search has no results list")
                entries.extend(data["results"])
                url, params = response.links.get("next", {}).get("url"), None
            accessions = set()
            exact_entries = active_references = 0
            for entry in entries:
                if entry.get("sequence", {}).get("value") != sequence:
                    continue
                upi = entry.get("uniParcId")
                if not isinstance(upi, str) or re.fullmatch(r"UPI[A-Fa-f0-9]{10}", upi) is None:
                    raise ValueError("UniParc returned an invalid entry identifier")
                full, _ = _lookup_get(UNIPROT_REST + "/uniparc/" + upi, sources_dir, cache, receipts)
                if full.get("uniParcId") != upi or full.get("sequence", {}).get("value") != sequence:
                    raise ValueError("UniParc entry does not match the queried sequence")
                exact_entries += 1
                references = full.get("uniParcCrossReferences")
                if not isinstance(references, list):
                    raise ValueError("UniParc entry has no cross-reference list")
                for reference in references:
                    accession = reference.get("id", "")
                    if (reference.get("active") is True
                            and str(reference.get("database", "")).startswith("UniProtKB/")
                            and isinstance(accession, str) and ACCESSION.fullmatch(accession)):
                        active_references += 1
                    if (reference.get("active") is not True or not str(reference.get("database", "")).startswith("UniProtKB/")
                            or ACCESSION.fullmatch(accession) is None
                            or taxon_id is not None and reference.get("organism", {}).get("taxonId") != taxon_id
                            or source["uniprot_accession"] and source["uniprot_accession"] != accession):
                        continue
                    accessions.add(accession)
            if len(accessions) > 100:
                raise ValueError("Too many active UniProt identities; specify an organism before resolving")
            for accession in sorted(accessions):
                current, _ = _lookup_get(UNIPROT_REST + "/uniprotkb/" + accession + ".fasta", sources_dir, cache, receipts, fasta=True)
                matches = [record for record in current.values()
                           if record["uniprot_accession"] == accession and record["sequence"] == sequence
                           and (taxon_id is None or record["taxon_id"] == taxon_id)]
                if len(matches) == 1:
                    candidates.append({"accession": accession, "sequence": sequence, "taxon_id": matches[0]["taxon_id"]})
            reason = (
                "sequence_not_found_in_uniparc" if exact_entries == 0 else
                "no_active_uniprot_reference" if active_references == 0 else
                "no_matching_active_uniprot_reference" if not accessions else
                "current_uniprot_sequence_mismatch"
            )
            mappings[identifier] = _identity(candidates, receipts, reason=reason)
            mappings[identifier]["lookup_details"] = {
                "exact_uniparc_entries": exact_entries,
                "active_uniprot_references": active_references,
                "candidate_accessions_after_filters": len(accessions),
                "current_exact_candidates": len(candidates),
                "requested_taxon_id": taxon_id,
            }
        except (requests.exceptions.RequestException, OSError, ValueError, KeyError, TypeError, AttributeError) as error:
            mappings[identifier] = _identity([], receipts)
            mappings[identifier].update(status="lookup_failed", reason=f"{type(error).__name__}: {error}")
    _write_json(sources_dir / "lookup_manifest.json", {"schema_version": 1, "taxon_id": taxon_id,
                "sources": [item[2] for item in cache.values()], "protein_status_counts": dict(Counter(m["status"] for m in mappings.values()))})
    return mappings


def _fdr_policy(level, cutoff):
    if level not in FDR_COLUMNS or isinstance(cutoff, bool) or not isinstance(cutoff, Real) or not math.isfinite(cutoff) or not 0 <= cutoff <= 1:
        raise ValueError("Choose an explicit FDR level and a finite fraction between 0 and 1")
    return {"level": level, "cutoff": float(cutoff), "boundary": "inclusive"}


def _protein_aligner():
    try:
        from Bio.Align import PairwiseAligner
    except ImportError as error:
        raise ImportError("Protein comparison requires the existing CLAUDIO dependency biopython") from error
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def _alignment_positions(alignment):
    mapping, blocks = {}, []
    for source_span, reference_span in zip(*alignment.aligned):
        start, end = map(int, source_span)
        ref_start, ref_end = map(int, reference_span)
        if end - start != ref_end - ref_start:
            raise ValueError("Unexpected alignment block lengths")
        blocks.append({"source": [start, end], "reference": [ref_start, ref_end]})
        mapping.update({i + 1: ref_start + i - start + 1 for i in range(start, end)})
    return mapping, blocks


def compare_protein_sequences(source, reference):
    """Describe a candidate alignment without deciding protein identity.

    Exact contiguous containment has an explicit coordinate map. Otherwise,
    coordinate maps are exposed only for a unique optimal global alignment.
    Identity includes all global alignment columns, including terminal gaps.
    """
    for sequence in (source, reference):
        if not isinstance(sequence, str) or not sequence or re.fullmatch(r"[A-Z]+", sequence) is None:
            raise ValueError("Comparison requires two nonempty, ungapped uppercase protein sequences")
    aligner = _protein_aligner()
    alignments = aligner.align(source, reference)
    alignment = alignments[0]
    try:
        count = len(alignments)
    except OverflowError:
        count = None
    representative, blocks = _alignment_positions(alignment)
    identical = sum(source[i - 1] == reference[j - 1] for i, j in representative.items())
    columns = int(alignment.shape[1])
    source_matches = [m.start() for m in re.finditer("(?=" + re.escape(source) + ")", reference)]
    reference_matches = [m.start() for m in re.finditer("(?=" + re.escape(reference) + ")", source)]
    if len(source_matches) == 1:
        offset = source_matches[0]
        mapping = {i + 1: i + offset + 1 for i in range(len(source))}
        method = "full_source_exact_subsequence"
    elif len(reference_matches) == 1:
        offset = reference_matches[0]
        mapping = {i + offset + 1: i + 1 for i in range(len(reference))}
        method = "full_reference_exact_subsequence"
    elif count == 1:
        mapping, method = representative, "unique_global_alignment"
    else:
        mapping, method = {}, "ambiguous_global_alignment"
    usable = bool(mapping)
    mismatches = [{"source_position": i, "reference_position": j,
                   "source_residue": source[i - 1], "reference_residue": reference[j - 1]}
                  for i, j in sorted(mapping.items()) if source[i - 1] != reference[j - 1]]
    return {
        "mapping_method": method, "coordinate_mapping_usable": usable,
        "source_length": len(source), "reference_length": len(reference),
        "alignment_columns": columns, "identical_residues": identical,
        "identity_percent": 100 * identical / columns,
        "identity_scope": "global_alignment_including_terminal_gaps",
        "source_coverage_percent": 100 * len(representative) / len(source),
        "reference_coverage_percent": 100 * len(representative) / len(reference),
        "optimal_alignment_count": count, "representative_alignment_blocks_zero_based": blocks,
        "aligned_position_map": {str(i): j for i, j in sorted(mapping.items())},
        "identical_position_map": {str(i): j for i, j in sorted(mapping.items()) if source[i - 1] == reference[j - 1]},
        "mismatches": mismatches,
        "source_unmapped_positions": sorted(set(range(1, len(source) + 1)) - set(mapping)) if usable else None,
        "reference_unmapped_positions": sorted(set(range(1, len(reference) + 1)) - set(mapping.values())) if usable else None,
        "alignment_parameters": {"mode": "global", "match": 2, "mismatch": -1, "gap_open": -10, "gap_extend": -0.5},
    }


def search_uniprot_candidates(proteins, sources_dir, *, taxon_id, gene_hints):
    """Find same-species gene candidates and retain sequence comparison evidence.

    Gene hints are caller-supplied search clues, not accepted identities. Every
    candidate stays in a separate review report; exact mappings are unchanged.
    """
    taxon_id = _taxon_id(taxon_id)
    if taxon_id is None or not isinstance(gene_hints, dict):
        raise ValueError("A taxon_id and explicit gene-hint dictionary are required")
    for identifier in proteins:
        hint = gene_hints.get(identifier)
        if hint is not None and (not isinstance(hint, str) or re.fullmatch(r"[A-Za-z0-9_.-]+", hint) is None):
            raise ValueError(f"Invalid gene query hint for {identifier}")
    sources_dir = Path(sources_dir)
    sources_dir.mkdir(parents=True, exist_ok=False)
    reports, cache = {}, {}
    for identifier, source in proteins.items():
        hint = gene_hints.get(identifier)
        report = {"status": "gene_hint_required", "reason": "no_gene_hint", "gene_hint": hint,
                  "taxon_id": taxon_id, "candidates": [], "sources": [], "identity_accepted": False}
        reports[identifier] = report
        if hint is None:
            continue
        try:
            url = UNIPROT_REST + "/uniprotkb/search"
            params = {"query": f"(gene_exact:{hint}) AND (organism_id:{taxon_id})", "format": "json", "size": 100}
            seen, records = set(), []
            while url:
                if url in seen or len(seen) >= 10:
                    raise ValueError("UniProt candidate pagination is repeated or incomplete")
                seen.add(url)
                data, response = _lookup_get(url, sources_dir, cache, report["sources"], params=params)
                if not isinstance(data.get("results"), list):
                    raise ValueError("UniProt candidate search has no results list")
                records.extend(data["results"])
                url, params = response.links.get("next", {}).get("url"), None
            used = set()
            for record in records:
                if record.get("organism", {}).get("taxonId") != taxon_id:
                    continue
                accession = record.get("primaryAccession")
                sequence = record.get("sequence", {}).get("value")
                if not isinstance(accession, str) or ACCESSION.fullmatch(accession) is None:
                    raise ValueError("UniProt candidate has no valid accession")
                if accession in used:
                    raise ValueError("UniProt candidate search repeated an accession")
                used.add(accession)
                report["candidates"].append({"accession": accession, "sequence": sequence, "taxon_id": taxon_id,
                    "entry_type": record.get("entryType"), "gene_annotations": record.get("genes", []),
                    "protein_description": record.get("proteinDescription", {}),
                    "comparison": compare_protein_sequences(source["sequence"], sequence)})
            report["status"] = "candidates_for_review" if report["candidates"] else "no_candidates"
            report["reason"] = "" if report["candidates"] else "no_candidate_for_gene_and_taxon"
        except (requests.exceptions.RequestException, OSError, ValueError, TypeError, KeyError, AttributeError, ImportError) as error:
            report.update(status="lookup_failed", reason=f"{type(error).__name__}: {error}")
    _write_json(sources_dir / "candidate_report.json", reports)
    _write_json(sources_dir / "candidate_manifest.json", {"schema_version": 1, "taxon_id": taxon_id,
                "identity_acceptance_performed": False, "status_counts": dict(Counter(r["status"] for r in reports.values())),
                "sources": [item[2] for item in cache.values()]})
    return reports


def _projection_coordinate_options(source, reference, maximum):
    for sequence in (source, reference):
        if not isinstance(sequence, str) or not sequence or re.fullmatch(r"[A-Z]+", sequence) is None:
            raise ValueError("Projection requires two nonempty, ungapped protein sequences")
    source_hits = [m.start() for m in re.finditer("(?=" + re.escape(source) + ")", reference)]
    reference_hits = [m.start() for m in re.finditer("(?=" + re.escape(reference) + ")", source)]
    if len(source_hits) == 1:
        offset = source_hits[0]
        mapping = {i + 1: i + offset + 1 for i in range(len(source))}
        blocks = [[{"source": [0, len(source)], "reference": [offset, offset + len(source)]}]]
        return [mapping], {"method": "full_source_exact_subsequence", "optimal_alignment_count": None,
                           "reason": "", "mapping_option_count": 1, "alignment_blocks_zero_based": blocks}
    if len(reference_hits) == 1:
        offset = reference_hits[0]
        mapping = {i + offset + 1: i + 1 for i in range(len(reference))}
        blocks = [[{"source": [offset, offset + len(reference)], "reference": [0, len(reference)]}]]
        return [mapping], {"method": "full_reference_exact_subsequence", "optimal_alignment_count": None,
                           "reason": "", "mapping_option_count": 1, "alignment_blocks_zero_based": blocks}
    alignments = _protein_aligner().align(source, reference)
    try:
        count = len(alignments)
    except OverflowError:
        count = None
    evidence = {"method": "optimal_global_alignments", "optimal_alignment_count": count,
                "reason": "", "mapping_option_count": 0, "alignment_blocks_zero_based": []}
    if count is None or count > maximum:
        evidence["reason"] = "optimal_alignment_limit_exceeded"
        return [], evidence
    options = []
    for alignment in alignments:
        mapping, blocks = _alignment_positions(alignment)
        options.append(mapping)
        evidence["alignment_blocks_zero_based"].append(blocks)
    evidence["mapping_option_count"] = len(options)
    return options, evidence


def _project_peptide(endpoint, profile, coordinate_options):
    evidence = profile["coordinate_evidence"]
    result = {"profile_id": profile["profile_id"], "accession": profile["accession"],
              "identity_scope": profile["identity_scope"], "mapping_method": evidence["method"],
              "source_protein_position": endpoint["protein_position"],
              "source_peptide_start": endpoint["peptide_start"], "link_position": endpoint["link_position"],
              "reference_protein_position": None, "reference_peptide_start": None,
              "reference_peptide_end": None, "status": "not_evaluable", "reason": "",
              "optimal_alignment_count": evidence["optimal_alignment_count"],
              "mapping_option_count": len(coordinate_options), "compatible_alignment_count": 0,
              "distinct_projection_count": 0, "projection_options": [], "rejected_alignment_reasons": {}}
    if evidence["reason"]:
        result["reason"] = evidence["reason"]
        return result
    peptide = endpoint["peptide"]
    start, link, source_site = endpoint["peptide_start"], endpoint["link_position"], endpoint["protein_position"]
    if (endpoint["status"] != "valid" or any(type(v) is not int or v < 1 for v in (start, link, source_site))
            or not isinstance(peptide, str) or not 1 <= link <= len(peptide)
            or source_site != start + link - 1
            or profile["source_sequence"][start - 1:start - 1 + len(peptide)] != peptide):
        result["reason"] = "source_peptide_validation_failed"
        return result
    projections, rejected = Counter(), Counter()
    reference = profile["reference_sequence"]
    for mapping in coordinate_options:
        positions = [mapping.get(i) for i in range(start, start + len(peptide))]
        if any(position is None for position in positions):
            rejected["peptide_overlaps_unmapped_source_residues"] += 1
        elif any(position != positions[0] + offset for offset, position in enumerate(positions)):
            rejected["peptide_spans_reference_insertion"] += 1
        elif reference[positions[0] - 1:positions[-1]] != peptide:
            rejected["peptide_sequence_differs_from_reference"] += 1
        else:
            projections[(positions[0], positions[-1], positions[link - 1])] += 1
    result["compatible_alignment_count"] = sum(projections.values())
    result["distinct_projection_count"] = len(projections)
    result["rejected_alignment_reasons"] = dict(rejected)
    result["projection_options"] = [
        {"peptide_start": key[0], "peptide_end": key[1], "protein_position": key[2], "supporting_mapping_count": value}
        for key, value in sorted(projections.items())
    ]
    if len(projections) == 1:
        peptide_start, peptide_end, site = next(iter(projections))
        result.update(status="projected", reference_protein_position=site,
                      reference_peptide_start=peptide_start, reference_peptide_end=peptide_end)
    elif projections:
        result["reason"] = "ambiguous_reference_position"
    else:
        result["reason"] = rejected.most_common(1)[0][0] if rejected else "no_coordinate_mapping"
    return result


def _selected_candidate_reports(candidate_reports):
    history = candidate_reports if isinstance(candidate_reports, list) else [candidate_reports]
    selected, selected_indices, contexts = {}, {}, {}
    completed = {"candidates_for_review", "no_candidates"}
    for index, snapshot in enumerate(history):
        if not isinstance(snapshot, dict):
            raise ValueError("Candidate report history must contain dictionaries")
        for identifier, report in snapshot.items():
            if not isinstance(report, dict):
                raise ValueError("Candidate report entries must be dictionaries")
            context = (report.get("gene_hint"), report.get("taxon_id"))
            if all(value is not None for value in context):
                if identifier in contexts and contexts[identifier] != context:
                    raise ValueError("Candidate report query contexts differ for the same protein")
                contexts[identifier] = context
            if (identifier not in selected or report.get("status") in completed
                    or selected[identifier].get("status") not in completed):
                selected[identifier], selected_indices[identifier] = report, index
    return selected, selected_indices, history


def project_xifdr_candidates(
    bundle, exact_mappings, candidate_reports, *,
    max_optimal_alignments=64,
    export_unverified_candidates=False,
):
    """Project intact peptides onto each reference hypothesis; retain every PSM.

    A peptide projection must have one distinct, contiguous, sequence-identical
    position among the compatible mappings. Multiple protein identities stay
    separate hypotheses. This stage neither filters FDR nor accepts identities.
    The alignment limit bounds computation; exceeding it leaves an unresolved
    audit record rather than selecting the first alignment.
    Candidate reports may be a dictionary or snapshots in retrieval order;
    the last completed result is reused while failed attempts stay in evidence.
    """
    if type(export_unverified_candidates) is not bool:
        raise ValueError("export_unverified_candidates must be a boolean")
    if type(max_optimal_alignments) is not int or max_optimal_alignments < 1:
        raise ValueError("max_optimal_alignments must be a positive integer")
    candidate_reports, report_indices, report_history = _selected_candidate_reports(candidate_reports)
    if export_unverified_candidates:
        candidate_prepared = _prepare_plink_candidate_input(
            bundle, exact_mappings
        )
        candidate_entries = {
            entry["source_record_id"]: entry
            for entry in candidate_prepared["audit"]
        }
        identity_reasons = {
            "protein_identity_unresolved_a",
            "protein_identity_unresolved_b",
        }
        candidate_eligible = {
            identifier
            for identifier, entry in candidate_entries.items()
            if not any(
                reason not in identity_reasons
                for reason in entry["input_reasons"]
            )
        }
        used = {
            record["endpoints"][site]["source_fasta_id"]
            for record in bundle["records"]
            if record["source_record_id"] in candidate_eligible
            for site in ("a", "b")
        }
        proteins = {
            identifier: source
            for identifier, source in bundle["proteins"].items()
            if identifier in used
        }
    else:
        proteins = target_proteins(bundle)
    profiles, options, lookup_states = {}, {}, {}
    reference_sequences = {}
    for identifier, source in proteins.items():
        exact = exact_mappings.get(identifier, {})
        report = candidate_reports.get(identifier, {})
        if exact.get("status") == "mapped":
            if exact.get("sequence") != source["sequence"]:
                raise ValueError("An exact identity no longer matches its source sequence")
            candidates = [{"accession": exact["accession"], "sequence": exact["sequence"], "taxon_id": exact["taxon_id"]}]
            scope, sources = "exact_sequence", exact.get("sources", [])
        else:
            candidates = report.get("candidates", []) if report.get("status") == "candidates_for_review" else []
            scope, sources = "candidate_hypothesis", report.get("sources", [])
        lookup_states[identifier] = {"exact_mapping_status": exact.get("status", "not_attempted"),
                                     "candidate_search_status": report.get("status", "not_attempted"),
                                     "candidate_search_reason": report.get("reason", ""),
                                     "selected_candidate_report_index": report_indices.get(identifier)}
        profiles[identifier] = []
        used = set()
        for candidate in candidates:
            accession, sequence = candidate.get("accession"), candidate.get("sequence")
            if not isinstance(accession, str) or ACCESSION.fullmatch(accession) is None or accession in used:
                raise ValueError("Invalid or repeated reference candidate accession")
            used.add(accession)
            if accession in reference_sequences and reference_sequences[accession] != sequence:
                raise ValueError("The same reference accession has conflicting sequence snapshots")
            reference_sequences[accession] = sequence
            coordinates, evidence = _projection_coordinate_options(source["sequence"], sequence, max_optimal_alignments)
            digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            profile_id = identifier + ":" + accession + ":" + digest
            profile = {"profile_id": profile_id, "accession": accession, "identity_scope": scope,
                       "source_fasta_id": identifier, "source_sequence": source["sequence"],
                       "reference_sequence": sequence, "reference_sequence_sha256": digest,
                       "taxon_id": candidate.get("taxon_id"), "gene_hint": report.get("gene_hint"),
                       "coordinate_evidence": evidence, "sources": deepcopy(sources)}
            profiles[identifier].append(profile)
            options[profile_id] = coordinates
    records = []
    for source_record in bundle["records"]:
        fields, endpoints = source_record["source_fields"], source_record["endpoints"]
        row = {"source_record_id": source_record["source_record_id"], "source_record_number": source_record["source_record_number"],
               "source_psm_id": fields.get("PSMID", ""), "source_run": fields.get("run", ""),
               "status": "not_evaluable", "reason": "", "endpoint_candidates": {"a": [], "b": []},
               "pairs": [], "discarded_cross_identity_pair_count": 0}
        records.append(row)
        if source_record["decoy_status"] not in ("target_target", "unknown"):
            row.update(status="excluded_decoy", reason="decoy_record")
            continue
        if export_unverified_candidates:
            identifier = source_record["source_record_id"]
            if identifier not in candidate_eligible:
                blocking = [
                    reason
                    for reason in candidate_entries[identifier]["input_reasons"]
                    if reason not in identity_reasons
                ]
                row["reason"] = (
                    blocking[0] if blocking else "source_record_unresolved"
                )
                continue
        elif (
            source_record["record_status"] != "valid"
            or source_record["decoy_status"] != "target_target"
        ):
            row["reason"] = (
                source_record["record_reason"] or "invalid_decoy_flags"
            )
            continue
        for site in ("a", "b"):
            endpoint = endpoints[site]
            for profile in profiles.get(endpoint["source_fasta_id"], []):
                row["endpoint_candidates"][site].append(_project_peptide(endpoint, profile, options[profile["profile_id"]]))
        for first in row["endpoint_candidates"]["a"]:
            if first["status"] != "projected":
                continue
            for second in row["endpoint_candidates"]["b"]:
                if second["status"] != "projected":
                    continue
                if endpoints["a"]["source_fasta_id"] == endpoints["b"]["source_fasta_id"] and first["accession"] != second["accession"]:
                    row["discarded_cross_identity_pair_count"] += 1
                    continue
                key = [source_record["source_record_id"], first["profile_id"], first["reference_protein_position"],
                       second["profile_id"], second["reference_protein_position"]]
                pair = {"projection_candidate_id": "projection_" + hashlib.sha256(json.dumps(key).encode()).hexdigest(),
                        "unip_id_a": first["accession"], "unip_id_b": second["accession"],
                        "profile_id_a": first["profile_id"], "profile_id_b": second["profile_id"],
                        "reference_position_a": first["reference_protein_position"], "reference_position_b": second["reference_protein_position"],
                        "reference_peptide_start_a": first["reference_peptide_start"], "reference_peptide_start_b": second["reference_peptide_start"],
                        "identity_scope": "exact_sequence" if first["identity_scope"] == second["identity_scope"] == "exact_sequence" else "candidate_hypothesis"}
                row["pairs"].append(pair)
        if row["pairs"]:
            row["status"] = "projected_candidates"
        else:
            row["reason"] = "no_compatible_reference_pair"
    return {"schema_version": 1, "kind": "xifdr_reference_projection", "provenance": deepcopy(bundle["provenance"]),
            "projection_rule_version": "intact-peptide-projection-v1",
            "max_optimal_alignments": max_optimal_alignments,
            "export_unverified_candidates": export_unverified_candidates,
            "protein_identity_acceptance_performed": False, "model_assessment_performed": False, "fdr_filter_applied": False,
            "protein_candidates": profiles, "protein_lookup_states": lookup_states, "records": records,
            "identity_evidence": {"exact_mappings": deepcopy(exact_mappings), "candidate_reports": deepcopy(candidate_reports),
                                  "candidate_report_history": deepcopy(report_history),
                                  "report_selection_policy": "last_completed_report_else_last_outcome"},
            "counts": {"source_records": len(records), "target_records": sum(r["decoy_status"] == "target_target" for r in bundle["records"]),
                       "records_with_projected_candidates": sum(r["status"] == "projected_candidates" for r in records),
                       "projection_pairs": sum(len(r["pairs"]) for r in records),
                       "record_status_counts": dict(Counter(r["status"] for r in records))}}


def write_projection_outputs(bundle, projection, output_dir):
    """Save candidate projections, all original PSMs, and a flat review table."""
    if (projection["provenance"] != bundle["provenance"]
            or [r["source_record_id"] for r in projection["records"]] != [r["source_record_id"] for r in bundle["records"]]):
        raise ValueError("Projection outputs do not match their original PSM records")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    _write_json(output_dir / "projection_report.json", projection)
    _write_json(output_dir / "source_psms.json", bundle)
    columns = ["source_record_id", "source_record_number", "source_psm_id", "source_run", "dataset_accession",
               "source_protein_id_a", "source_protein_id_b", "source_position_a", "source_position_b",
               "projection_candidate_id", "unip_id_a", "unip_id_b", "profile_id_a", "profile_id_b",
               "reference_position_a", "reference_position_b", "reference_peptide_start_a", "reference_peptide_start_b", "identity_scope"]
    with (output_dir / "projection_pairs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for original, projected in zip(bundle["records"], projection["records"]):
            for pair in projected["pairs"]:
                row = {key: projected[key] for key in ("source_record_id", "source_record_number", "source_psm_id", "source_run")}
                row.update(dataset_accession=bundle["provenance"]["dataset_accession"])
                for site in ("a", "b"):
                    endpoint = original["endpoints"][site]
                    row[f"source_protein_id_{site}"] = endpoint["source_protein_id"]
                    row[f"source_position_{site}"] = endpoint["protein_position"]
                writer.writerow({**row, **pair})
    names = ["projection_report.json", "source_psms.json", "projection_pairs.csv"]
    _write_json(output_dir / "manifest.json", {"schema_version": 1, "kind": projection["kind"], "status": "complete",
                "provenance": projection["provenance"], "counts": projection["counts"],
                "projection_rule_version": projection["projection_rule_version"],
                "max_optimal_alignments": projection["max_optimal_alignments"],
                "protein_identity_acceptance_performed": False, "model_assessment_performed": False, "fdr_filter_applied": False,
                "output_sha256": {name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest() for name in names}})
    return output_dir


def _input_row(bundle, record, identities, policy):
    """Build shared source metadata without deciding identity or coordinate eligibility."""
    author = policy.get("kind") == "author_filtered_psm_table"
    unverified = policy.get("kind") == "unverified_candidates"
    level = None if author or unverified else policy["level"]
    fields = record["source_fields"]
    plink = bundle.get("format") == "plink_psm"
    row = {"source_record_id": record["source_record_id"], "dataset_accession": bundle["provenance"]["dataset_accession"],
           "source_psm_sha256": bundle["provenance"]["psm_sha256"], "source_fasta_sha256": bundle["provenance"]["fasta_sha256"],
           "source_record_number": record["source_record_number"], "source_psm_id": fields.get("Order" if plink else "PSMID", ""),
           "source_run": fields.get("run", ""), "source_scan": fields.get("scan", ""),
           "source_crosslinker": fields.get("Linker" if plink else "Crosslinker", ""), "source_fdr_group": fields.get("fdrGroup", ""),
           "source_fields_json": json.dumps(fields, ensure_ascii=False),
           "selected_fdr_level": level,
           "selected_fdr_value": (
               None if author or unverified
               else record["fdr"][level]["value"]
           ),
           "selected_fdr_cutoff": (
               None if author or unverified
               else float(policy["cutoff"])
           )}
    if unverified:
        row.update(
            input_scope="unverified_candidate",
            fdr_verification_status="unverified",
            source_decoy_status=record["decoy_status"],
            selection_basis="sequence_and_site_checks_only",
        )
    if author:
        row.update(selection_basis="author_filtered_psm_table", source_filter_level=policy["level"],
                   source_filter_cutoff=policy["reported_cutoff"],
                   source_filter_evidence_sha256=bundle["provenance"]["author_filter_sha256"])
    if plink:
        row.update(source_spectrum_title=fields.get("Title", ""), source_modifications=fields.get("Modifications", ""))
    for site in ("a", "b"):
        endpoint, identity = record["endpoints"][site], identities[site]
        row.update({f"unip_id_{site}": identity["accession"], f"seq_{site}": identity["sequence"],
                    f"pep_{site}": endpoint["peptide"], f"pos_{site}": endpoint["protein_position"],
                    f"res_pos_{site}": endpoint["link_position"], f"position_status_{site}": "valid",
                    f"source_protein_id_{site}": endpoint["source_protein_id"],
                    f"source_protein_position_{site}": endpoint["protein_position"],
                    f"source_modified_peptide_{site}": endpoint["source_peptide"]})
    return row


def prepare_af_input(bundle, mappings, *, fdr_level=None, fdr_cutoff=None):
    """Project validated target PSMs onto exact UniProt identities, retaining all audit rows."""
    author = "author_filter_policy" in bundle
    if author:
        if fdr_level is not None or fdr_cutoff is not None:
            raise ValueError("Author filtering cannot be combined with per-record FDR selection")
        _author_filter_source_bodies(bundle)
        policy = {"kind": "author_filtered_psm_table", "level": bundle["author_filter_policy"]["filter_level"],
                  "reported_cutoff": bundle["author_filter_policy"]["reported_fdr_cutoff"],
                  "per_record_fdr_applied": False}
    else:
        policy = _fdr_policy(fdr_level, fdr_cutoff)
    rows, audit = [], []
    for record in bundle["records"]:
        entry, reasons = deepcopy(record), []
        if record["decoy_status"] == "unknown":
            reasons.append(record.get("decoy_reason") or "invalid_decoy_flags")
        elif record["decoy_status"] != "target_target":
            reasons.append("decoy_record")
        if record["record_status"] != "valid":
            reasons.append(record["record_reason"])
        if author:
            if record.get("source_filter_status") != "reported_accepted":
                reasons.append("author_filter_acceptance_unresolved")
        else:
            confidence = record["fdr"][fdr_level]
            if confidence["status"] != "valid":
                reasons.append("selected_fdr_" + confidence["status"])
            elif confidence["value"] > fdr_cutoff:
                reasons.append("selected_fdr_above_cutoff")
        selected = {}
        for site, endpoint in record["endpoints"].items():
            identity = mappings.get(endpoint["source_fasta_id"])
            source = bundle["proteins"].get(endpoint["source_fasta_id"])
            if (not identity or identity["status"] != "mapped" or not source
                    or not identity.get("accession") or ACCESSION.fullmatch(identity["accession"]) is None
                    or identity.get("sequence") != source["sequence"]):
                reasons.append("protein_identity_unresolved_" + site)
            else:
                selected[site] = identity
        entry.update(input_status="not_evaluable" if reasons else "ready", input_reasons=list(dict.fromkeys(reasons)))
        audit.append(entry)
        if reasons:
            continue
        rows.append(_input_row(bundle, record, selected, policy))
    return {"rows": rows, "audit": audit, "fdr_policy": policy, "provenance": deepcopy(bundle["provenance"])}

def _verify_reference_accession_review(profile, report, source):
    """Verify a source-bound review of a different reference accession."""
    review = report.get("reference_accession_review")
    if not isinstance(review, dict):
        raise ValueError(
            "Reference projection requires the declared source UniProt accession"
        )

    source_sequence = source.get("sequence")
    reference_sequence = profile.get("reference_sequence")
    if any(
        not isinstance(sequence, str)
        or re.fullmatch(r"[A-Z]+", sequence) is None
        for sequence in (source_sequence, reference_sequence)
    ):
        raise ValueError("Reference accession review requires valid sequences")

    source_digest = hashlib.sha256(
        source_sequence.encode("ascii")
    ).hexdigest()
    reference_digest = hashlib.sha256(
        reference_sequence.encode("ascii")
    ).hexdigest()

    if (
        profile.get("source_fasta_id") != source.get("id")
        or profile.get("source_sequence") != source_sequence
        or profile.get("reference_sequence_sha256") != reference_digest
    ):
        raise ValueError("Reference accession review profile differs from its source")

    expected = {
        "schema_version": 1,
        "review_status": "reviewed",
        "source_fasta_id": source["id"],
        "source_fasta_sha256": source.get("fasta_sha256"),
        "source_sequence_sha256": source_digest,
        "declared_accession": source.get("uniprot_accession"),
        "reference_accession": profile["accession"],
        "reference_sequence_sha256": reference_digest,
        "taxon_id": profile["taxon_id"],
    }
    if (
        type(review.get("schema_version")) is not int
        or any(review.get(key) != value for key, value in expected.items())
        or not isinstance(expected["source_fasta_sha256"], str)
        or re.fullmatch(
            r"[0-9a-f]{64}", expected["source_fasta_sha256"]
        ) is None
        or not isinstance(review.get("reason"), str)
        or not review["reason"].strip()
    ):
        raise ValueError(
            "Reference accession review does not match the source and reference"
        )

    occurrences = sum(
        1 for _ in re.finditer(
            "(?=" + re.escape(reference_sequence) + ")", source_sequence
        )
    )
    if occurrences != 1:
        raise ValueError(
            "Reviewed accession changes require one complete reference "
            "occurrence in the source sequence"
        )

    receipts = review.get("sources")
    if not isinstance(receipts, list) or not receipts:
        raise ValueError("Reference accession review requires evidence sources")

    for receipt in receipts:
        if (
            not isinstance(receipt, dict)
            or receipt.get("kind") != "reference_accession_review"
            or not isinstance(receipt.get("absolute_path"), str)
            or not receipt["absolute_path"]
            or not Path(receipt["absolute_path"]).is_absolute()
            or not isinstance(receipt.get("sha256"), str)
            or re.fullmatch(r"[0-9a-f]{64}", receipt["sha256"]) is None
        ):
            raise ValueError("Invalid reference accession review evidence")

        raw = Path(receipt["absolute_path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != receipt["sha256"]:
            raise ValueError("Reference accession review evidence checksum changed")

    return review
def _reference_evidence(
    profile, report, source, *, require_declared_accession=True
):
    """Verify reference evidence and any explicitly reviewed accession change."""
    declared = source.get("uniprot_accession")
    if require_declared_accession and not declared:
        raise ValueError(
            "Reference projection requires the declared source UniProt accession"
        )
    if declared and profile["accession"] != declared:
        _verify_reference_accession_review(profile, report, source)

    taxon = profile["taxon_id"]
    if (type(taxon) is not int or taxon <= 0 or report.get("taxon_id") != taxon
            or source.get("taxon_id") is not None and source["taxon_id"] != taxon):
        raise ValueError("Reference projection organism differs from its source or query")
    matched = False
    for receipt in report.get("sources", []):
        raw = Path(receipt["absolute_path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != receipt["sha256"]:
            raise ValueError("Reference evidence checksum changed after projection")
        kind = receipt.get("kind")
        if kind == "uniprot_api":
            uri = urlsplit(receipt.get("response_url", receipt.get("request_url", "")))
            if (uri.scheme != "https" or uri.netloc != "rest.uniprot.org"
                    or not uri.path.startswith("/uniprotkb/") or receipt.get("status") != 200):
                continue
        if kind == "reference_fasta" or kind == "uniprot_api" and uri.path.endswith(".fasta"):
            entries = _fasta_records(raw.decode("utf-8-sig")).values()
            matched |= any(entry["uniprot_accession"] == profile["accession"]
                           and entry["taxon_id"] == taxon and entry["sequence"] == profile["reference_sequence"]
                           for entry in entries)
        elif kind == "uniprot_api":
            payload = json.loads(raw)
            entries = payload.get("results", [payload])
            if not isinstance(entries, list):
                raise ValueError("Invalid reference evidence response")
            matched |= any(isinstance(entry, dict) and entry.get("primaryAccession") == profile["accession"]
                           and entry.get("organism", {}).get("taxonId") == taxon
                           and entry.get("sequence", {}).get("value") == profile["reference_sequence"]
                           for entry in entries)
    if not matched:
        raise ValueError("Reference sequence and organism are not confirmed by the saved evidence")


def _prepare_reference_input(bundle, mappings, prepared, projection):
    """Recompute candidate coordinates, retaining the original identity and filtering audit."""
    if (not isinstance(projection, dict) or projection.get("kind") != "xifdr_reference_projection"
            or projection.get("projection_rule_version") != "intact-peptide-projection-v1"):
        raise ValueError("Unsupported reference projection report")
    candidate_mode = (
        prepared["fdr_policy"].get("kind") == "unverified_candidates"
    )
    projection_mode = projection.get("export_unverified_candidates", False)
    if type(projection_mode) is not bool or projection_mode != candidate_mode:
        raise ValueError(
            "Reference projection candidate mode differs from the export mode"
        )
    source_context = {k: v for k, v in bundle["provenance"].items() if not k.endswith("_path")}
    projection_context = {k: v for k, v in projection.get("provenance", {}).items() if not k.endswith("_path")}
    if projection_context != source_context:
        raise ValueError("Reference projection does not match the input source checksums or policies")
    for kind in ("psm", "fasta"):
        raw = Path(bundle["provenance"][kind + "_path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != bundle["provenance"][kind + "_sha256"]:
            raise ValueError("An input source changed after reference projection")
    evidence = projection.get("identity_evidence", {})
    reports = evidence.get("candidate_reports")
    if not isinstance(reports, dict):
        raise ValueError("Reference projection has no candidate evidence")
    checked = project_xifdr_candidates(
        bundle, mappings, reports,
        max_optimal_alignments=projection.get("max_optimal_alignments"),
        export_unverified_candidates=candidate_mode,
    )
    if checked["records"] != projection.get("records"):
        raise ValueError("Reference projection record order or coordinates differ from recomputed results")
    profiles = {profile["profile_id"]: profile
                for values in checked["protein_candidates"].values() for profile in values}
    original_profiles = {profile["profile_id"]: profile
                         for values in projection.get("protein_candidates", {}).values() for profile in values}
    if (profiles.keys() != original_profiles.keys()
            or any({k: v for k, v in profile.items() if k != "sources"}
                   != {k: v for k, v in original_profiles[key].items() if k != "sources"}
                   for key, profile in profiles.items())):
        raise ValueError("Reference projection profiles differ from recomputed results")
    for profile in profiles.values():
        if profile["identity_scope"] != "exact_sequence":
            identifier = profile["source_fasta_id"]
            _reference_evidence(profile, reports[identifier], bundle["proteins"][identifier])

    rows = []
    for record, entry, projected in zip(bundle["records"], prepared["audit"], checked["records"]):
        # Projection cannot repair source, decoy, author-filter or numeric FDR failures.
        reasons = [reason for reason in entry["input_reasons"]
                   if not reason.startswith("protein_identity_unresolved_")]
        if reasons:
            continue
        if len(projected["pairs"]) != 1:
            entry["input_reasons"].append("reference_projection_unresolved")
            entry["input_status"] = "not_evaluable"
            continue
        pair = projected["pairs"][0]
        selected = {site: profiles[pair[f"profile_id_{site}"]] for site in ("a", "b")}
        identities = {site: {"accession": profile["accession"], "sequence": profile["reference_sequence"]}
                      for site, profile in selected.items()}
        row = _input_row(bundle, record, identities, prepared["fdr_policy"])
        row.update(projection_candidate_id=pair["projection_candidate_id"], identity_scope=pair["identity_scope"])
        reference_endpoints = {}
        for site, profile in selected.items():
            candidate = next(candidate for candidate in projected["endpoint_candidates"][site]
                             if candidate["profile_id"] == profile["profile_id"])
            reference_endpoints[site] = deepcopy(candidate)
            row.update({f"pos_{site}": candidate["reference_protein_position"],
                        f"source_peptide_start_{site}": candidate["source_peptide_start"],
                        f"reference_peptide_start_{site}": candidate["reference_peptide_start"],
                        f"reference_mapping_method_{site}": candidate["mapping_method"],
                        f"reference_identity_scope_{site}": candidate["identity_scope"]})
        entry.update(input_status="ready", input_reasons=[], reference_endpoints=reference_endpoints,
                     projection_candidate_id=pair["projection_candidate_id"], identity_scope=pair["identity_scope"])
        rows.append(row)
    projection_policy = {
        "scope": "declared_source_accession",
        "intact_peptide_required": True,
        "unique_coordinates_required": True,
        "identity_acceptance_performed": False,
    }
    accession_changes = []
    for profile in profiles.values():
        identifier = profile["source_fasta_id"]
        source = bundle["proteins"][identifier]
        declared = source.get("uniprot_accession")
        if declared and declared != profile["accession"]:
            review = reports[identifier]["reference_accession_review"]
            accession_changes.append({
                "source_fasta_id": identifier,
                "profile_id": profile["profile_id"],
                "declared_accession": declared,
                "reference_accession": profile["accession"],
                "source_fasta_sha256": review["source_fasta_sha256"],
                "source_sequence_sha256": review["source_sequence_sha256"],
                "reference_sequence_sha256": review[
                    "reference_sequence_sha256"
                ],
                "reason": review["reason"],
                "review_evidence_sha256": [
                    receipt["sha256"] for receipt in review["sources"]
                ],
            })
    if accession_changes:
        projection_policy.update({
            "scope": "declared_or_reviewed_source_accession",
            "reviewed_accession_changes": accession_changes,
        })

    prepared.update(
        rows=rows,
        reference_projection=deepcopy(projection),
        reference_projection_policy=projection_policy,
    )
    return prepared

def _prepare_plink_candidate_input(bundle, mappings):
    """Prepare sequence-validated candidates without accepting FDR or decoy status."""
    if bundle.get("format") != "plink_psm":
        raise ValueError("Unverified candidate export currently supports pLink PSMs only")
    if (
        "author_filter_policy" in bundle
        or "native_filter_policy" in bundle
    ):
        raise ValueError(
            "Unverified candidate export cannot be combined with an existing "
            "author or native filter policy"
        )

    for kind in ("psm", "fasta"):
        path = Path(bundle["provenance"][kind + "_path"])
        expected = bundle["provenance"][kind + "_sha256"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(
                f"Input source checksum changed before candidate export: {kind}"
            )

    policy = {
        "kind": "unverified_candidates",
        "level": None,
        "cutoff": None,
        "per_record_fdr_applied": False,
        "fdr_verification_status": "unverified",
        "decoy_acceptance_performed": False,
        "identity_acceptance_performed": False,
    }
    rows, audit = [], []

    for record in bundle["records"]:
        entry = deepcopy(record)
        reasons = []
        decoy = record.get("decoy_status")

        if decoy in {"target_decoy", "decoy_decoy"}:
            reasons.append("decoy_record")
        elif decoy == "unknown":
            if (
                record.get("decoy_reason")
                != "decoy_policy_not_provided"
            ):
                reasons.append("decoy_status_unresolved")
        elif decoy != "target_target":
            reasons.append("invalid_decoy_flags")

        missing_policy_only = (
            record.get("record_reason") == "decoy_policy_not_provided"
            and decoy == "unknown"
            and record.get("decoy_reason") == "decoy_policy_not_provided"
        )
        if record.get("record_status") != "valid" and not missing_policy_only:
            reasons.append(
                record.get("record_reason") or "source_record_unresolved"
            )

        endpoints = record.get("endpoints", {})
        selected = {}
        for site in ("a", "b"):
            endpoint = endpoints.get(site)
            if not isinstance(endpoint, dict) or endpoint.get("status") != "valid":
                reasons.append("source_endpoint_unresolved_" + site)
                continue

            identifier = endpoint.get("source_fasta_id")
            source = bundle["proteins"].get(identifier)
            identity = mappings.get(identifier)
            if (
                not source
                or not isinstance(identity, dict)
                or identity.get("status") != "mapped"
                or not isinstance(identity.get("accession"), str)
                or ACCESSION.fullmatch(identity["accession"]) is None
                or identity.get("sequence") != source["sequence"]
                or (
                    source.get("uniprot_accession")
                    and identity["accession"] != source["uniprot_accession"]
                )
                or (
                    source.get("taxon_id") is not None
                    and identity.get("taxon_id") != source["taxon_id"]
                )
            ):
                reasons.append("protein_identity_unresolved_" + site)
            else:
                selected[site] = identity

        entry.update(
            input_status="not_evaluable" if reasons else "candidate_ready",
            input_reasons=list(dict.fromkeys(reasons)),
        )
        audit.append(entry)

        if not reasons:
            rows.append(_input_row(bundle, record, selected, policy))

    return {
        "rows": rows,
        "audit": audit,
        "fdr_policy": policy,
        "provenance": deepcopy(bundle["provenance"]),
    }
def write_claudio_input(bundle, mappings, output_dir, *, fdr_level=None, fdr_cutoff=None,
                        allowed_crosslinkers, allowed_residues, source_protein_ids=None,
                        reference_projection=None, export_unverified_candidates=False):
    """Export exact identities or explicitly selected reference projections for the full CLI.

    The CSV uses M1's default input column names.
    Selection is explicit and every original PSM remains in the audit. The
    original PSM key is named xlms_source_record_id because M1 creates its own
    source_record_id. Use the AF pipeline's source-preserving preprocessing.
    Reference projections keep identity hypotheses separate from exact mappings.
    This exports inputs only; it does not retrieve or assess an AF model.
    """
    selections = {}
    for name, values in (("allowed_crosslinkers", allowed_crosslinkers),
                         ("allowed_residues", allowed_residues),
                         ("source_protein_ids", source_protein_ids)):
        if name == "source_protein_ids" and values is None:
            selections[name] = None
            continue
        if (not isinstance(values, (set, frozenset, list, tuple)) or not values
                or any(not isinstance(value, str) or not value.strip() for value in values)):
            raise ValueError(f"{name} must be a nonempty collection of explicit labels")
        values = set(values)
        if name == "allowed_residues" and any(re.fullmatch(r"[A-Z]", value) is None for value in values):
            raise ValueError("allowed_residues must contain one-letter residue codes")
        if name == "source_protein_ids" and not values <= bundle["proteins"].keys():
            raise ValueError("Unknown source_protein_ids in the search FASTA")
        selections[name] = values

    if type(export_unverified_candidates) is not bool:
        raise ValueError("export_unverified_candidates must be a boolean")

    if export_unverified_candidates:
        if fdr_level is not None or fdr_cutoff is not None:
            raise ValueError(
                "Unverified candidate export cannot be combined with "
                "per-record FDR selection"
            )
        prepared = _prepare_plink_candidate_input(bundle, mappings)
    else:
        prepared = prepare_af_input(
            bundle, mappings,
            fdr_level=fdr_level,
            fdr_cutoff=fdr_cutoff,
        )
    if reference_projection is not None:
        prepared = _prepare_reference_input(bundle, mappings, prepared, reference_projection)
    eligible = {row["source_record_id"]: row for row in prepared["rows"]}
    selected = []
    for record in prepared["audit"]:
        row = eligible.get(record["source_record_id"])
        if row is None:
            continue
        reasons = []
        scope = selections["source_protein_ids"]
        if scope is not None and any(endpoint["source_fasta_id"] not in scope for endpoint in record["endpoints"].values()):
            reasons.append("outside_selected_source_proteins")
        if row["source_crosslinker"] not in selections["allowed_crosslinkers"]:
            reasons.append("crosslinker_not_selected")
        for site in ("a", "b"):
            residue = row[f"seq_{site}"][row[f"pos_{site}"] - 1]
            if residue not in selections["allowed_residues"]:
                reasons.append("residue_not_selected_" + site)
        if reasons:
            record.update(input_status="not_evaluable", input_reasons=reasons)
        else:
            selected.append(row)
    prepared["rows"] = selected
    if export_unverified_candidates:
        selected_ids = {row["source_record_id"] for row in selected}
        for record in prepared["audit"]:
            if record["source_record_id"] in selected_ids:
                record.update(
                    input_status="candidate_ready",
                    input_scope="unverified_candidate",
                    fdr_verification_status="unverified",
                )
    policy = {key: sorted(values) if values is not None else None for key, values in selections.items()}
    return write_xifdr_outputs(bundle, mappings, prepared, output_dir, pipeline_policy=policy)


def write_xifdr_outputs(bundle, mappings, prepared, output_dir, *, pipeline_policy=None):
    """Write source/audit snapshots and the eligible input table to a new directory."""
    output_dir = Path(output_dir)
    _verify_msannika_sources(bundle)
    _verify_plink_native_sources(bundle)
    author_sources = _author_filter_source_bodies(bundle)
    copied = deepcopy(mappings)
    source_bodies = dict(prepared.get("extra_source_bodies", {}))
    for identity in copied.values():
        for source in identity["sources"]:
            raw = Path(source["absolute_path"]).read_bytes()
            if hashlib.sha256(raw).hexdigest() != source["sha256"]:
                raise ValueError("A protein-identity source changed after lookup")
            name = "identity_sources/" + source["sha256"] + Path(source["absolute_path"]).suffix
            source_bodies[name] = raw
            source["original_path"], source["path"] = source.pop("absolute_path"), name
    projection = deepcopy(prepared.get("reference_projection"))
    if projection is not None:
        def copy_sources(value):
            if isinstance(value, list):
                for item in value:
                    copy_sources(item)
            elif isinstance(value, dict):
                if "kind" in value and "sha256" in value and "absolute_path" in value:
                    raw = Path(value["absolute_path"]).read_bytes()
                    if hashlib.sha256(raw).hexdigest() != value["sha256"]:
                        raise ValueError("A reference projection evidence source changed")
                    name = "reference_sources/" + value["sha256"] + Path(value["absolute_path"]).suffix
                    source_bodies[name] = raw
                    value["original_path"], value["path"] = value.pop("absolute_path"), name
                for child in value.values():
                    copy_sources(child)
        copy_sources(projection)
    output_dir.mkdir(parents=True, exist_ok=False)
    if source_bodies:
        for name, raw in source_bodies.items():
            destination = output_dir / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(raw)
    _write_source_bodies(output_dir, author_sources)
    _write_json(output_dir / "psm_audit.json", prepared["audit"])
    _write_json(output_dir / "search_proteins.json", bundle["proteins"])
    _write_json(output_dir / "protein_mapping.json", copied)
    _write_json(output_dir / "source_columns.json", bundle["source_columns"])
    rows = prepared["rows"]
    pipeline = pipeline_policy is not None
    input_name = "claudio_input.csv" if pipeline else "af_input.csv"
    core_columns = ["pep_a", "pep_b", "pos_a", "pos_b", "res_pos_a", "res_pos_b", "unip_id_a", "unip_id_b"]
    if pipeline:
        default_columns = ["peptide1", "peptide2", "position1", "position2",
                           "k_pos1", "k_pos2", "entry1", "entry2"]
        rename = dict(zip(core_columns, default_columns))
        # M1 owns these fields and retrieves current UniProt sequences itself.
        generated = {"source_record_id", "original_pos_a", "original_pos_b", "candidate_id",
                     "seq_a", "seq_b", "position_status_a", "position_status_b"}
        rows = [{"xlms_source_record_id": row["source_record_id"],
                 **{rename.get(key, key): value for key, value in row.items() if key not in generated}}
                for row in rows]
        core_columns = default_columns
        extra = [key for key in rows[0] if key not in core_columns] if rows else ["xlms_source_record_id"]
        columns = core_columns + extra
    else:
        columns = list(rows[0]) if rows else ["source_record_id", "unip_id_a", "unip_id_b", "pos_a", "pos_b"]
    with (output_dir / input_name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    names = ["psm_audit.json", "search_proteins.json", "protein_mapping.json", "source_columns.json", input_name,
             *source_bodies, *author_sources]
    if projection is not None:
        _write_json(output_dir / "reference_projection.json", projection)
        names.append("reference_projection.json")
    format_name = {
        "plink_psm": "plink",
        "merox_csm": "merox",
        "msannika_csm": "msannika",
    }.get(bundle.get("format"), "xifdr")
    if bundle.get("format") == "scout_csm":
        format_name = "scout"
    if bundle.get("format") == "xlinkx_csm":
        format_name = "xlinkx"
    manifest = {"schema_version": 1, "kind": format_name + ("_claudio_input" if pipeline else "_prepared_input"), "status": "complete",
                "provenance": bundle["provenance"], "fdr_policy": prepared["fdr_policy"],
                "counts": {"source_records": len(bundle["records"]),
                           "pipeline_input_rows" if pipeline else "af_input_rows": len(rows),
                           "protein_mapping_statuses": dict(Counter(m["status"] for m in mappings.values()))},
                "assessment_performed": False,
                "output_sha256": {name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest() for name in names}}
    if projection is not None:
        manifest["reference_projection_policy"] = deepcopy(prepared["reference_projection_policy"])
        manifest["counts"]["reference_projected_input_rows"] = sum(
            row["identity_scope"] == "candidate_hypothesis" for row in rows)
    if "parser" in bundle:
        manifest["parser"] = deepcopy(bundle["parser"])
    if "native_filter_policy" in bundle:
        manifest["native_filter_policy"] = deepcopy(bundle["native_filter_policy"])
    if author_sources:
        manifest["author_filter_policy"] = deepcopy(bundle["author_filter_policy"])
        manifest["author_filter_evidence_snapshot"] = "author_filter_sources/" + Path(bundle["provenance"]["author_filter_path"]).name
    if pipeline:
        manifest["selection_policy"] = pipeline_policy
        manifest["pipeline_input"] = {
            "path": input_name, "projections": ",".join(core_columns),
            "source_key": "xlms_source_record_id", "requires_source_preservation": True,
        }
    if "author_filter_policy" in prepared:
        manifest["author_filter_policy"] = deepcopy(
            prepared["author_filter_policy"]
        )
        manifest["author_filter_evidence_snapshot"] = prepared[
            "author_filter_evidence_snapshot"
        ]
    _write_json(output_dir / "manifest.json", manifest)
    return output_dir


def _af_backend():
    try:
        import pandas as pd
        from claudio.module02.src_structure_search.algorithm.af_complex_search import retrieve_discovered_complexes
        from claudio.module02.src_distance_reevaluation.algorithm.af_complex_distances import calculate_complex_site_dists
        from claudio.module04.src.algorithm.distance_assessment import annotate_distance_assessments
        from claudio.module04.src.algorithm.relationship_classification import classify_model_candidates
    except ImportError as error:
        raise RuntimeError("The AF assessment branch is required for model retrieval and distances; the dataset catalog branch can still parse and prepare inputs") from error
    return pd, retrieve_discovered_complexes, calculate_complex_site_dists, classify_model_candidates, annotate_distance_assessments


def _check_af_sources(table, expected):
    if "source_record_id" not in table.columns or set(table["source_record_id"]) != expected:
        raise ValueError("AF analysis lost or introduced source records")


def assess_af_complexes(prepared, output_dir, *, xl_residues, allowed_crosslinkers, plddt_cutoff,
                       topolink_bin, euclidean_policy, topological_policy, topolink_seed=None):
    """Use existing strict M2 mapping and M4 assessments; no duplicate AF implementation.

    This entry point assesses complex candidates only. Chemistry, confidence and
    distance policies are explicit inputs. Excluded PSM rows stay in the audit.
    """
    if not prepared["rows"]:
        raise ValueError("No eligible PSMs are available for AF complex assessment")
    if isinstance(allowed_crosslinkers, str) or not allowed_crosslinkers:
        raise ValueError("An explicit set of allowed crosslinker labels is required")
    observed = {r["source_crosslinker"] for r in prepared["rows"]}
    if not observed <= set(allowed_crosslinkers) or "" in observed:
        raise ValueError("Input crosslinker labels do not match the explicit crosslinker policy")
    if not isinstance(xl_residues, str) or not xl_residues.strip():
        raise ValueError("Explicit residue/atom rules are required")
    for policy in (euclidean_policy, topological_policy):
        if not isinstance(policy, dict) or set(policy) != {"minimum", "maximum", "boundary"}:
            raise ValueError("Distance policies require minimum, maximum and boundary")
        minimum, maximum = policy["minimum"], policy["maximum"]
        if (any(isinstance(v, bool) or not isinstance(v, Real) or not math.isfinite(v) for v in (minimum, maximum))
                or not 0 <= minimum < maximum or policy["boundary"] not in ("inclusive", "exclusive")):
            raise ValueError("Invalid distance policy")
    if isinstance(plddt_cutoff, bool) or not isinstance(plddt_cutoff, Real) or not math.isfinite(plddt_cutoff) or not 0 <= plddt_cutoff <= 100:
        raise ValueError("Invalid pLDDT cutoff")
    pd, retrieve, distances, classify, assess = _af_backend()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    structures, work = output_dir / "structures", output_dir / "work"
    structures.mkdir()
    work.mkdir()
    data = pd.DataFrame(prepared["rows"])
    expected = {row["source_record_id"] for row in prepared["rows"]}
    candidates = retrieve(data, str(structures))
    _check_af_sources(candidates, expected)
    calculated = distances(candidates, str(work), xl_residues, plddt_cutoff, topolink_bin, 0, topolink_seed=topolink_seed)
    _check_af_sources(calculated, expected)
    # M2's low-level calculation returns *_tplk; M4 consumes *_dist.
    calculated["eucl_dist"] = calculated.pop("eucl_dist_tplk")
    calculated["topo_dist"] = calculated.pop("topo_dist_tplk")
    classified = classify(calculated)
    result = assess(classified, euclidean_policy=euclidean_policy,
                    topological_policy=topological_policy, plddt_cutoff=plddt_cutoff)
    _check_af_sources(result, expected)
    result.to_csv(output_dir / "model_candidates.csv", index=False)
    _write_json(output_dir / "psm_input_audit.json", prepared["audit"])
    _write_json(output_dir / "manifest.json", {"schema_version": 1, "kind": "xifdr_af_complex_assessment", "status": "complete",
                "provenance": prepared["provenance"], "fdr_policy": prepared["fdr_policy"],
                "xl_residues": xl_residues, "crosslinkers": sorted(observed), "plddt_cutoff": plddt_cutoff,
                "euclidean_policy": euclidean_policy, "topological_policy": topological_policy,
                "topolink_seed": topolink_seed, "input_record_count": len(data), "model_candidate_count": len(result),
                "output_sha256": {name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
                                  for name in ("model_candidates.csv", "psm_input_audit.json")}})
    return result
def _merox_endpoint(fields, site, proteins):
    keys = {
        "a": (
            "Peptide 1", "Protein 1", "From", "To",
            "best linkage position peptide 1",
        ),
        "b": (
            "Peptide2", "Protein 2", "From.1", "To.1",
            "best linkage position peptide 2",
        ),
    }
    pep_key, protein_key, start_key, end_key, link_key = keys[site]

    identifiers = set(re.findall(
        r"(?:^|>)\s*([^\s@()>]+)", fields[protein_key]
    ))
    if len(identifiers) != 1 or not identifiers <= proteins.keys():
        raise ValueError(
            "Protein annotation does not identify one search FASTA entry"
        )

    identifier = next(iter(identifiers))
    sequence = proteins[identifier]["sequence"]
    raw = fields[pep_key]
    if not raw or raw[0] not in "[{" or raw[-1] not in "]}":
        raise ValueError("Unsupported MeroX peptide boundaries")

    encoded = raw[1:-1]
    peptide = encoded.replace("B", "C").replace("m", "M")
    if (
        not peptide
        or re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", peptide) is None
    ):
        raise ValueError("Unsupported MeroX residue or modification symbol")

    start = int(fields[start_key])
    end = int(fields[end_key])
    if raw.startswith("{"):
        if start != 0:
            raise ValueError("Unexpected MeroX N-terminal coordinate")
        start = 1
    if raw.endswith("}"):
        if end != len(sequence) + 1:
            raise ValueError("Unexpected MeroX C-terminal coordinate")
        end = len(sequence)

    if (
        not 1 <= start <= end <= len(sequence)
        or sequence[start - 1:end] != peptide
    ):
        raise ValueError("MeroX peptide span does not match the search FASTA")

    match = re.fullmatch(r"([A-Za-z{}])(-?\d+)", fields[link_key])
    if match is None:
        raise ValueError("Unsupported MeroX linkage annotation")
    label, position = match[1], int(match[2])

    if label == "{":
        if not raw.startswith("{") or position != 0 or start != 1:
            raise ValueError("Inconsistent MeroX N-terminal site")
        position, kind = 1, "protein_n_terminus"
    elif label == "}":
        if not raw.endswith("}") or position != len(peptide) + 1:
            raise ValueError("Inconsistent MeroX C-terminal site")
        position, kind = len(peptide), "protein_c_terminus"
    else:
        residue = label.replace("B", "C").replace("m", "M")
        if (
            not 1 <= position <= len(peptide)
            or peptide[position - 1] != residue
        ):
            raise ValueError("MeroX linkage residue does not match the peptide")
        kind = "residue"

    fields.update({
        pep_key: "[" + encoded + "]",
        protein_key: ">" + identifier,
        start_key: start,
        end_key: end,
        link_key: label + str(position),
    })
    return {
        "source_protein_id": identifier,
        "source_fasta_id": identifier,
        "source_peptide": raw,
        "peptide": peptide,
        "peptide_start": start,
        "peptide_end": end,
        "link_position": position,
        "protein_position": start + position - 1,
        "kind": kind,
        "status": "valid",
        "reason": "",
    }
def _merox_record_kind(fields):
    """Recognize paired peptides and the observed water dead-end encoding."""
    keys = (
        (
            "Peptide 1", "Protein 1", "From", "To",
            "best linkage position peptide 1",
        ),
        (
            "Peptide2", "Protein 2", "From.1", "To.1",
            "best linkage position peptide 2",
        ),
    )
    sides = [
        tuple(fields.get(key) for key in side_keys)
        for side_keys in keys
    ]
    placeholder = ("0", "H2O", "0", "0", "x0")
    placeholders = [side == placeholder for side in sides]

    peptides = [
        isinstance(side[0], str)
        and len(side[0]) > 2
        and side[0][0] in "[{"
        and side[0][-1] in "]}"
        and isinstance(side[1], str)
        and bool(side[1].strip())
        and side[1] != "H2O"
        for side in sides
    ]

    if any(placeholders):
        if placeholders[0] and peptides[1]:
            return "dead_end"
        if placeholders[1] and peptides[0]:
            return "dead_end"
        return "unresolved"

    if any(
        side[0] == "0" or side[1] == "H2O" or side[4] == "x0"
        for side in sides
    ):
        return "unresolved"

    return "paired_crosslink_candidate" if all(peptides) else "unresolved"
def _merox_localization(fields, proteins):
    """Validate paired native candidates without selecting a localization.

    Candidate order, duplicate entries and uninterpreted native fields are
    preserved. A valid candidate only means its source coordinates agree
    with the search FASTA; it does not establish localization confidence.
    """
    raw = fields.get("All linkage positions", "")
    tokens = raw.split("#") if raw.strip() else []
    best_pair = (
        fields.get("best linkage position peptide 1"),
        fields.get("best linkage position peptide 2"),
    )
    candidates = []
    best_numbers = []
    dead_end = _merox_record_kind(fields) == "dead_end"
    for number, token in enumerate(tokens, 1):
        parts = token.split("~")
        is_source_best = (
            len(parts) == 5
            and tuple(parts[1:3]) == best_pair
        )
        candidate = {
            "source_candidate_number": number,
            "raw": token,
            "source_fields": parts,
            "is_source_best": is_source_best,
            "endpoints": {},
            "status": "invalid",
            "reason": "",
        }
        candidates.append(candidate)

        if is_source_best:
            best_numbers.append(number)
        if dead_end:
            candidate.update(
                status="not_applicable",
                reason="dead_end_not_pairwise_crosslink",
            )
            continue
        if len(parts) != 5:
            candidate["reason"] = "unsupported_candidate_field_count"
            continue

        # Normalize only a private view; retain the original source fields.
        view = dict(fields)
        view["best linkage position peptide 1"] = parts[1]
        view["best linkage position peptide 2"] = parts[2]

        try:
            for site in ("a", "b"):
                candidate["endpoints"][site] = _merox_endpoint(
                    view, site, proteins
                )
        except (ValueError, KeyError, TypeError, IndexError) as error:
            candidate["reason"] = str(error)
            continue

        candidate["status"] = "valid"

    if dead_end:
        status = "not_applicable"
        reason = "dead_end_not_pairwise_crosslink"
    elif not candidates:
        status = "unresolved"
        reason = "no_localization_candidates"
    elif any(item["status"] != "valid" for item in candidates):
        status = "unresolved"
        reason = "invalid_localization_candidate"
    elif not best_numbers:
        status = "unresolved"
        reason = "source_best_pair_not_found"
    else:
        status = (
            "single_candidate"
            if len(candidates) == 1
            else "multiple_candidates"
        )
        reason = ""

    return {
        "raw": raw,
        "status": status,
        "reason": reason,
        "source_best_candidate_numbers": best_numbers,
        "candidates": candidates,
    }
def _merox_confidence_evidence(settings, native_local_fdr):
    """Compare native local FDR with saved filter settings.

    Configuration and numerical agreement do not prove that filtering
    was applied, establish an FDR level, or resolve target/decoy identity.
    """
    cutoff = _fdr_value(settings.get("FDRCUTOFF"))
    local_fdr = _fdr_value(native_local_fdr.get("raw"))

    evidence = {
        "status": "unresolved",
        "reason": "",
        "fdr_level": "native_unspecified",
        "cutoff": cutoff["value"],
        "value": local_fdr["value"],
        "application_verified": False,
    }

    if settings.get("APPLYFDR") == "0":
        evidence.update(
            status="filter_disabled",
            reason="native_fdr_filter_disabled",
        )
    elif settings.get("APPLYFDR") != "1":
        evidence["reason"] = "native_filter_flag_unresolved"
    elif settings.get("PERFORMDECOY") != "1":
        evidence["reason"] = "native_decoy_analysis_not_enabled"
    elif cutoff["status"] != "valid":
        evidence["reason"] = "native_fdr_cutoff_unresolved"
    elif (
        native_local_fdr.get("status") != "valid"
        or local_fdr["status"] != "valid"
    ):
        evidence["reason"] = "native_local_fdr_unresolved"
    elif local_fdr["value"] <= cutoff["value"]:
        evidence["status"] = "within_saved_cutoff"
    else:
        evidence.update(
            status="above_saved_cutoff",
            reason="native_local_fdr_above_saved_cutoff",
        )

    return evidence
def _merox_native_settings(raw):
    """Read saved settings and the explicitly selected crosslinker block."""
    global_keys = {
        "FDRCUTOFF", "APPLYFDR", "PERFORMDECOY", "DECOYTYPE",
        "SAVEDECOYS", "DECOYIDENTIFIER", "USEDCROSSLINKER",
        "NAME", "SITE1", "SITE2",
    }
    linker_keys = {"NAME", "SITE1", "SITE2"}
    settings, blocks = {}, []
    current = None

    def store(destination, key, value):
        if key in destination:
            raise ValueError("Repeated native MeroX setting: " + key)
        destination[key] = value

    for original in raw.splitlines():
        line = original.split(b"//", 1)[0].strip()
        if not line:
            continue

        key_raw, separator, value_raw = line.partition(b"=")
        key = key_raw.strip().decode("ascii", errors="backslashreplace")
        value = value_raw.strip().decode(
            "ascii", errors="backslashreplace"
        )

        if key == "CROSSLINKER":
            current = {}
            if separator:
                if not value:
                    raise ValueError("Empty native crosslinker name")
                current["NAME"] = value
            blocks.append(current)
            continue

        if line == b"END":
            current = None
            continue

        if not separator:
            continue

        if current is not None:
            if key in linker_keys:
                store(current, key, value)
        elif key in global_keys:
            store(settings, key, value)

    selected = settings.get("USEDCROSSLINKER")
    if blocks:
        names = [block.get("NAME") for block in blocks]
        if any(not name for name in names):
            raise ValueError("Native crosslinker block has no name")
        if len(set(names)) != len(names):
            raise ValueError("Repeated native crosslinker block name")

        if selected is None:
            if len(blocks) != 1:
                raise ValueError(
                    "Multiple native crosslinkers without an explicit selection"
                )
            active = blocks[0]
        else:
            matches = [
                block for block in blocks
                if block["NAME"] == selected
            ]
            if len(matches) != 1:
                raise ValueError(
                    "Selected native crosslinker block was not found"
                )
            active = matches[0]

        settings.update(active)
    elif selected is not None:
        raise ValueError("Selected native crosslinker block was not found")

    return settings


def _merox_spectrum_identity(source_identifier, embedded_identifier):
    """Validate spectrum identity before preparing a pyXLMS scan reference."""
    if re.fullmatch(r"[0-9]+", source_identifier):
        scan = _positive_integer(source_identifier)
        if _positive_integer(embedded_identifier) != scan:
            raise ValueError(
                "Embedded spectrum scan differs from its source row"
            )
        return {
            "scan_number": scan,
            "original_file": None,
            "native_id": None,
            "source_identifier": source_identifier,
        }

    if source_identifier != embedded_identifier:
        raise ValueError(
            "Embedded spectrum identifier differs from its source row"
        )

    label_match = re.fullmatch(
        r"([0-9]+)~([^\r\n]+)", source_identifier
    )
    if label_match is not None:
        label = label_match[2]
        if not label.strip():
            raise ValueError("Native MeroX source label is empty")
        return {
            "scan_number": _positive_integer(label_match[1]),
            "original_file": None,
            "native_id": None,
            "source_identifier": source_identifier,
            "source_file_label": label,
        }

    match = re.fullmatch(
        r'File:"([^"]+)",NativeID:"(.+)"',
        source_identifier,
    )
    if match is None:
        raise ValueError("Unsupported native MeroX spectrum identifier")

    original_file, native_id = match.groups()
    scan_match = re.match(r"scan=([0-9]+)(?:,|$)", native_id)
    if scan_match is None:
        raise ValueError("Native spectrum identifier has no supported scan")

    return {
        "scan_number": _positive_integer(scan_match[1]),
        "original_file": original_file,
        "native_id": native_id,
        "source_identifier": source_identifier,
    }
def read_merox_csms(psm_path, fasta_path, *, dataset_accession):
    """Audit supported 27/33-column MeroX native results.

    Terminal positions are typed residue anchors. Native localization strings,
    confidence values and source fields remain unmodified in the audit.
    """
    from importlib.metadata import version

    try:
        from pyXLMS import parser
        from pyXLMS.parser._parser_xldbse_merox import MEROX_COLNAMES
    except ImportError as error:
        raise RuntimeError("Reading MeroX requires pyXLMS==2.0.7") from error

    if (
        not isinstance(dataset_accession, str)
        or re.fullmatch(r"PXD\d{6,}", dataset_accession) is None
    ):
        raise ValueError("A ProteomeXchange PXD accession is required")

    columns = list(MEROX_COLNAMES)
    if (
        len(columns) < 33
        or columns[14] != "Scan number"
        or columns[23] != "Spectrum UUID"
    ):
        raise RuntimeError("Unsupported pyXLMS MeroX column layout")

    psm_path, fasta_path = Path(psm_path), Path(fasta_path)
    raw = psm_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    proteins = read_search_fasta(fasta_path)
    records, prepared, spectra, settings = [], {}, {}, {}

    marker = "__claudio_merox_source_id"
    stream = io.StringIO()
    writer = csv.DictWriter(
        stream, fieldnames=columns + [marker], delimiter=";"
    )
    writer.writeheader()

    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        members = Counter(archive.namelist())
        if any(
            members[name] != 1
            for name in ("Result.csv", "properties.mxf")
        ):
            raise ValueError(
                "MeroX must contain one Result.csv and one properties.mxf"
            )

        properties = archive.read("properties.mxf")
        settings = _merox_native_settings(properties)
        result_text = archive.read("Result.csv").decode("utf-8-sig")

        widths = {
            len(row)
            for row in csv.reader(
                io.StringIO(result_text), delimiter=";"
            )
            if any(value.strip() for value in row)
        }
        source_width = 33 if 33 in widths else 27
        reader = csv.reader(io.StringIO(result_text), delimiter=";")
        previous_line = 0

        for number, values in enumerate(reader, 1):
            line_start, line_end = previous_line + 1, reader.line_num
            previous_line = line_end
            if not any(value.strip() for value in values):
                continue

            fields = (
                dict(zip(columns[:source_width], values))
                if len(values) == source_width
                else {}
            )
            record = {
                "source_record_id": (
                    f"{dataset_accession}:{digest}:{number:06d}"
                ),
                "source_record_number": number,
                "source_line_start": line_start,
                "source_line_end": line_end,
                "source_values": values,
                "source_fields": fields,
                "endpoints": {},
                "spectrum": None,
                "decoy_status": "unknown",
                "decoy_reason": "native_decoy_policy_not_reviewed",
                "fdr": {
                    level: _fdr_value(None) for level in FDR_COLUMNS
                },
                "record_status": "not_evaluable",
                "record_reason": "merox_assessment_policies_pending",
                "record_kind": _merox_record_kind(fields),
            }
            records.append(record)
            record["localization"] = _merox_localization(
                fields, proteins
            )
            if len(values) != source_width:
                record["record_reason"] = "unsupported_merox_row_width"
                continue

            record["native_local_fdr"] = _fdr_value(
                fields["local False discovery rate"]
            )
            record["native_confidence"] = _merox_confidence_evidence(
                settings, record["native_local_fdr"]
            )
            if record["record_kind"] == "unresolved":
                record["record_reason"] = "merox_record_kind_unresolved"
                continue
            view = fields.copy()
            native_decoy_sides = []
            if record["record_kind"] == "paired_crosslink_candidate":
                native_decoy_sides = [
                    site
                    for site, key in (
                        ("a", "Protein 1"),
                        ("b", "Protein 2"),
                    )
                    if fields[key].startswith("DEC_>")
                ]

            if native_decoy_sides:
                record.update(
                    decoy_status=(
                        "decoy_decoy"
                        if len(native_decoy_sides) == 2
                        else "target_decoy"
                    ),
                    decoy_reason="native_decoy_annotation",
                    decoy_evidence="native_merox_DEC_prefix",
                    native_decoy_sites=native_decoy_sides,
                )

            try:
                if record["record_kind"] == "dead_end":
                    peptide_sites = (
                        ("b",)
                        if fields["Peptide 1"] == "0"
                        else ("a",)
                    )
                else:
                    peptide_sites = ("a", "b")

                if not native_decoy_sides:
                    for site in peptide_sites:
                        record["endpoints"][site] = _merox_endpoint(
                            view, site, proteins
                        )

                uuid = fields["Spectrum UUID"]
                if re.fullmatch(r"[A-Za-z0-9_-]+", uuid) is None:
                    raise ValueError("Invalid spectrum UUID")

                member = f"spectra/{uuid}.txt"
                if members[member] != 1:
                    raise ValueError(
                        "Missing or repeated embedded spectrum"
                    )

                if member not in spectra:
                    body = archive.read(member)
                    header = body.splitlines()[0].decode("utf-8-sig")
                    spectra[member] = (
                        header,
                        hashlib.sha256(body).hexdigest(),
                    )

                embedded_identifier, member_hash = spectra[member]
                spectrum_identity = _merox_spectrum_identity(
                    fields["Scan number"], embedded_identifier
                )
                scan = spectrum_identity["scan_number"]
            except (
                ValueError, IndexError, KeyError, zipfile.BadZipFile
            ) as error:
                record.update(
                    record_reason="merox_source_unresolved",
                    source_error=str(error),
                )
                continue

            reference = f"merox-embedded:sha256:{digest}!/{member}"
            record["spectrum"] = {
                "uuid": uuid,
                **spectrum_identity,
                "reference_kind": "merox_embedded_member",
                "reference": reference,
                "member": member,
                "member_sha256": member_hash,
            }
            if native_decoy_sides:
                record["record_reason"] = "decoy_record"
                continue
            if record["record_kind"] == "dead_end":
                record["record_reason"] = "dead_end_not_pairwise_crosslink"
                continue
            view.update({
                "Scan number": f"{scan}~{reference}",
                marker: record["source_record_id"],
            })
            writer.writerow(view)
            prepared[record["source_record_id"]] = record

    parser_info = {
        "name": "pyXLMS",
        "version": version("pyXLMS"),
        "entrypoint": "read_merox",
        "format_variant": (
            "legacy_27_column_native"
            if source_width == 27
            else "native_33_column"
        ),
        "parsed_records": 0,
        "modification_masses_assessed": False,
    }

    if prepared and not settings.get("NAME"):
        for record in prepared.values():
            record["record_reason"] = "merox_crosslinker_not_declared"
    elif prepared:
        try:
            result = parser.read_merox(
                io.BytesIO(stream.getvalue().encode("utf-8")),
                crosslinker=settings["NAME"],
                parse_modifications=False,
                sep=";",
            )
            by_source = {}
            for csm in result["crosslink-spectrum-matches"] or []:
                source = (csm["additional_information"] or {}).get(
                    "source", {}
                )
                source_id = source.get(marker)
                if source_id not in prepared or source_id in by_source:
                    raise ValueError(
                        "Unknown or repeated MeroX parser source ID"
                    )
                by_source[source_id] = csm

            if set(by_source) != set(prepared):
                raise ValueError("MeroX parser omitted source records")

            for source_id, record in prepared.items():
                csm = by_source[source_id]
                expected = [
                    (
                        end["peptide"],
                        end["link_position"],
                        (end["source_fasta_id"],),
                        (end["protein_position"],),
                        (end["peptide_start"],),
                    )
                    for end in (
                        record["endpoints"]["a"],
                        record["endpoints"]["b"],
                    )
                ]
                actual = [
                    (
                        csm[f"{side}_peptide"],
                        csm[f"{side}_peptide_crosslink_position"],
                        tuple(csm[f"{side}_proteins"] or []),
                        tuple(
                            csm[f"{side}_proteins_crosslink_positions"] or []
                        ),
                        tuple(
                            csm[f"{side}_proteins_peptide_positions"] or []
                        ),
                    )
                    for side in ("alpha", "beta")
                ]
                spectrum = record["spectrum"]
                if (
                    (actual != expected and actual != expected[::-1])
                    or csm["scan_nr"] != spectrum["scan_number"]
                    or csm["spectrum_file"] != spectrum["reference"]
                ):
                    raise ValueError(
                        "MeroX parser changed endpoints or spectrum identity"
                    )

                record["parser_orientation"] = (
                    "same" if actual == expected else "swapped"
                )

            parser_info["parsed_records"] = len(by_source)
        except (
            ValueError, TypeError, KeyError, RuntimeError,
            IndexError, AttributeError, OverflowError,
        ) as error:
            for record in prepared.values():
                record.pop("parser_orientation", None)
                record.update(
                    record_reason="pyxlms_parse_failed",
                    parser_error=f"{type(error).__name__}: {error}",
                )

    return {
        "format": "merox_csm",
        "format_version": 1,
        "source_columns": columns[:source_width],
        "records": records,
        "proteins": proteins,
        "parser": parser_info,
        "native_settings": settings,
        "scientific_review": "pending",
        "provenance": {
            "dataset_accession": dataset_accession,
            "psm_path": str(psm_path.resolve()),
            "psm_sha256": digest,
            "fasta_path": str(fasta_path.resolve()),
            "fasta_sha256": next(
                iter(proteins.values())
            )["fasta_sha256"],
            "native_properties_sha256": hashlib.sha256(
                properties
            ).hexdigest(),
        },
    }
def _verify_merox_construct_review(source, specification):
    """Verify a source-bound declaration for a construct reference."""
    review = specification.get("construct_review")
    if not isinstance(review, dict):
        raise ValueError("Construct alignment requires an explicit review")

    reference = specification["reference_sequence"]
    expected = {
        "kind": "merox_construct_reference_review",
        "schema_version": 1,
        "review_status": "reviewed",
        "source_fasta_id": source["id"],
        "source_fasta_sha256": source["fasta_sha256"],
        "source_sequence_sha256": hashlib.sha256(
            source["sequence"].encode("ascii")
        ).hexdigest(),
        "reference_accession": specification["reference_accession"],
        "reference_taxon_id": specification["reference_taxon_id"],
        "reference_sequence_sha256": hashlib.sha256(
            reference.encode("ascii")
        ).hexdigest(),
    }
    if (
        type(review.get("schema_version")) is not int
        or any(review.get(key) != value for key, value in expected.items())
        or not isinstance(review.get("reason"), str)
        or not review["reason"].strip()
    ):
        raise ValueError("Construct review differs from its source or reference")

    interval = review.get("reference_interval")
    maximum = review.get("max_optimal_alignments")
    if (
        not isinstance(interval, (list, tuple))
        or len(interval) != 2
        or any(type(value) is not int for value in interval)
        or not 1 <= interval[0] <= interval[1] <= len(reference)
        or type(maximum) is not int
        or not 1 <= maximum <= 64
        or review.get("require_all_optimal_alignments") is not True
    ):
        raise ValueError("Invalid construct review projection policy")

    receipts = review.get("sources")
    if not isinstance(receipts, list) or not receipts:
        raise ValueError("Construct review requires saved evidence")
    declared = {
        key: value for key, value in review.items() if key != "sources"
    }
    for receipt in receipts:
        if (
            not isinstance(receipt, dict)
            or receipt.get("kind") != "merox_construct_reference_review"
            or not isinstance(receipt.get("absolute_path"), str)
            or not Path(receipt["absolute_path"]).is_absolute()
            or not isinstance(receipt.get("sha256"), str)
            or re.fullmatch(r"[0-9a-f]{64}", receipt["sha256"]) is None
        ):
            raise ValueError("Invalid construct review evidence")
        raw = Path(receipt["absolute_path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != receipt["sha256"]:
            raise ValueError("Construct review evidence checksum changed")
        if json.loads(raw.decode("utf-8-sig")) != declared:
            raise ValueError("Construct review document differs from its declaration")
    return review


def _merox_multi_reference_profiles(bundle, reference_mapping):
    """Prepare independently verified references for each source protein."""
    if (
        type(reference_mapping.get("schema_version")) is not int
        or reference_mapping["schema_version"] != 1
        or reference_mapping.get("source_fasta_sha256")
        != bundle["provenance"]["fasta_sha256"]
    ):
        raise ValueError("Multi-reference mapping differs from the search FASTA")

    references = reference_mapping.get("references")
    if not isinstance(references, dict) or not references:
        raise ValueError("Multi-reference mapping requires protein references")

    profiles, options_by_source, mappings = {}, {}, {}
    reference_sequences = {}
    for identifier, specification in references.items():
        source = bundle["proteins"].get(identifier)
        if source is None or not isinstance(specification, dict):
            raise ValueError("Unknown source protein in multi-reference mapping")

        accession = specification.get("reference_accession")
        taxon = specification.get("reference_taxon_id")
        sequence = specification.get("reference_sequence")
        receipts = specification.get("sources")
        if (
            not isinstance(accession, str)
            or ACCESSION.fullmatch(accession) is None
            or type(taxon) is not int or taxon <= 0
            or not isinstance(sequence, str)
            or re.fullmatch(r"[A-Z]+", sequence) is None
            or not isinstance(receipts, list) or not receipts
        ):
            raise ValueError("Invalid per-protein reference declaration")
        if accession in reference_sequences and reference_sequences[accession] != sequence:
            raise ValueError("Conflicting sequences for the same reference accession")
        reference_sequences[accession] = sequence

        profile = {
            "profile_id": identifier + ":" + accession,
            "source_fasta_id": identifier,
            "accession": accession,
            "taxon_id": taxon,
            "source_sequence": source["sequence"],
            "reference_sequence": sequence,
            "reference_sequence_sha256": hashlib.sha256(
                sequence.encode("ascii")
            ).hexdigest(),
        }
        _reference_evidence(
            profile,
            {"taxon_id": taxon, "sources": receipts},
            source,
            require_declared_accession=False,
        )

        review = None
        maximum = 64
        if specification.get("construct_review") is not None:
            review = _verify_merox_construct_review(source, specification)
            maximum = review["max_optimal_alignments"]

        options, coordinate_evidence = _projection_coordinate_options(
            source["sequence"], sequence, maximum
        )
        if review is None:
            if (
                coordinate_evidence["method"] != "full_source_exact_subsequence"
                or len(options) != 1
            ):
                raise ValueError(
                    "Construct alignment requires a reviewed reference mapping"
                )
        elif not options:
            raise ValueError("Construct review has no bounded coordinate mapping")

        exact = (
            review is None
            and source["sequence"] == sequence
            and source.get("uniprot_accession") == accession
        )
        scope = "exact_sequence" if exact else "candidate_hypothesis"
        profile.update(
            identity_scope=scope,
            coordinate_evidence=coordinate_evidence,
            require_all_optimal_alignments=True,
            reference_interval=(
                list(review["reference_interval"]) if review else None
            ),
        )
        profiles[identifier] = profile
        options_by_source[identifier] = options
        sources = deepcopy(receipts)
        if review:
            sources.extend(deepcopy(review["sources"]))
        mappings[identifier] = {
            "status": "mapped" if exact else "reference_projected",
            "accession": accession,
            "sequence": sequence,
            "taxon_id": taxon,
            "identity_scope": scope,
            "coordinate_evidence": deepcopy(coordinate_evidence),
            "sources": sources,
        }
    return profiles, options_by_source, mappings


def write_merox_claudio_input(
    bundle,
    reference_mapping,
    output_dir,
    *,
    allowed_crosslinkers,
    allowed_residues,
):
    """Export source-reported best residue sites using a verified reference.

    Input readiness means compatibility with CLAUDIO's input format.
    Source decoy status, localization uncertainty and FDR scope remain
    unchanged. Terminal sites require a separate atom-selection policy.
    """
    if bundle.get("format") != "merox_csm":
        raise ValueError("Expected a MeroX import bundle")

    selections = {}
    for name, values in (
        ("allowed_crosslinkers", allowed_crosslinkers),
        ("allowed_residues", allowed_residues),
    ):
        if (
            not isinstance(values, (set, frozenset, list, tuple))
            or not values
            or any(
                not isinstance(value, str) or not value.strip()
                for value in values
            )
        ):
            raise ValueError(f"{name} requires explicit nonempty labels")
        selections[name] = set(values)

    if any(
        re.fullmatch(r"[A-Z]", residue) is None
        for residue in selections["allowed_residues"]
    ):
        raise ValueError("Residue labels must be single uppercase letters")

    for kind in ("psm", "fasta"):
        path = Path(bundle["provenance"][kind + "_path"])
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != bundle["provenance"][kind + "_sha256"]:
            raise ValueError(f"Source {kind} checksum changed")

    multi_reference = (
        reference_mapping.get("kind") == "merox_multi_reference_mapping"
    )
    if multi_reference:
        profiles, options_by_source, mappings = _merox_multi_reference_profiles(
            bundle, reference_mapping
        )
    else:
        identifier = reference_mapping.get("source_fasta_id")
        source = bundle["proteins"].get(identifier)
        if (
            source is None
            or reference_mapping.get("source_fasta_sha256")
            != bundle["provenance"]["fasta_sha256"]
        ):
            raise ValueError("Reference mapping differs from the search FASTA")

        accession = reference_mapping.get("reference_accession")
        taxon = reference_mapping.get("reference_taxon_id")
        reference_sequence = reference_mapping.get("reference_sequence")
        if (
            not isinstance(accession, str)
            or ACCESSION.fullmatch(accession) is None
            or type(taxon) is not int
            or taxon <= 0
        ):
            raise ValueError("Invalid reference accession or organism")
        if (
            source.get("taxon_id") is not None
            and source["taxon_id"] != taxon
        ):
            raise ValueError("Reference organism differs from the search FASTA")

        receipts = reference_mapping.get("sources", [])
        confirmed = False
        for receipt in receipts:
            raw = Path(receipt["absolute_path"]).read_bytes()
            if hashlib.sha256(raw).hexdigest() != receipt.get("sha256"):
                raise ValueError("Reference evidence checksum changed")

            uri = urlsplit(
                receipt.get("response_url", receipt.get("request_url", ""))
            )
            if (
                receipt.get("kind") != "uniprot_api"
                or receipt.get("status") != 200
                or uri.scheme != "https"
                or uri.netloc != "rest.uniprot.org"
                or not uri.path.startswith("/uniprotkb/")
            ):
                continue

            payload = json.loads(raw)
            if not isinstance(payload, dict):
                continue
            if (
                payload.get("primaryAccession") == accession
                and payload.get("organism", {}).get("taxonId") == taxon
                and payload.get("sequence", {}).get("value")
                == reference_sequence
            ):
                confirmed = True

        if not confirmed:
            raise ValueError("Reference sequence lacks matching saved evidence")

        options, coordinate_evidence = _projection_coordinate_options(
            source["sequence"], reference_sequence, 64
        )
        if (
            coordinate_evidence["method"] != "full_source_exact_subsequence"
            or len(options) != 1
        ):
            raise ValueError(
                "MeroX export currently requires unique exact source containment"
            )

        profile = {
            "profile_id": identifier + ":" + accession,
            "accession": accession,
            "identity_scope": "candidate_hypothesis",
            "source_sequence": source["sequence"],
            "reference_sequence": reference_sequence,
            "coordinate_evidence": coordinate_evidence,
        }
        mappings = {
            identifier: {
                "status": "reference_projected",
                "accession": accession,
                "sequence": reference_sequence,
                "taxon_id": taxon,
                "identity_scope": "candidate_hypothesis",
                "coordinate_evidence": deepcopy(coordinate_evidence),
                "sources": deepcopy(receipts),
            }
        }

        profiles = {identifier: profile}
        options_by_source = {identifier: options}

    settings = bundle["native_settings"]
    crosslinker = settings.get("NAME", "")
    cutoff = _fdr_value(settings.get("FDRCUTOFF"))
    rows, audit = [], []

    for record in bundle["records"]:
        entry = deepcopy(record)
        if record.get("record_kind") == "dead_end":
            entry.update(
                input_status="not_evaluable",
                input_reasons=["dead_end_not_pairwise_crosslink"],
            )
            audit.append(entry)
            continue
        reasons = []
        fields = record["source_fields"]
        localization = record.get("localization", {})
        best = [
            candidate
            for candidate in localization.get("candidates", [])
            if candidate.get("is_source_best")
        ]
        confidence = _merox_confidence_evidence(
            settings,
            _fdr_value(fields.get("local False discovery rate")),
        )

        if record.get("parser_orientation") not in ("same", "swapped"):
            reasons.append("merox_parser_unresolved")
        if record.get("decoy_status") not in ("unknown", "target_target"):
            reasons.append("decoy_record")
        if confidence["status"] != "within_saved_cutoff":
            reasons.append("native_confidence_not_within_saved_cutoff")
        if crosslinker not in selections["allowed_crosslinkers"]:
            reasons.append("crosslinker_not_selected")
        if (
            localization.get("status")
            not in ("single_candidate", "multiple_candidates")
            or len(best) != 1
            or best[0].get("status") != "valid"
        ):
            reasons.append("source_best_localization_unresolved")

        endpoints = best[0].get("endpoints", {}) if len(best) == 1 else {}
        if set(endpoints) != {"a", "b"}:
            reasons.append("source_endpoint_unresolved")

        for site, endpoint in endpoints.items():
            if endpoint.get("source_fasta_id") not in profiles:
                reasons.append("reference_mapping_unavailable_" + site)
            if endpoint.get("kind") != "residue":
                reasons.append("terminal_atom_policy_unresolved")

        projected = {}
        if not reasons:
            for site in ("a", "b"):
                endpoint = endpoints[site]
                identifier = endpoint["source_fasta_id"]
                site_profile = profiles[identifier]
                source_options = options_by_source[identifier]
                result = _project_peptide(
                    endpoint, site_profile, source_options
                )
                projected[site] = result
                if result["status"] != "projected":
                    reasons.append("reference_projection_unresolved_" + site)
                    continue
                if (
                    site_profile.get("require_all_optimal_alignments", False)
                    and result["compatible_alignment_count"] != len(source_options)
                ):
                    reasons.append(
                        "reference_projection_not_supported_by_all_alignments_" + site
                    )
                    continue
                interval = site_profile.get("reference_interval")
                if interval is not None and not (
                    interval[0] <= result["reference_peptide_start"]
                    and result["reference_peptide_end"] <= interval[1]
                ):
                    reasons.append(
                        "reference_peptide_outside_reviewed_core_" + site
                    )
                    continue
                residue = site_profile["reference_sequence"][
                    result["reference_protein_position"] - 1
                ]
                if residue not in selections["allowed_residues"]:
                    reasons.append("residue_not_selected_" + site)

        entry.update(
            input_status="not_evaluable" if reasons else "ready",
            input_reasons=list(dict.fromkeys(reasons)),
            reference_endpoints=deepcopy(projected),
        )
        audit.append(entry)
        if reasons:
            continue

        row = {
            "source_record_id": record["source_record_id"],
            "dataset_accession": bundle["provenance"]["dataset_accession"],
            "source_record_number": record["source_record_number"],
            "source_psm_sha256": bundle["provenance"]["psm_sha256"],
            "source_fasta_sha256": bundle["provenance"]["fasta_sha256"],
            "source_crosslinker": crosslinker,
            "source_decoy_status": record["decoy_status"],
            "localization_scope": "source_reported_best",
            "localization_status": localization["status"],
            "localization_candidate_count": len(localization["candidates"]),
            "source_best_candidate_number": best[0]["source_candidate_number"],
            "identity_scope": (
                "exact_sequence"
                if all(
                    profiles[endpoint["source_fasta_id"]]["identity_scope"]
                    == "exact_sequence"
                    for endpoint in endpoints.values()
                )
                else "candidate_hypothesis"
            ),
            "native_local_fdr_value": confidence["value"],
            "native_fdr_cutoff": confidence["cutoff"],
            "native_fdr_scope": "native_unspecified",
        }
        for site in ("a", "b"):
            endpoint = endpoints[site]
            site_profile = profiles[endpoint["source_fasta_id"]]
            projection = projected[site]
            row.update({
                f"unip_id_{site}": site_profile["accession"],
                f"seq_{site}": site_profile["reference_sequence"],
                f"pep_{site}": endpoint["peptide"],
                f"pos_{site}": projection["reference_protein_position"],
                f"res_pos_{site}": endpoint["link_position"],
                f"source_protein_id_{site}": endpoint["source_protein_id"],
                f"source_protein_position_{site}": endpoint["protein_position"],
                f"source_peptide_start_{site}": endpoint["peptide_start"],
                f"reference_peptide_start_{site}": (
                    projection["reference_peptide_start"]
                ),
                f"source_modified_peptide_{site}": endpoint["source_peptide"],
            })
        rows.append(row)

    prepared = {
        "rows": rows,
        "audit": audit,
        "fdr_policy": {
            "kind": "native_merox_saved_filter",
            "level": "native_unspecified",
            "cutoff": cutoff["value"],
            "application_verified": False,
            "target_decoy_validation_performed": False,
        },
    }
    policy = {
        "allowed_crosslinkers": sorted(selections["allowed_crosslinkers"]),
        "allowed_residues": sorted(selections["allowed_residues"]),
        "localization": "source_reported_best",
        "terminal_policy": "exclude_pending_atom_configuration",
        "decoy_policy": "preserve_source_status_without_inference",
    }
    if multi_reference:
        prepared["reference_projection"] = deepcopy(reference_mapping)
        prepared["reference_projection_policy"] = {
            "scope": "per_source_protein_reference",
            "intact_peptide_required": True,
            "all_optimal_alignments_required": True,
            "identity_acceptance_performed": False,
        }

    return write_xifdr_outputs(
        bundle, mappings, prepared, output_dir, pipeline_policy=policy
    )
def _scout_endpoint(fields, site, proteins):
    """Find FASTA candidates using the complete peptide and reported position."""
    label = {"a": "Alpha", "b": "Beta"}[site]
    raw_peptide = fields[f"{label} peptide"]
    peptide = raw_peptide.strip()
    raw_mapping = fields[f"{label} protein mapping(s)"]

    result = {
        "source_protein_id": raw_mapping,
        "source_fasta_id": None,
        "source_fasta_candidates": [],
        "source_mapping_scope": "sequence_and_reported_position_candidate",
        "source_peptide": raw_peptide,
        "source_modifications": fields.get(f"{label} modification(s)", ""),
        "peptide": peptide,
        "link_position": None,
        "protein_position": None,
        "peptide_start": None,
        "peptide_end": None,
        "kind": "residue",
        "status": "not_evaluable",
        "reason": "",
    }

    try:
        if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", peptide) is None:
            raise ValueError("Unsupported Scout peptide sequence")

        link = _positive_integer(fields[f"{label} peptide position"])
        result["link_position"] = link
        if link > len(peptide):
            raise ValueError("Scout linkage position exceeds peptide length")

        mappings = [value.strip() for value in raw_mapping.split(";")]
        positions = fields[f"{label} protein(s) position(s)"].split(";")
        if len(mappings) != 1 or len(positions) != 1:
            result["reason"] = "multiple_reported_mappings_require_review"
            return result
        if not mappings[0]:
            raise ValueError("Missing Scout protein mapping")

        position = _positive_integer(positions[0])
        start = position - link + 1
        end = start + len(peptide) - 1
        result.update(
            source_protein_id=mappings[0],
            protein_position=position,
            peptide_start=start,
            peptide_end=end,
        )
        if start < 1:
            raise ValueError("Scout peptide starts before the protein sequence")
    except (ValueError, TypeError, KeyError) as error:
        result["reason"] = str(error)
        return result

    matches = sorted(
        identifier
        for identifier, protein in proteins.items()
        if protein["sequence"][start - 1:end] == peptide
    )
    result["source_fasta_candidates"] = matches
    if len(matches) == 1:
        result.update(source_fasta_id=matches[0], status="valid")
    elif matches:
        result.update(
            status="ambiguous",
            reason="multiple_fasta_entries_match_sequence_and_position",
        )
    else:
        result.update(
            status="unmapped",
            reason="no_fasta_sequence_and_position_match",
        )
    return result


def read_scout_csms(psm_path, fasta_path, *, dataset_accession, crosslinker):
    """Import Scout CSM exports while preserving all original source rows.

    FASTA matches validate source coordinates, not UniProt identity.
    Loop-links remain separate single-peptide records. Source scores and
    localization percentages are never converted into FDR values.
    """
    from importlib.metadata import version

    try:
        from pyXLMS import parser
    except ImportError as error:
        raise RuntimeError("Reading Scout requires pyXLMS==2.0.7") from error

    if (
        not isinstance(dataset_accession, str)
        or re.fullmatch(r"PXD\d{6,}", dataset_accession) is None
    ):
        raise ValueError("A ProteomeXchange PXD accession is required")
    if not isinstance(crosslinker, str) or not crosslinker.strip():
        raise ValueError("An explicit crosslinker label is required")

    psm_path, fasta_path = Path(psm_path), Path(fasta_path)
    raw = psm_path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    proteins = read_search_fasta(fasta_path)
    if not proteins:
        raise ValueError("Search FASTA has no proteins")

    reader = csv.reader(io.StringIO(raw.decode("utf-8-sig")))
    columns = next(reader, None)
    required = {
        "Scan", "Link-Type", "Alpha peptide", "Beta peptide",
        "Alpha peptide position", "Beta peptide position",
        "Alpha protein(s) position(s)", "Beta protein(s) position(s)",
        "Alpha protein mapping(s)", "Beta protein mapping(s)",
        "Precursor charge", "IsDecoy", "Score", "File",
    }
    if not columns or len(columns) != len(set(columns)):
        raise ValueError("Scout CSV requires a nonempty, unique header")
    if not required <= set(columns):
        raise ValueError("Unsupported Scout CSM export columns")

    marker = "__claudio_scout_source_id"
    if marker in columns:
        raise ValueError("Reserved Scout source marker already exists")
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=columns + [marker])
    writer.writeheader()

    records, prepared = [], {}
    previous_line = reader.line_num
    for number, values in enumerate(reader, 1):
        line_start, line_end = previous_line + 1, reader.line_num
        previous_line = line_end
        if not any(value.strip() for value in values):
            continue

        fields = (
            dict(zip(columns, values))
            if len(values) == len(columns)
            else {}
        )
        link_type = fields.get("Link-Type", "").strip()
        kind = (
            "loop_link" if link_type == "Looplink"
            else "paired_crosslink_candidate"
            if link_type in ("Intralink", "Interlink", "Homomeric-inter")
            else "unresolved"
        )
        source_decoy = _boolean(fields.get("IsDecoy"))
        record = {
            "source_record_id": f"{dataset_accession}:{digest}:{number:06d}",
            "source_record_number": number,
            "source_line_start": line_start,
            "source_line_end": line_end,
            "source_values": values,
            "source_fields": fields,
            "record_kind": kind,
            "endpoints": {},
            "spectrum": None,
            "source_decoy_flag": source_decoy,
            "decoy_status": (
                "target_target" if source_decoy is False else "unknown"
            ),
            "decoy_reason": (
                "" if source_decoy is False
                else "source_decoy_endpoint_assignment_unavailable"
                if source_decoy is True
                else "invalid_source_decoy_flag"
            ),
            "fdr": {level: _fdr_value(None) for level in FDR_COLUMNS},
            "record_status": "not_evaluable",
            "record_reason": "scout_source_unresolved",
        }
        records.append(record)

        if len(values) != len(columns):
            record["record_reason"] = "unsupported_scout_row_width"
            continue

        try:
            scan = _positive_integer(fields["Scan"])
            charge = _positive_integer(fields["Precursor charge"])
            if not fields["File"].strip():
                raise ValueError("Missing Scout spectrum filename")
            record["spectrum"] = {
                "scan_number": scan,
                "original_file": fields["File"],
                "charge": charge,
            }

            if kind == "unresolved":
                record["record_reason"] = "unsupported_scout_record_type"
                continue

            sites = ("a",) if kind == "loop_link" else ("a", "b")
            for site in sites:
                record["endpoints"][site] = _scout_endpoint(
                    fields, site, proteins
                )
        except (ValueError, TypeError, KeyError) as error:
            record["source_error"] = str(error)
            continue

        if kind == "loop_link":
            record["record_reason"] = "loop_link_requires_single_peptide_support"
            continue
        if any(
            endpoint["status"] != "valid"
            for endpoint in record["endpoints"].values()
        ):
            record["record_reason"] = "source_endpoint_unresolved"
            continue
        if source_decoy is None:
            record["record_reason"] = "invalid_source_decoy_flag"
            continue

        writer.writerow({**fields, marker: record["source_record_id"]})
        prepared[record["source_record_id"]] = record

    parser_info = {
        "name": "pyXLMS",
        "version": version("pyXLMS"),
        "entrypoint": "read_scout",
        "format_variant": "scout_filtered_csm_columns",
        "parsed_records": 0,
        "modification_masses_assessed": False,
    }
    if prepared:
        try:
            parsed = parser.read_scout(
                io.BytesIO(stream.getvalue().encode("utf-8")),
                crosslinker=crosslinker,
                parse_modifications=False,
                sep=",",
                dtype=str,
                keep_default_na=False,
            )
            by_source = {}
            for csm in parsed["crosslink-spectrum-matches"] or []:
                source = (csm["additional_information"] or {}).get("source", {})
                source_id = source.get(marker)
                if source_id not in prepared or source_id in by_source:
                    raise ValueError("Unknown or repeated Scout parser source ID")
                by_source[source_id] = csm
            if set(by_source) != set(prepared):
                raise ValueError("Scout parser omitted source records")

            orientations = {}
            for source_id, record in prepared.items():
                csm = by_source[source_id]
                expected = [
                    (
                        endpoint["peptide"],
                        endpoint["link_position"],
                        (endpoint["source_protein_id"],),
                        (endpoint["protein_position"],),
                        (endpoint["peptide_start"],),
                    )
                    for endpoint in (
                        record["endpoints"]["a"], record["endpoints"]["b"]
                    )
                ]
                actual = [
                    (
                        csm[f"{side}_peptide"],
                        csm[f"{side}_peptide_crosslink_position"],
                        tuple(csm[f"{side}_proteins"] or []),
                        tuple(csm[f"{side}_proteins_crosslink_positions"] or []),
                        tuple(csm[f"{side}_proteins_peptide_positions"] or []),
                    )
                    for side in ("alpha", "beta")
                ]
                if actual != expected and actual != expected[::-1]:
                    raise ValueError("Scout parser changed source endpoints")
                if (
                    csm["scan_nr"] != record["spectrum"]["scan_number"]
                    or csm["spectrum_file"]
                    != record["spectrum"]["original_file"].strip()
                ):
                    raise ValueError("Scout parser changed spectrum identity")
                orientations[source_id] = (
                    "same" if actual == expected else "swapped"
                )

            for source_id, record in prepared.items():
                record.update(
                    record_status="valid",
                    record_reason="",
                    parser_orientation=orientations[source_id],
                )
            parser_info["parsed_records"] = len(by_source)
        except (
            ValueError, TypeError, KeyError, RuntimeError,
            IndexError, AttributeError, OverflowError,
        ) as error:
            for record in prepared.values():
                record.update(
                    record_status="not_evaluable",
                    record_reason="pyxlms_parse_failed",
                    parser_error=f"{type(error).__name__}: {error}",
                )

    return {
        "format": "scout_csm",
        "format_version": 1,
        "source_columns": columns,
        "records": records,
        "proteins": proteins,
        "parser": parser_info,
        "crosslinker": crosslinker,
        "crosslinker_label_source": "caller_declared",
        "provenance": {
            "dataset_accession": dataset_accession,
            "psm_path": str(psm_path.resolve()),
            "psm_sha256": digest,
            "fasta_path": str(fasta_path.resolve()),
            "fasta_sha256": next(iter(proteins.values()))["fasta_sha256"],
        },
    }
def _scout_export_evidence(bundle, manifest_path):
    """Verify the reviewed supplementary table and portable reference evidence."""
    from openpyxl import load_workbook

    manifest_path = Path(manifest_path).resolve()
    root = manifest_path.parent
    raw_manifest = manifest_path.read_bytes()
    policy = json.loads(raw_manifest)

    if (
        policy.get("dataset_accession")
        != bundle["provenance"]["dataset_accession"]
        or policy.get("source_csv_sha256")
        != bundle["provenance"]["psm_sha256"]
        or policy.get("search_fasta_sha256")
        != bundle["provenance"]["fasta_sha256"]
    ):
        raise ValueError("Scout author evidence source checksum mismatch")

    cutoff = policy.get("reported_fdr_cutoff")
    if (
        isinstance(cutoff, bool)
        or not isinstance(cutoff, (int, float))
        or not math.isfinite(cutoff)
        or not 0 <= cutoff <= 1
        or policy.get("per_record_fdr_available") is not False
        or policy.get("reported_fdr_level")
        not in {"unspecified", *FDR_COLUMNS}
    ):
        raise ValueError("Unsupported Scout author filtering declaration")

    bodies = {manifest_path.name: raw_manifest}
    checksums = policy.get("output_sha256")
    if not isinstance(checksums, dict):
        raise ValueError("Scout author evidence has no source checksums")

    for name, expected in checksums.items():
        relative = Path(name)
        source = (root / relative).resolve()
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or root not in source.parents
            or relative.as_posix() == manifest_path.name
        ):
            raise ValueError("Scout evidence path leaves its directory")
        raw = source.read_bytes()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError(f"Scout evidence checksum mismatch: {name}")
        bodies[relative.as_posix()] = raw

    required = {
        "author_results.xlsx",
        "membership_report.json",
        "candidate_report.json",
    }
    if not required <= bodies.keys():
        raise ValueError("Scout author or reference evidence is incomplete")

    with Path(bundle["provenance"]["psm_path"]).open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        source_rows = list(reader)

    if (
        columns != bundle["source_columns"]
        or source_rows
        != [record["source_fields"] for record in bundle["records"]]
    ):
        raise ValueError("Scout source records differ from their CSV")

    fields = [
        "Scan", "File", "Link-Type", "Precursor charge",
        "Alpha peptide", "Beta peptide",
        "Alpha peptide position", "Beta peptide position",
        "Alpha protein(s) position(s)", "Beta protein(s) position(s)",
        "Alpha protein mapping(s)", "Beta protein mapping(s)",
        "Alpha modification(s)", "Beta modification(s)",
        "Alpha gene mapping(s)", "Beta gene mapping(s)", "IsDecoy",
    ]

    def membership_key(row):
        values = []
        for field in fields:
            value = row[field]
            if value is None:
                text = ""
            elif isinstance(value, bool):
                text = str(value).lower()
            elif isinstance(value, float) and value.is_integer():
                text = str(int(value))
            else:
                text = str(value).strip()
            if field == "IsDecoy":
                text = text.lower()
            if (
                row["Link-Type"] == "Looplink"
                and field.startswith("Beta")
                and text in ("", "NA")
            ):
                text = ""
            values.append(text)
        return tuple(values)

    workbook = load_workbook(
        io.BytesIO(bodies["author_results.xlsx"]),
        read_only=True,
        data_only=False,
    )
    try:
        rows = workbook[policy["worksheet"]].iter_rows(values_only=True)
        headers = list(next(rows))
        if (
            len(headers) != len(set(headers))
            or not set(fields) <= set(headers)
        ):
            raise ValueError("Unsupported Scout supplementary columns")
        author_rows = [
            dict(zip(headers, values))
            for values in rows
            if any(value is not None for value in values)
        ]
    finally:
        workbook.close()

    source_keys = [membership_key(row) for row in source_rows]
    author_keys = [membership_key(row) for row in author_rows]
    if Counter(source_keys) != Counter(author_keys):
        raise ValueError("Scout records differ from the author results table")

    membership = json.loads(bodies["membership_report.json"])
    workbook_hash = hashlib.sha256(
        bodies["author_results.xlsx"]
    ).hexdigest()
    if (
        membership.get("source_csv_sha256")
        != bundle["provenance"]["psm_sha256"]
        or membership.get("supplement_xlsx_sha256") != workbook_hash
        or membership.get("worksheet") != policy["worksheet"]
        or membership.get("comparison_fields") != fields
        or membership.get("source_records") != len(source_rows)
        or membership.get("supplement_records") != len(author_rows)
        or membership.get("matched_records") != len(source_rows)
        or membership.get("missing_records") != 0
        or membership.get("extra_records") != 0
        or membership.get("same_record_order")
        != (source_keys == author_keys)
    ):
        raise ValueError("Scout supplementary membership review differs")

    if policy.get("expected_source_records") != len(source_rows):
        raise ValueError("Scout author evidence record count differs")

    reports = json.loads(bodies["candidate_report.json"])
    if not isinstance(reports, dict):
        raise ValueError("Invalid Scout reference candidate report")

    for report in reports.values():
        for receipt in report.get("sources", []):
            relative = Path(receipt["path"])
            name = relative.as_posix()
            if name not in bodies:
                raise ValueError("Scout reference evidence is missing")
            if hashlib.sha256(bodies[name]).hexdigest() != receipt["sha256"]:
                raise ValueError("Scout reference evidence checksum mismatch")
            receipt["absolute_path"] = str((root / relative).resolve())

    return policy, reports, bodies


def write_scout_claudio_input(
    bundle,
    reference_projection,
    output_dir,
    *,
    author_filter_path,
    allowed_crosslinkers,
    allowed_residues,
):
    """Export verified Scout pairs while retaining author and reference evidence."""
    if bundle.get("format") != "scout_csm":
        raise ValueError("Scout input export requires a Scout CSM import")
    if Path(output_dir).exists():
        raise ValueError("Output directory already exists")

    selections = {}
    for name, values in (
        ("crosslinkers", allowed_crosslinkers),
        ("residues", allowed_residues),
    ):
        if (
            not isinstance(values, (set, frozenset, list, tuple))
            or not values
            or any(
                not isinstance(value, str)
                or not value
                or value != value.strip()
                for value in values
            )
        ):
            raise ValueError(f"Explicit nonempty {name} are required")
        selections[name] = set(values)

    if any(
        re.fullmatch(r"[A-Z]", residue) is None
        for residue in selections["residues"]
    ):
        raise ValueError("Residues must be uppercase one-letter codes")

    for kind in ("psm", "fasta"):
        raw = Path(bundle["provenance"][kind + "_path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != \
                bundle["provenance"][kind + "_sha256"]:
            raise ValueError(f"Scout {kind} source checksum mismatch")

    proteins = read_search_fasta(bundle["provenance"]["fasta_path"])
    if (
        proteins.keys() != bundle["proteins"].keys()
        or any(
            protein["sequence"] != bundle["proteins"][identifier]["sequence"]
            for identifier, protein in proteins.items()
        )
    ):
        raise ValueError("Scout search FASTA sequences changed")

    policy, reports, bodies = _scout_export_evidence(
        bundle, author_filter_path
    )

    if (
        reference_projection.get("kind") != "xifdr_reference_projection"
        or reference_projection.get("projection_rule_version")
        != "intact-peptide-projection-v1"
    ):
        raise ValueError("Unsupported Scout reference projection")

    def source_context(provenance):
        return {
            key: value for key, value in provenance.items()
            if not key.endswith("_path")
        }

    if source_context(reference_projection["provenance"]) != \
            source_context(bundle["provenance"]):
        raise ValueError("Scout projection source checksums differ")

    mappings = reference_projection.get(
        "identity_evidence", {}
    ).get("exact_mappings", {})
    checked = project_xifdr_candidates(
        bundle,
        mappings,
        reports,
        max_optimal_alignments=reference_projection.get(
            "max_optimal_alignments"
        ),
    )
    if checked["records"] != reference_projection.get("records"):
        raise ValueError("Scout projection coordinates or source order changed")

    def profile_index(projection):
        profiles = {}
        for values in projection["protein_candidates"].values():
            for profile in values:
                identifier = profile["profile_id"]
                if identifier in profiles:
                    raise ValueError("Repeated Scout reference projection profile")
                profiles[identifier] = profile
        return profiles

    profiles = profile_index(checked)
    original_profiles = profile_index(reference_projection)
    if (
        profiles.keys() != original_profiles.keys()
        or any(
            {key: value for key, value in profile.items() if key != "sources"}
            != {
                key: value
                for key, value in original_profiles[identifier].items()
                if key != "sources"
            }
            for identifier, profile in profiles.items()
        )
    ):
        raise ValueError("Scout reference projection profiles changed")

    for profile in profiles.values():
        identifier = profile["source_fasta_id"]
        report = reports.get(identifier)
        if report is None:
            report = {
                "taxon_id": profile["taxon_id"],
                "sources": mappings.get(identifier, {}).get("sources", []),
            }
        _reference_evidence(
            profile,
            report,
            proteins[identifier],
            require_declared_accession=False,
        )

    working = deepcopy(bundle)
    declaration_path = Path(author_filter_path).resolve()
    declaration_hash = hashlib.sha256(
        declaration_path.read_bytes()
    ).hexdigest()
    working["provenance"].update(
        author_filter_path=str(declaration_path),
        author_filter_sha256=declaration_hash,
    )
    row_policy = {
        "kind": "author_filtered_psm_table",
        "level": policy["reported_fdr_level"],
        "reported_cutoff": policy["reported_fdr_cutoff"],
        "per_record_fdr_applied": False,
    }

    rows, audit = [], []
    for record, projected in zip(working["records"], checked["records"]):
        entry = deepcopy(record)
        reasons = []
        entry["source_filter_status"] = "author_reported"
        entry["input_status"] = "not_evaluable"
        audit.append(entry)

        if record.get("record_kind") != "paired_crosslink_candidate":
            reasons.append(
                record.get("record_reason")
                or "unsupported_scout_record_type"
            )
        elif record["record_status"] != "valid":
            reasons.append(record["record_reason"] or "scout_source_unresolved")

        if (
            record["decoy_status"] != "target_target"
            or _boolean(record["source_fields"].get("IsDecoy")) is not False
        ):
            reasons.append("source_target_status_unresolved")

        if reasons:
            entry["input_reasons"] = list(dict.fromkeys(reasons))
            continue

        for site in ("a", "b"):
            endpoint = _scout_endpoint(
                record["source_fields"], site, proteins
            )
            if (
                endpoint["status"] != "valid"
                or endpoint != record["endpoints"][site]
            ):
                reasons.append("source_endpoint_changed_" + site)

        if len(projected["pairs"]) != 1:
            reasons.append("reference_projection_unresolved")

        if reasons:
            entry["input_reasons"] = list(dict.fromkeys(reasons))
            continue

        pair = projected["pairs"][0]
        selected = {
            site: profiles[pair["profile_id_" + site]]
            for site in ("a", "b")
        }
        identities = {
            site: {
                "accession": profile["accession"],
                "sequence": profile["reference_sequence"],
            }
            for site, profile in selected.items()
        }
        row = _input_row(working, record, identities, row_policy)
        row.update(
            source_psm_id=record["source_fields"]["Scan"],
            source_scan=record["source_fields"]["Scan"],
            source_run=record["source_fields"]["File"],
            source_crosslinker=bundle["crosslinker"],
            source_link_type=record["source_fields"]["Link-Type"],
            selection_basis="author_reported_csm_table",
            projection_candidate_id=pair["projection_candidate_id"],
            identity_scope=pair["identity_scope"],
        )
        reference_endpoints = {}
        for site, label in (("a", "Alpha"), ("b", "Beta")):
            candidate = next(
                endpoint
                for endpoint in projected["endpoint_candidates"][site]
                if endpoint["profile_id"] == pair["profile_id_" + site]
            )
            reference_endpoints[site] = deepcopy(candidate)
            row.update({
                "pos_" + site: candidate["reference_protein_position"],
                "source_peptide_start_" + site:
                    candidate["source_peptide_start"],
                "reference_peptide_start_" + site:
                    candidate["reference_peptide_start"],
                "reference_mapping_method_" + site:
                    candidate["mapping_method"],
                "reference_identity_scope_" + site:
                    candidate["identity_scope"],
                "source_modifications_" + site:
                    record["source_fields"].get(label + " modification(s)", ""),
            })
            residue = row["seq_" + site][row["pos_" + site] - 1]
            if residue not in selections["residues"]:
                reasons.append("residue_not_selected_" + site)

        if bundle["crosslinker"] not in selections["crosslinkers"]:
            reasons.append("crosslinker_not_selected")

        entry.update(
            reference_endpoints=reference_endpoints,
            projection_candidate_id=pair["projection_candidate_id"],
            identity_scope=pair["identity_scope"],
            input_reasons=list(dict.fromkeys(reasons)),
        )
        if not reasons:
            entry["input_status"] = "ready"
            rows.append(row)

    prepared = {
        "rows": rows,
        "audit": audit,
        "provenance": deepcopy(working["provenance"]),
        "fdr_policy": {
            **row_policy,
            "kind": "author_reported_csm_table",
        },
        "reference_projection": checked,
        "reference_projection_policy": {
            "scope": "evidence_supported_reference_candidates",
            "intact_peptide_required": True,
            "unique_coordinates_required": True,
            "identity_acceptance_performed": False,
        },
        "author_filter_policy": deepcopy(policy),
        "author_filter_evidence_snapshot":
            "author_filter_sources/" + declaration_path.name,
        "extra_source_bodies": {
            "author_filter_sources/" + name: raw
            for name, raw in bodies.items()
        },
    }
    return write_xifdr_outputs(
        working,
        mappings,
        prepared,
        output_dir,
        pipeline_policy={
            "allowed_crosslinkers": sorted(selections["crosslinkers"]),
            "allowed_residues": sorted(selections["residues"]),
            "source_protein_ids": None,
        },
    )
def _xlinkx_endpoint(fields, site, proteins):
    """Validate one source-reported leading mapping against the search FASTA."""
    suffix = "A" if site == "a" else "B"
    result = {
        "source_protein_id": fields.get("ProteinAccession" + suffix),
        "source_peptide": fields.get("Sequence" + suffix),
        "source_modifications": fields.get("Modifications" + suffix),
        "source_all_protein_annotations": fields.get(
            "ProteinAccessions" + suffix
        ),
        "source_mapping_scope": "source_reported_leading_mapping",
        "source_fasta_id": None,
        "source_fasta_candidates": [],
        "peptide": None,
        "link_position": None,
        "protein_position": None,
        "peptide_start": None,
        "peptide_end": None,
        "kind": "residue",
        "status": "not_evaluable",
        "reason": "",
    }
    try:
        raw_peptide = fields.get("Sequence" + suffix)
        if not isinstance(raw_peptide, str):
            raise ValueError("Native peptide format requires review")
        peptide = raw_peptide.strip()
        if re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", peptide) is None:
            raise ValueError("Native peptide notation requires review")

        result["peptide"] = peptide
        link = _positive_integer(
            str(fields.get("CrosslinkerPosition" + suffix))
        )
        if link > len(peptide):
            raise ValueError("Native peptide link position is outside the peptide")
        result["link_position"] = link

        identifiers = [
            value.strip()
            for value in str(
                fields.get("ProteinAccession" + suffix, "")
            ).split(";")
        ]
        positions = [
            _positive_integer(value)
            for value in str(
                fields.get("LeadingProteinPosition" + suffix, "")
            ).split(";")
        ]
        result["source_protein_candidates"] = identifiers
        result["source_position_candidates"] = positions
        if not all(identifiers) or len(identifiers) != len(positions):
            raise ValueError("Native protein and position alternatives differ")
        if len(identifiers) != 1:
            result.update(
                status="ambiguous",
                reason="multiple_reported_mappings_require_review",
            )
            return result

        identifier, position = identifiers[0], positions[0]
        start = position - link + 1
        end = start + len(peptide) - 1
        result.update(
            source_protein_id=identifier,
            protein_position=position,
            peptide_start=start,
            peptide_end=end,
        )
        if start < 1:
            raise ValueError("Native peptide span starts before the protein")

        aliases = [
            (key, protein)
            for key, protein in proteins.items()
            if identifier in (key, protein.get("uniprot_accession"))
        ]
        if not aliases:
            raise ValueError("Source protein is not identified in the search FASTA")

        matches = [
            key for key, protein in aliases
            if protein["sequence"][start - 1:end] == peptide
        ]
        result["source_fasta_candidates"] = matches
        if len(matches) == 1:
            result.update(
                source_fasta_id=matches[0], status="valid", reason=""
            )
        elif matches:
            result.update(
                status="ambiguous",
                reason="multiple_fasta_sequence_and_position_matches",
            )
        else:
            result["reason"] = "source_peptide_position_mismatch"
    except (ValueError, TypeError, KeyError) as error:
        result["reason"] = str(error)
    return result


def read_xlinkx_csms(psm_path, fasta_path, *, dataset_accession):
    """Read native XlinkX CSM tables without mixing result levels or inferring FDR."""
    import sqlite3
    from importlib.metadata import version

    try:
        from pyXLMS import parser
    except ImportError as error:
        raise RuntimeError(
            "Reading XlinkX CSMs requires pyXLMS==2.0.7"
        ) from error

    if (
        not isinstance(dataset_accession, str)
        or re.fullmatch(r"PXD\d{6,}", dataset_accession) is None
    ):
        raise ValueError("A ProteomeXchange PXD accession is required")

    psm_path = Path(psm_path).resolve()
    fasta_path = Path(fasta_path).resolve()
    digest = hashlib.sha256()
    with psm_path.open("rb") as handle:
        signature = handle.read(16)
        if signature != b"SQLite format 3\x00":
            raise ValueError("Expected a native SQLite XlinkX pdResult file")
        digest.update(signature)
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    source_hash = digest.hexdigest()
    proteins = read_search_fasta(fasta_path)

    column_mapping = {
        "SequenceA": "Sequence A",
        "SequenceB": "Sequence B",
        "ModificationsA": "Modifications A",
        "ModificationsB": "Modifications B",
        "CrosslinkerPositionA": "Crosslinker Position A",
        "CrosslinkerPositionB": "Crosslinker Position B",
        "ProteinAccessionA": "Protein Accession A",
        "ProteinAccessionB": "Protein Accession B",
        "LeadingProteinPositionA": "Leading Protein Position A",
        "LeadingProteinPositionB": "Leading Protein Position B",
        "IsDecoy": "Is Decoy",
        "XlinkXScore": "XlinkX Score",
        "SpectrumFileName": "Spectrum File",
        "FirstScan": "First Scan",
        "ChargeState": "Charge",
        "RetentionTime": "RT [min]",
    }
    required = set(column_mapping) | {"WorkflowId", "ID", "Crosslinker"}
    tables = ("CSMs", "DecoyCSMs", "Crosslinks")
    native_rows, table_columns = {}, {}

    connection = sqlite3.connect(
        psm_path.as_uri() + "?mode=ro&immutable=1", uri=True
    )
    connection.row_factory = sqlite3.Row
    try:
        available = {
            row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        if not set(tables) <= available:
            raise ValueError("Unsupported XlinkX native table layout")
        for table in tables:
            cursor = connection.execute(
                f'SELECT * FROM "{table}" ORDER BY WorkflowId, ID'
            )
            columns = [item[0] for item in cursor.description]
            if table != "Crosslinks" and not required <= set(columns):
                missing = sorted(required - set(columns))
                raise ValueError(
                    f"Unsupported XlinkX {table} columns: {missing}"
                )
            table_columns[table] = columns
            native_rows[table] = [
                dict(row) for row in cursor.fetchall()
            ]
    finally:
        connection.close()

    def source_value(value):
        if isinstance(value, bytes):
            return {"sqlite_blob_hex": value.hex()}
        if isinstance(value, float) and not math.isfinite(value):
            return {"sqlite_non_finite_float": repr(value)}
        return value

    def source_fields(row):
        return {
            key: source_value(value) for key, value in row.items()
        }

    records, prepared, seen = [], {}, set()
    marker = "__claudio_xlinkx_source_id"
    for table in ("CSMs", "DecoyCSMs"):
        for native in native_rows[table]:
            fields = source_fields(native)
            source_id = (
                f"{dataset_accession}:{source_hash}:{table}:"
                f"{native['WorkflowId']}:{native['ID']}"
            )
            if source_id in seen:
                raise ValueError("Native XlinkX CSM row keys are not unique")
            seen.add(source_id)

            flag = _boolean(native.get("IsDecoy"))
            conflict = table == "DecoyCSMs" and flag is False
            if conflict:
                decoy_reason = "native_decoy_table_flag_conflict"
            elif flag is True:
                decoy_reason = "source_decoy_endpoint_assignment_unavailable"
            elif flag is None:
                decoy_reason = "invalid_source_decoy_flag"
            else:
                decoy_reason = ""

            record = {
                "source_record_id": source_id,
                "source_record_number": len(records) + 1,
                "source_table": table,
                "source_native_key": {
                    "WorkflowId": fields["WorkflowId"],
                    "ID": fields["ID"],
                },
                "source_fields": fields,
                "record_kind": "paired_crosslink_candidate",
                "crosslinker": fields.get("Crosslinker"),
                "source_decoy_flag": flag,
                "decoy_status": (
                    "target_target"
                    if flag is False and not conflict
                    else "unknown"
                ),
                "decoy_reason": decoy_reason,
                "fdr": {
                    level: _fdr_value(None) for level in FDR_COLUMNS
                },
                "modification_status": "source_strings_retained",
                "endpoints": {},
                "spectrum": None,
                "record_status": "not_evaluable",
                "record_reason": "xlinkx_source_unresolved",
            }
            records.append(record)

            try:
                filename = native.get("SpectrumFileName")
                if not isinstance(filename, str) or not filename.strip():
                    raise ValueError("Missing native spectrum filename")
                record["spectrum"] = {
                    "scan_number": _positive_integer(
                        str(native.get("FirstScan"))
                    ),
                    "original_file": filename,
                    "charge": _positive_integer(
                        str(native.get("ChargeState"))
                    ),
                }
            except (ValueError, TypeError) as error:
                record.update(
                    record_reason="invalid_spectrum_identity",
                    source_error=str(error),
                )
                continue

            record["endpoints"] = {
                site: _xlinkx_endpoint(fields, site, proteins)
                for site in ("a", "b")
            }
            if any(
                endpoint["status"] != "valid"
                for endpoint in record["endpoints"].values()
            ):
                record["record_reason"] = "source_endpoint_unresolved"
                continue
            if flag is None or conflict:
                record["record_reason"] = decoy_reason
                continue
            prepared[source_id] = (record, native)

    parser_info = {
        "name": "pyXLMS",
        "version": version("pyXLMS"),
        "entrypoint": "read_xlinkx",
        "format_variant": "xlinkx_pdresult_csm_tables",
        "parsed_records": 0,
        "modification_masses_assessed": False,
        "structured_modifications_parsed": False,
    }

    if prepared:
        stream = io.StringIO()
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(list(column_mapping.values()) + [marker])
        for source_id, (_, native) in prepared.items():
            writer.writerow(
                [native[column] for column in column_mapping]
                + [source_id]
            )

        try:
            result = parser.read_xlinkx(
                io.BytesIO(stream.getvalue().encode("utf-8")),
                format="csv",
                sep="\t",
                dtype=str,
                keep_default_na=False,
                parse_modifications=False,
                ignore_errors=False,
            )
            parsed = result["crosslink-spectrum-matches"] or []
            by_source = {}
            for csm in parsed:
                source = (csm["additional_information"] or {}).get(
                    "source", {}
                )
                source_id = source.get(marker)
                if source_id not in prepared or source_id in by_source:
                    raise ValueError(
                        "XlinkX parser source associations differ"
                    )
                by_source[source_id] = csm
            if by_source.keys() != prepared.keys():
                raise ValueError("XlinkX parser omitted source records")

            def expected_endpoint(endpoint):
                return (
                    endpoint["peptide"],
                    endpoint["link_position"],
                    (endpoint["source_protein_id"],),
                    (endpoint["protein_position"],),
                    (endpoint["peptide_start"],),
                )

            def parsed_endpoint(csm, side):
                return (
                    csm[side + "_peptide"],
                    csm[side + "_peptide_crosslink_position"],
                    tuple(csm[side + "_proteins"] or []),
                    tuple(csm[side + "_proteins_crosslink_positions"] or []),
                    tuple(csm[side + "_proteins_peptide_positions"] or []),
                )

            updates = []
            for source_id, (record, _) in prepared.items():
                csm = by_source[source_id]
                expected = (
                    expected_endpoint(record["endpoints"]["a"]),
                    expected_endpoint(record["endpoints"]["b"]),
                )
                actual = (
                    parsed_endpoint(csm, "alpha"),
                    parsed_endpoint(csm, "beta"),
                )
                if actual == expected:
                    orientation = "same"
                elif actual == expected[::-1]:
                    orientation = "swapped"
                else:
                    raise ValueError(
                        f"XlinkX parser endpoints differ: {source_id}"
                    )

                spectrum = record["spectrum"]
                if (
                    csm["scan_nr"] != spectrum["scan_number"]
                    or csm["charge"] != spectrum["charge"]
                    or csm["spectrum_file"]
                    != spectrum["original_file"].strip()
                ):
                    raise ValueError(
                        f"XlinkX parser spectrum identity differs: {source_id}"
                    )
                if record["source_decoy_flag"] is False and (
                    csm["alpha_decoy"] is not False
                    or csm["beta_decoy"] is not False
                ):
                    raise ValueError(
                        f"XlinkX parser target flags differ: {source_id}"
                    )
                updates.append((record, orientation))

            for record, orientation in updates:
                record.update(
                    record_status="valid",
                    record_reason="",
                    parser_orientation=orientation,
                )
            parser_info["parsed_records"] = len(updates)
        except (
            ValueError, TypeError, KeyError, RuntimeError,
            IndexError, AttributeError, OverflowError,
        ) as error:
            parser_info["error"] = f"{type(error).__name__}: {error}"
            for record, _ in prepared.values():
                record.update(
                    record_status="not_evaluable",
                    record_reason="pyxlms_parse_failed",
                    parser_error=parser_info["error"],
                )
                record.pop("parser_orientation", None)

    return {
        "format": "xlinkx_csm",
        "format_version": 1,
        "source_columns": table_columns["CSMs"],
        "source_table_columns": table_columns,
        "native_table_counts": {
            table: len(native_rows[table]) for table in tables
        },
        "native_crosslinks": [
            source_fields(row) for row in native_rows["Crosslinks"]
        ],
        "records": records,
        "proteins": proteins,
        "parser": parser_info,
        "provenance": {
            "dataset_accession": dataset_accession,
            "psm_path": str(psm_path),
            "psm_sha256": source_hash,
            "fasta_path": str(fasta_path),
            "fasta_sha256": hashlib.sha256(
                fasta_path.read_bytes()
            ).hexdigest(),
        },
        "fdr_filter_applied": False,
        "identity_mapping_applied": False,
        "input_exported": False,
        "assessment_performed": False,
    }
def _xlinkx_native_validation(bundle):
    """Read native confidence, workflow settings and execution evidence."""
    import sqlite3
    from xml.etree import ElementTree as ET

    def encoded(value):
        if isinstance(value, dict):
            return {key: encoded(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [encoded(item) for item in value]
        if isinstance(value, bytes):
            return {"sqlite_blob_hex": value.hex()}
        if isinstance(value, float) and not math.isfinite(value):
            return {"sqlite_non_finite_float": repr(value)}
        return value

    native = Path(bundle["provenance"]["psm_path"]).resolve()
    connection = sqlite3.connect(
        native.as_uri() + "?mode=ro&immutable=1", uri=True
    )
    connection.row_factory = sqlite3.Row
    try:
        available = {
            row[0] for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }

        def read_table(table, ordered=False):
            if table not in available:
                return []
            sql = f'SELECT * FROM "{table}"'
            if ordered:
                sql += " ORDER BY WorkflowId, ID"
            return [
                dict(row) for row in connection.execute(sql)
            ]

        tables = {
            table: read_table(
                table,
                ordered=table in {"CSMs", "DecoyCSMs", "Crosslinks"},
            )
            for table in (
                "CSMs", "DecoyCSMs", "Crosslinks", "CrosslinksCSMs",
                "DataTypes", "DataTypesColumns", "EnumDataTypeValues",
                "Workflows", "WorkflowMessages",
            )
        }
    finally:
        connection.close()

    expected = [
        (
            table,
            {"WorkflowId": row["WorkflowId"], "ID": row["ID"]},
            encoded(row),
        )
        for table in ("CSMs", "DecoyCSMs")
        for row in tables[table]
    ]
    actual = [
        (
            record["source_table"],
            record["source_native_key"],
            record["source_fields"],
        )
        for record in bundle["records"]
    ]
    if actual != expected:
        raise ValueError("XlinkX source records changed after import")
    if bundle["native_crosslinks"] != encoded(tables["Crosslinks"]):
        raise ValueError("XlinkX native crosslinks changed after import")

    labels, definitions = {}, []
    for table, field in (
        ("CSMs", "MatchConfidence"),
        ("DecoyCSMs", "MatchConfidence"),
        ("Crosslinks", "Confidence"),
    ):
        type_ids = {
            row["DataTypeID"] for row in tables["DataTypes"]
            if row["TableName"] == table
        }
        columns = [
            row for row in tables["DataTypesColumns"]
            if row["DataTypeID"] in type_ids
            and row["DBColumnName"] == field
            and row["SpecialValueType"] == "Enum"
        ]
        labels[table] = {}
        if len(columns) == 1:
            column = columns[0]
            values = [
                row for row in tables["EnumDataTypeValues"]
                if row["EnumID"] == column["SpecialValueTypeID"]
            ]
            labels[table] = {
                str(row["Value"]): row["DisplayName"]
                for row in values
            }
            definitions.append({
                "table": table,
                "field": field,
                "column": column,
                "enum_values": values,
            })

    def local_tag(element):
        return str(element.tag).rsplit("}", 1)[-1]

    profiles, warnings = {}, []
    workflow_ids = {
        record["source_native_key"]["WorkflowId"]
        for record in bundle["records"]
    }
    for workflow_id in sorted(workflow_ids):
        profile = {
            "workflow_id": workflow_id,
            "node_number": None,
            "node_state": None,
            "configured_fdr_cutoff": None,
            "configured_strategy": None,
            "observed_strategy": "not_established",
            "fdr_verification_status": "unverified",
        }
        workflows = [
            row for row in tables["Workflows"]
            if row["WorkflowID"] == workflow_id
        ]
        nodes = []
        if len(workflows) == 1:
            try:
                root = ET.fromstring(workflows[0]["WorkflowXML"])
                nodes = [
                    node for node in root.iter()
                    if local_tag(node) == "WorkflowNode"
                    and node.get("ProcessingNodeName") == "HlxlValidator"
                ]
            except (ET.ParseError, TypeError, KeyError) as error:
                profile["metadata_error"] = str(error)

        if len(nodes) == 1:
            node = nodes[0]
            number = node.get("ProcessingNodeNumber")
            profile.update(
                node_number=number,
                node_state=node.get("NodeState"),
            )
            parameters = [
                parameter
                for container in node
                if local_tag(container) == "ProcessingNodeParameters"
                for parameter in container
                if local_tag(parameter) == "ProcessingNodeParameter"
                and parameter.get("IsValueSet", "").lower() == "true"
            ]

            def parameter_values(name):
                return [
                    (parameter.text or "").strip()
                    for parameter in parameters
                    if parameter.get("Name") == name
                ]

            cutoffs = parameter_values("FdrThreshold")
            if len(cutoffs) == 1:
                try:
                    cutoff = float(cutoffs[0])
                    if math.isfinite(cutoff) and 0 <= cutoff <= 1:
                        profile["configured_fdr_cutoff"] = cutoff
                except ValueError:
                    pass

            strategies = parameter_values("FdrStrategy")
            if len(strategies) == 1:
                profile["configured_strategy"] = strategies[0]

            messages = []
            for message in tables["WorkflowMessages"]:
                match = re.match(
                    r"^\(([^)]+)\):",
                    str(message.get("ProcessingNodeName", "")),
                )
                if (
                    message["WorkflowID"] == workflow_id
                    and match is not None
                    and match.group(1) == number
                    and isinstance(message.get("Message"), str)
                ):
                    messages.append(message)

            if any(
                "reverting to simple validation"
                in message["Message"].lower()
                for message in messages
            ):
                profile["observed_strategy"] = "simple_validation"

            warnings.extend(
                message for message in messages
                if re.search(
                    r"warning|error|exception|no decoy|revert|unable",
                    message["Message"],
                    re.IGNORECASE,
                )
            )
        else:
            profile["metadata_error"] = (
                "A unique XlinkX validator node was not identified"
            )

        profiles[str(workflow_id)] = profile

    crosslinks = {}
    for row in tables["Crosslinks"]:
        key = (row["WorkflowId"], row["ID"])
        if key in crosslinks:
            raise ValueError("Repeated native crosslink keys")
        crosslinks[key] = row

    parents = {}
    for relation in tables["CrosslinksCSMs"]:
        child = (relation["CSMsWorkflowId"], relation["CSMsID"])
        parent = (
            relation["CrosslinksWorkflowId"],
            relation["CrosslinksID"],
        )
        parents.setdefault(child, set()).add(parent)

    record_evidence = {}
    for record in bundle["records"]:
        key = record["source_native_key"]
        child = (key["WorkflowId"], key["ID"])
        parent_keys = (
            parents.get(child, set())
            if record["source_table"] == "CSMs"
            else set()
        )
        info = {
            "csm_confidence": labels.get(
                record["source_table"], {}
            ).get(
                str(record["source_fields"].get("MatchConfidence")),
                "unavailable",
            ),
            "parent_keys": [
                {"WorkflowId": parent[0], "ID": parent[1]}
                for parent in sorted(parent_keys)
            ],
            "parent_status": "missing",
            "parent_confidence": "unavailable",
            "validator": profiles.get(str(key["WorkflowId"]), {}),
        }
        if len(parent_keys) > 1:
            info["parent_status"] = "multiple"
        elif len(parent_keys) == 1:
            parent_key = next(iter(parent_keys))
            parent = crosslinks.get(parent_key)
            if parent is None:
                info["parent_status"] = "invalid"
            else:
                info.update(
                    parent_status="unique",
                    parent_confidence=labels["Crosslinks"].get(
                        str(parent["Confidence"]), "unavailable"
                    ),
                )
        record_evidence[record["source_record_id"]] = info

    def common_value(name, default=None):
        values = [profile.get(name) for profile in profiles.values()]
        if values and all(value == values[0] for value in values):
            return values[0]
        return default

    policy = {
        "kind": "native_reported_candidates",
        "configured_fdr_cutoff": common_value("configured_fdr_cutoff"),
        "configured_strategy": common_value("configured_strategy"),
        "observed_strategy": common_value(
            "observed_strategy", "not_established"
        ),
        "fdr_verification_status": "unverified",
        "per_record_fdr_applied": False,
        "validator_profiles": list(profiles.values()),
        "runtime_warnings": warnings,
    }
    evidence = {
        "kind": "xlinkx_native_validation_evidence",
        "native_sha256": bundle["provenance"]["psm_sha256"],
        "confidence_definitions": definitions,
        "crosslinks": tables["Crosslinks"],
        "crosslink_csm_relationships": tables["CrosslinksCSMs"],
        "workflows": tables["Workflows"],
        "workflow_messages": tables["WorkflowMessages"],
        "policy": policy,
    }
    body = (
        json.dumps(
            encoded(evidence),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        ) + "\n"
    ).encode("utf-8")
    return policy, record_evidence, body


def write_xlinkx_claudio_input(
    bundle,
    mappings,
    output_dir,
    *,
    allowed_crosslinkers,
    allowed_residues,
):
    """Export source-reported candidates with exact identities and explicit FDR status."""
    if bundle.get("format") != "xlinkx_csm":
        raise ValueError("XlinkX input export requires an XlinkX CSM import")
    if Path(output_dir).exists():
        raise ValueError("Output directory already exists")

    selections = {}
    for name, values in (
        ("crosslinkers", allowed_crosslinkers),
        ("residues", allowed_residues),
    ):
        if (
            not isinstance(values, (set, frozenset, list, tuple))
            or not values
            or any(
                not isinstance(value, str)
                or not value
                or value != value.strip()
                for value in values
            )
        ):
            raise ValueError(f"Explicit nonempty {name} are required")
        selections[name] = set(values)
    if any(
        re.fullmatch(r"[A-Z]", residue) is None
        for residue in selections["residues"]
    ):
        raise ValueError("Residues must be uppercase one-letter codes")

    for kind in ("psm", "fasta"):
        source = Path(bundle["provenance"][kind + "_path"])
        if hashlib.sha256(source.read_bytes()).hexdigest() != \
                bundle["provenance"][kind + "_sha256"]:
            raise ValueError(f"XlinkX {kind} source checksum mismatch")

    proteins = read_search_fasta(bundle["provenance"]["fasta_path"])
    identities = {}
    for identifier in target_proteins(bundle):
        source = proteins.get(identifier)
        identity = mappings.get(identifier)
        if (
            source is None
            or not isinstance(identity, dict)
            or identity.get("status") != "mapped"
            or identity.get("sequence") != source["sequence"]
            or not isinstance(identity.get("accession"), str)
            or ACCESSION.fullmatch(identity["accession"]) is None
        ):
            continue
        _reference_evidence(
            {
                "accession": identity["accession"],
                "taxon_id": identity.get("taxon_id"),
                "reference_sequence": identity["sequence"],
            },
            {
                "taxon_id": identity.get("taxon_id"),
                "sources": identity.get("sources", []),
            },
            source,
            require_declared_accession=False,
        )
        identities[identifier] = identity

    policy, native_evidence, evidence_body = _xlinkx_native_validation(
        bundle
    )
    working = deepcopy(bundle)
    working["proteins"] = proteins
    working["native_filter_policy"] = policy
    rows, audit = [], []

    for original in working["records"]:
        entry = deepcopy(original)
        info = native_evidence[original["source_record_id"]]
        validator = info["validator"]
        entry.update(
            input_status="not_evaluable",
            input_reasons=[],
            native_validation_evidence=deepcopy(info),
            fdr_verification_status="unverified",
        )
        audit.append(entry)
        reasons = entry["input_reasons"]

        if original.get("record_kind") != "paired_crosslink_candidate":
            reasons.append("unsupported_xlinkx_record_type")
        if original["record_status"] != "valid":
            reasons.append(original["record_reason"] or "source_unresolved")
        if (
            original["decoy_status"] != "target_target"
            or original["source_decoy_flag"] is not False
        ):
            reasons.append("source_target_status_unresolved")

        parent_reason = {
            "missing": "native_parent_crosslink_missing",
            "multiple": "multiple_native_parent_crosslinks",
            "invalid": "native_parent_crosslink_invalid",
        }.get(info["parent_status"])
        if parent_reason:
            reasons.append(parent_reason)
        if reasons:
            continue

        selected = {}
        for site in ("a", "b"):
            endpoint = _xlinkx_endpoint(
                original["source_fields"], site, proteins
            )
            if (
                endpoint["status"] != "valid"
                or endpoint != original["endpoints"][site]
            ):
                reasons.append("source_endpoint_changed_" + site)
                continue

            identity = identities.get(endpoint["source_fasta_id"])
            if identity is None:
                reasons.append("protein_identity_unresolved_" + site)
                continue
            selected[site] = identity
            residue = identity["sequence"][
                endpoint["protein_position"] - 1
            ]
            if residue not in selections["residues"]:
                reasons.append("residue_not_selected_" + site)

        if original["crosslinker"] not in selections["crosslinkers"]:
            reasons.append("crosslinker_not_selected")
        if reasons:
            continue

        row = {
            "source_record_id": original["source_record_id"],
            "dataset_accession": working["provenance"]["dataset_accession"],
            "source_psm_sha256": working["provenance"]["psm_sha256"],
            "source_fasta_sha256": working["provenance"]["fasta_sha256"],
            "source_record_number": original["source_record_number"],
            "source_table": original["source_table"],
            "source_native_key_json": json.dumps(
                original["source_native_key"]
            ),
            "source_run": original["spectrum"]["original_file"],
            "source_scan": original["spectrum"]["scan_number"],
            "source_crosslinker": original["crosslinker"],
            "source_fields_json": json.dumps(
                original["source_fields"], ensure_ascii=False
            ),
            "identity_scope": "exact_sequence",
            "selection_basis": "native_reported_candidate",
            "fdr_verification_status": "unverified",
            "source_validation_strategy": validator.get(
                "observed_strategy", "not_established"
            ),
            "source_configured_fdr_cutoff": validator.get(
                "configured_fdr_cutoff"
            ),
            "source_native_csm_confidence": info["csm_confidence"],
            "source_native_parent_confidence": info["parent_confidence"],
            "selected_fdr_level": None,
            "selected_fdr_value": None,
            "selected_fdr_cutoff": None,
        }
        for site in ("a", "b"):
            endpoint = original["endpoints"][site]
            identity = selected[site]
            row.update({
                "unip_id_" + site: identity["accession"],
                "seq_" + site: identity["sequence"],
                "pep_" + site: endpoint["peptide"],
                "pos_" + site: endpoint["protein_position"],
                "res_pos_" + site: endpoint["link_position"],
                "source_protein_id_" + site: endpoint["source_protein_id"],
                "source_protein_position_" + site:
                    endpoint["protein_position"],
                "source_peptide_start_" + site: endpoint["peptide_start"],
                "source_modifications_" + site:
                    endpoint["source_modifications"],
            })
        entry.update(input_status="ready", identity_scope="exact_sequence")
        rows.append(row)

    prepared = {
        "rows": rows,
        "audit": audit,
        "provenance": deepcopy(working["provenance"]),
        "fdr_policy": {
            "kind": "native_reported_candidates",
            "fdr_verification_status": "unverified",
            "per_record_fdr_applied": False,
        },
        "extra_source_bodies": {
            "native_validation_sources/native_validation.json":
                evidence_body,
        },
    }
    return write_xifdr_outputs(
        working,
        mappings,
        prepared,
        output_dir,
        pipeline_policy={
            "allowed_crosslinkers": sorted(selections["crosslinkers"]),
            "allowed_residues": sorted(selections["residues"]),
            "source_protein_ids": None,
        },
    )
def _verify_msannika_sources(bundle):
    """Verify every source used by the MS Annika importer."""
    if bundle.get("format") != "msannika_csm":
        return
    for receipt in bundle["source_receipts"]:
        raw = Path(receipt["path"]).read_bytes()
        if hashlib.sha256(raw).hexdigest() != receipt["sha256"]:
            raise ValueError("An MS Annika source checksum changed after parsing")


def read_msannika_csms(
    psm_path, fasta_path, *, dataset_accession,
    crosslinks_path=None, contaminant_fasta_path=None, backend="pyxlms",
):
    """Import XLSX CSMs while retaining raw rows and separate crosslink records."""
    import importlib
    from importlib.metadata import version
    from unittest.mock import patch

    import openpyxl
    import pandas as pd

    if backend not in {"native", "pyxlms"}:
        raise ValueError("backend must be native or pyxlms")
    if re.fullmatch(r"PXD\d{6,}", dataset_accession) is None:
        raise ValueError("A ProteomeXchange PXD accession is required")

    receipts = []

    def source_receipt(path, kind):
        path = Path(path).resolve()
        receipt = {
            "kind": kind, "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
        receipts.append(receipt)
        return receipt

    def table(path, required, kind):
        receipt = source_receipt(path, kind)
        workbook = openpyxl.load_workbook(
            path, read_only=True, data_only=False
        )
        candidates = []
        try:
            for sheet in workbook:
                nonempty = [
                    (number, cells)
                    for number, cells in enumerate(sheet.iter_rows(), 1)
                    if any(cell.value is not None for cell in cells)
                ]
                if not nonempty:
                    continue
                header_row, cells = nonempty[0]
                header = [cell.value for cell in cells]
                while header and header[-1] is None:
                    header.pop()
                if not required <= set(header):
                    continue
                if any(not isinstance(value, str) for value in header):
                    raise ValueError("Unsupported MS Annika column label")
                if len(header) != len(set(header)):
                    raise ValueError("Duplicate MS Annika column labels")
                records = []
                for number, cells in nonempty[1:]:
                    values = [cell.value for cell in cells]
                    if len(values) > len(header) and not any(
                        value is not None for value in values[len(header):]
                    ):
                        values = values[:len(header)]
                    fields = (
                        dict(zip(header, values))
                        if len(values) == len(header) else {}
                    )
                    records.append({
                        "source_record_id": (
                            f"{dataset_accession}:{receipt['sha256']}:"
                            f"{sheet.title}:{number:06d}"
                        ),
                        "source_record_number": len(records) + 1,
                        "source_worksheet": sheet.title,
                        "source_worksheet_row": number,
                        "source_values": values,
                        "source_fields": fields,
                        "source_cell_types": [
                            cell.data_type for cell in cells[:len(values)]
                        ],
                    })
                candidates.append((header, records))
        finally:
            workbook.close()
        if len(candidates) != 1:
            raise ValueError("Expected one uniquely identifiable MS Annika worksheet")
        return receipt, candidates[0]

    def integer(value, minimum):
        if isinstance(value, bool):
            raise ValueError("Boolean position")
        if isinstance(value, float):
            if not math.isfinite(value) or not value.is_integer():
                raise ValueError("Noninteger position")
            value = int(value)
        text = str(value).strip()
        if re.fullmatch(r"\d+", text) is None or int(text) < minimum:
            raise ValueError("Invalid position")
        return int(text)

    required = {
        "Sequence A", "Sequence B", "Crosslinker Position A",
        "Crosslinker Position B", "Accession A", "Accession B",
        "A in protein", "B in protein", "Alpha T/D", "Beta T/D",
        "Spectrum File", "First Scan", "Charge", "RT [min]", "Crosslink Strategy",
    }
    psm_receipt, (columns, records) = table(psm_path, required, "csms")
    fasta_receipt = source_receipt(fasta_path, "search_fasta")
    proteins = read_search_fasta(fasta_path)
    if contaminant_fasta_path is not None:
        source_receipt(contaminant_fasta_path, "contaminant_fasta")
        for identifier, protein in read_search_fasta(
            contaminant_fasta_path
        ).items():
            if identifier in proteins:
                if proteins[identifier]["sequence"] != protein["sequence"]:
                    raise ValueError("Conflicting search FASTA identifiers")
            else:
                proteins[identifier] = protein

    aliases = {}
    for identifier, protein in proteins.items():
        labels = [identifier]
        if protein["uniprot_accession"]:
            labels.append(protein["uniprot_accession"])
        parts = identifier.split("|")
        if (
            len(parts) == 3 and parts[0] == "sg"
            and ACCESSION.fullmatch(parts[1])
        ):
            labels.append(parts[1])
        for label in labels:
            aliases.setdefault(label, set()).add(identifier)

    crosslink_records = []
    crosslink_columns = []
    if crosslinks_path is not None:
        _, (crosslink_columns, crosslink_records) = table(
            crosslinks_path,
            {"Sequence A", "Sequence B", "Position A", "Position B",
             "Accession A", "Accession B", "In protein A", "In protein B"},
            "crosslinks",
        )

    def endpoint(fields, side):
        accession = str(fields.get("Accession " + side, "")).strip()
        peptide = fields.get("Sequence " + side)
        result = {
            "source_protein_id": accession,
            "source_fasta_id": None,
            "source_peptide": peptide,
            "source_modifications": fields.get("Modifications " + side),
            "peptide": peptide,
            "peptide_start": None,
            "link_position": None,
            "protein_position": None,
            "modifications": None,
            "status": "not_evaluable",
            "reason": "",
        }
        try:
            start_zero = integer(fields.get(side + " in protein"), 0)
            link = integer(fields.get("Crosslinker Position " + side), 1)
            if not isinstance(peptide, str) or re.fullmatch(r"[A-Z]+", peptide) is None:
                raise ValueError("Unsupported peptide sequence")
            if link > len(peptide):
                raise ValueError("Link position outside peptide")
        except (TypeError, ValueError):
            result["reason"] = "invalid_sequence_or_position"
            return result
        result.update(
            peptide_start=start_zero + 1,
            link_position=link,
            protein_position=start_zero + link,
        )
        matches = aliases.get(accession, set())
        if len(matches) != 1:
            result["reason"] = (
                "protein_identity_ambiguous" if matches
                else "protein_identifier_unresolved"
            )
            return result
        identifier = next(iter(matches))
        result["source_fasta_id"] = identifier
        if proteins[identifier]["sequence"][
            start_zero:start_zero + len(peptide)
        ] != peptide:
            result["reason"] = "peptide_sequence_mismatch"
            return result
        result.update(status="valid", reason="")
        return result

    for record in records:
        fields = record["source_fields"]
        flags = [
            str(fields.get(column, "")).strip().upper()
            for column in ("Alpha T/D", "Beta T/D")
        ]
        decoy = (
            "unknown" if any(flag not in {"T", "D"} for flag in flags)
            else "target_target" if flags == ["T", "T"]
            else "decoy_decoy" if flags == ["D", "D"]
            else "target_decoy"
        )
        endpoints = {
            site: endpoint(fields, side)
            for site, side in (("a", "A"), ("b", "B"))
        }
        reason = (
            "row_width_mismatch" if not fields
            else "unsupported_workbook_cell"
            if any(kind in {"f", "e"} for kind in record["source_cell_types"])
            else "source_engine_unresolved"
            if re.fullmatch(
                r"MS Annika(?:\s*\([^)]*\))?",
                str(fields.get("Crosslink Strategy", "")).strip(),
            ) is None
            else "invalid_decoy_flags" if decoy == "unknown"
            else "source_endpoint_unresolved"
            if any(item["status"] != "valid" for item in endpoints.values())
            else ""
        )
        record.update(
            record_kind="paired_crosslink_candidate",
            decoy_status=decoy,
            fdr={level: _fdr_value(None) for level in FDR_COLUMNS},
            endpoints=endpoints,
            record_status="not_evaluable" if reason else "valid",
            record_reason=reason,
            compensation_voltage_status=(
                "not_present_in_source" if "Compensation Voltage" not in columns
                else "reported_column"
            ),
        )

    bundle = {
        "format": "msannika_csm", "format_version": 1,
        "source_columns": columns, "records": records,
        "proteins": proteins,
        "crosslink_source_columns": crosslink_columns,
        "crosslink_records": crosslink_records,
        "source_receipts": receipts,
        "provenance": {
            "dataset_accession": dataset_accession,
            "psm_path": psm_receipt["path"],
            "psm_sha256": psm_receipt["sha256"],
            "fasta_path": fasta_receipt["path"],
            "fasta_sha256": fasta_receipt["sha256"],
        },
        "parser": {
            "name": "pyXLMS" if backend == "pyxlms" else "native",
            "version": version("pyXLMS") if backend == "pyxlms" else None,
            "parsed_records": 0,
            "modification_masses_assessed": False,
            "missing_compensation_voltage": (
                "unknown" if "Compensation Voltage" not in columns
                else "source_column"
            ),
        },
    }

    eligible = [
        record for record in records if record["record_status"] == "valid"
    ]
    if backend == "pyxlms" and eligible:
        from pyXLMS import parser

        marker = "__claudio_source_record_id"
        while marker in columns:
            marker += "_"
        frame = pd.DataFrame(
            [record["source_values"] for record in eligible], columns=columns
        )
        frame[marker] = [record["source_record_id"] for record in eligible]
        if "Compensation Voltage" not in frame:
            frame["Compensation Voltage"] = float("nan")
        else:
            frame["Compensation Voltage"] = frame["Compensation Voltage"].map(
                lambda value: float("nan") if value is None else value
            )
        module = importlib.import_module(parser.read_msannika.__module__)
        expected = {record["source_record_id"]: record for record in eligible}

        def parsed_endpoint(csm, prefix):
            return (
                csm[prefix + "_peptide"],
                csm[prefix + "_peptide_crosslink_position"],
                tuple(csm[prefix + "_proteins"]),
                tuple(csm[prefix + "_proteins_crosslink_positions"]),
                tuple(csm[prefix + "_proteins_peptide_positions"]),
                csm[prefix + "_decoy"],
            )

        def source_endpoint(record, site, side):
            value = record["endpoints"][site]
            return (
                value["peptide"], value["link_position"],
                (value["source_protein_id"],),
                (value["protein_position"],), (value["peptide_start"],),
                str(record["source_fields"][
                    "Alpha T/D" if side == "A" else "Beta T/D"
                ]).strip().upper() == "D",
            )

        try:
            with patch.object(module.pd, "read_excel", return_value=frame):
                parsed = parser.read_msannika(
                    str(psm_path), format="xlsx",
                    parse_modifications=False, unsafe=False, verbose=2,
                ).model_dump(mode="python")["crosslink_spectrum_matches"]
            if len(parsed) != len(eligible):
                raise ValueError("MS Annika parser record count differs from source")
            seen = set()
            for csm in parsed:
                source_id = csm["additional_information"]["source"][marker]
                if source_id not in expected or source_id in seen:
                    raise ValueError("Unknown or repeated MS Annika source association")
                seen.add(source_id)
                record = expected[source_id]
                actual = tuple(
                    parsed_endpoint(csm, prefix) for prefix in ("alpha", "beta")
                )
                native = (
                    source_endpoint(record, "a", "A"),
                    source_endpoint(record, "b", "B"),
                )
                direct, swapped = actual == native, actual[::-1] == native
                if not (direct or swapped):
                    raise ValueError("MS Annika parsed endpoints differ from source")
                fields = record["source_fields"]
                if (
                    csm["spectrum_file"] != str(fields["Spectrum File"]).strip()
                    or csm["scan_nr"] != integer(fields["First Scan"], 1)
                    or csm["charge"] != integer(fields["Charge"], 1)
                    or not math.isclose(
                        csm["retention_time"], float(fields["RT [min]"]) * 60,
                        rel_tol=1e-9, abs_tol=1e-6,
                    )
                ):
                    raise ValueError("MS Annika spectrum identity differs from source")
                if "Compensation Voltage" not in columns:
                    for key, value in csm.items():
                        if "compensation" in key.lower() or key.lower().endswith("_cv"):
                            if value is not None and not (
                                isinstance(value, float) and not math.isfinite(value)
                            ):
                                raise ValueError("Missing compensation voltage became a measurement")
                record["parser_orientation"] = (
                    "indistinguishable" if direct and swapped
                    else "same" if direct else "swapped"
                )
            bundle["parser"]["parsed_records"] = len(parsed)
        except Exception as error:
            for record in eligible:
                record.update(
                    record_status="not_evaluable",
                    record_reason="pyxlms_parse_failed",
                    parser_error=f"{type(error).__name__}: {error}",
                )

    _verify_msannika_sources(bundle)
    return bundle
def _msannika_author_filter_evidence(bundle, evidence_path):
    """Validate and package the author's description of the exported CSM table."""
    if evidence_path is None:
        raise ValueError("MS Annika export requires author filter evidence")
    path = Path(evidence_path).resolve()
    raw_policy = path.read_bytes()
    manifest_path = path.parent / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("Author evidence requires a checksum manifest")
    manifest = json.loads(manifest_path.read_bytes())
    expected = manifest.get("output_sha256", {}).get(path.name)
    if expected != hashlib.sha256(raw_policy).hexdigest():
        raise ValueError("Author filter evidence checksum changed")

    policy = json.loads(raw_policy)
    if (
        policy.get("kind") != "msannika_author_reported_filter_evidence"
        or policy.get("schema_version") != 1
        or policy.get("review_status") != "reviewed"
        or policy.get("dataset_accession")
        != bundle["provenance"]["dataset_accession"]
        or policy.get("reported_fdr_level") != "unspecified"
        or policy.get("per_record_fdr_values_assigned") is not False
        or policy.get("parent_confidence_propagated") is not False
        or policy.get("target_decoy_validation_performed") is not False
    ):
        raise ValueError("Unsupported MS Annika author evidence declaration")
    cutoff = policy.get("reported_fdr_cutoff")
    if (
        isinstance(cutoff, bool)
        or not isinstance(cutoff, (int, float))
        or not math.isfinite(cutoff)
        or not 0 <= cutoff <= 1
    ):
        raise ValueError("Invalid author-reported FDR cutoff")

    receipts = {item["kind"]: item for item in bundle["source_receipts"]}
    bindings = {
        "psm_sha256": "csms",
        "fasta_sha256": "search_fasta",
        "crosslinks_sha256": "crosslinks",
        "contaminant_fasta_sha256": "contaminant_fasta",
    }
    for field, kind in bindings.items():
        receipt = receipts.get(kind)
        declared = policy.get(field)
        if receipt is None:
            if declared is not None:
                raise ValueError("Author evidence requires the declared source: " + kind)
        elif declared != receipt["sha256"]:
            raise ValueError("Author evidence source checksum differs: " + kind)

    sources = policy.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("Author evidence has no supporting description")
    packaged = deepcopy(policy)
    packaged["sources"] = []
    bodies = {}
    supporting_hashes = []
    matched = False

    for receipt in sources:
        if receipt.get("kind") != "author_result_description":
            raise ValueError("Unsupported author evidence source kind")
        relative = Path(receipt["path"])
        source = (path.parent / relative).resolve()
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or path.parent not in source.parents
        ):
            raise ValueError("Author evidence source must remain inside its package")
        raw = source.read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if digest != receipt["sha256"]:
            raise ValueError("Author description checksum changed")
        line = policy.get("supporting_description_line")
        filename = policy.get("source_file_name")
        if (
            isinstance(line, str)
            and isinstance(filename, str)
            and line in raw.decode("utf-8-sig").splitlines()
            and line.split()
            and line.split()[0] == filename
            and "CSMs" in line
        ):
            reported = re.search(
                r"\((\d+(?:\.\d+)?)%\s*FDR\)", line, re.I
            )
            if reported and math.isclose(
                float(reported.group(1)) / 100, cutoff,
                rel_tol=0, abs_tol=1e-12,
            ):
                matched = True
        name = "author_filter_sources/" + digest + source.suffix
        bodies[name] = raw
        copied = deepcopy(receipt)
        copied.pop("absolute_path", None)
        copied["path"] = name
        packaged["sources"].append(copied)
        supporting_hashes.append(digest)

    if not matched:
        raise ValueError("Author description does not confirm the declared CSM table and cutoff")
    return packaged, bodies, sorted(set(supporting_hashes))


def write_msannika_claudio_input(
    bundle, mappings, output_dir, *, author_filter_path=None,
    allowed_crosslinkers, allowed_residues,
):
    """Export source-coordinate candidates with explicit author-reported confidence."""
    if bundle.get("format") != "msannika_csm":
        raise ValueError("Expected an MS Annika CSM bundle")
    if Path(output_dir).exists():
        raise FileExistsError("Output directory already exists")
    _verify_msannika_sources(bundle)
    policy, bodies, description_hashes = _msannika_author_filter_evidence(
        bundle, author_filter_path
    )

    selections = {}
    for name, values in (
        ("allowed_crosslinkers", allowed_crosslinkers),
        ("allowed_residues", allowed_residues),
    ):
        if (
            not isinstance(values, (set, frozenset, list, tuple))
            or not values
            or any(not isinstance(value, str) or not value.strip() for value in values)
        ):
            raise ValueError(name + " requires explicit nonempty labels")
        selections[name] = set(values)
    if any(
        re.fullmatch(r"[A-Z]", value) is None
        for value in selections["allowed_residues"]
    ):
        raise ValueError("Residues must be one-letter uppercase codes")

    receipts = {item["kind"]: item for item in bundle["source_receipts"]}
    native = read_msannika_csms(
        bundle["provenance"]["psm_path"],
        bundle["provenance"]["fasta_path"],
        dataset_accession=bundle["provenance"]["dataset_accession"],
        crosslinks_path=(
            receipts["crosslinks"]["path"] if "crosslinks" in receipts else None
        ),
        contaminant_fasta_path=(
            receipts["contaminant_fasta"]["path"]
            if "contaminant_fasta" in receipts else None
        ),
        backend="native",
    )
    if (
        native["source_columns"] != bundle["source_columns"]
        or native["proteins"] != bundle["proteins"]
        or native["crosslink_records"] != bundle["crosslink_records"]
        or len(native["records"]) != len(bundle["records"])
    ):
        raise ValueError("MS Annika bundle differs from its bound source files")
    for kind in ("psm", "fasta"):
        if native["provenance"][kind + "_sha256"] != bundle["provenance"][kind + "_sha256"]:
            raise ValueError("MS Annika provenance checksum differs from source")
    for original, checked in zip(bundle["records"], native["records"]):
        for field in (
            "source_record_id", "source_record_number", "source_worksheet",
            "source_worksheet_row", "source_values", "source_fields",
            "source_cell_types", "decoy_status", "fdr", "endpoints",
        ):
            if original[field] != checked[field]:
                raise ValueError("MS Annika source record changed after import")
        if original["record_status"] == "valid" and checked["record_status"] != "valid":
            raise ValueError("MS Annika source validation differs from its workbook")

    usable = {}
    for identifier, identity in mappings.items():
        source = bundle["proteins"].get(identifier)
        if not source or identity.get("status") != "mapped":
            continue
        accession = identity.get("accession")
        if (
            not isinstance(accession, str)
            or ACCESSION.fullmatch(accession) is None
            or identity.get("sequence") != source["sequence"]
        ):
            raise ValueError("MS Annika mapped identity differs from the source sequence")
        declared = source.get("uniprot_accession")
        parts = identifier.split("|")
        if (
            declared is None and len(parts) == 3 and parts[0] == "sg"
            and ACCESSION.fullmatch(parts[1])
        ):
            declared = parts[1]
        if declared is None or accession != declared:
            raise ValueError("Exact identity requires the declared search accession alias")
        _reference_evidence(
            {
                "accession": accession,
                "taxon_id": identity.get("taxon_id"),
                "reference_sequence": identity["sequence"],
            },
            identity,
            source,
            require_declared_accession=False,
        )
        usable[identifier] = identity

    rows, audit = [], []
    for record in bundle["records"]:
        entry = deepcopy(record)
        reasons = []
        if record["record_status"] != "valid":
            reasons.append(record["record_reason"] or "source_record_unresolved")
        if record["decoy_status"] == "unknown":
            reasons.append("invalid_decoy_flags")
        elif record["decoy_status"] != "target_target":
            reasons.append("decoy_record")
        identities = {}
        for site in ("a", "b"):
            endpoint = record["endpoints"][site]
            identity = usable.get(endpoint["source_fasta_id"])
            if identity is None:
                reasons.append("protein_identity_unresolved_" + site)
            else:
                identities[site] = identity

        fields = record["source_fields"]
        if not reasons:
            if fields.get("Crosslinker") not in selections["allowed_crosslinkers"]:
                reasons.append("crosslinker_not_selected")
            for site in ("a", "b"):
                endpoint = record["endpoints"][site]
                residue = endpoint["peptide"][endpoint["link_position"] - 1]
                if residue not in selections["allowed_residues"]:
                    reasons.append("residue_not_selected_" + site)

        entry.update(
            input_status="not_evaluable" if reasons else "ready",
            input_reasons=list(dict.fromkeys(reasons)),
            fdr_verification_status="author_reported_only",
            source_reported_fdr_level="unspecified",
            source_reported_fdr_cutoff=policy["reported_fdr_cutoff"],
        )
        audit.append(entry)
        if reasons:
            continue

        row = {
            "source_record_id": record["source_record_id"],
            "identity_scope": "exact_sequence",
            "input_scope": "author_reported_candidate",
            "fdr_verification_status": "author_reported_only",
            "source_reported_fdr_level": "unspecified",
            "source_reported_fdr_cutoff": policy["reported_fdr_cutoff"],
            "selected_fdr_level": None,
            "selected_fdr_value": None,
            "selected_fdr_cutoff": None,
            "source_decoy_status": record["decoy_status"],
            "source_crosslinker": fields["Crosslinker"],
            "source_spectrum_file": fields["Spectrum File"],
            "source_scan": fields["First Scan"],
            "source_worksheet": record["source_worksheet"],
            "source_worksheet_row": record["source_worksheet_row"],
            "source_filter_description_sha256_json": json.dumps(description_hashes),
            "source_fields_json": json.dumps(fields, ensure_ascii=False),
        }
        for site, side in (("a", "A"), ("b", "B")):
            endpoint, identity = record["endpoints"][site], identities[site]
            row.update({
                "unip_id_" + site: identity["accession"],
                "seq_" + site: identity["sequence"],
                "pep_" + site: endpoint["peptide"],
                "pos_" + site: endpoint["protein_position"],
                "res_pos_" + site: endpoint["link_position"],
                "source_protein_id_" + site: endpoint["source_protein_id"],
                "source_fasta_id_" + site: endpoint["source_fasta_id"],
                "source_peptide_start_zero_based_" + site: fields[side + " in protein"],
                "source_modifications_" + site: fields.get("Modifications " + side),
            })
        rows.append(row)

    def encoded(value):
        return (
            json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        ).encode("utf-8")

    archived_receipts = []
    for receipt in bundle["source_receipts"]:
        path = Path(receipt["path"])
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != receipt["sha256"]:
            raise ValueError("MS Annika source checksum changed before output")
        name = "input_sources/" + receipt["sha256"] + path.suffix
        bodies[name] = raw
        archived_receipts.append({
            **receipt, "original_path": receipt["path"], "path": name,
        })
    bodies["input_source_receipts.json"] = encoded(archived_receipts)
    bodies["crosslink_audit.json"] = encoded(bundle["crosslink_records"])
    bodies["crosslink_source_columns.json"] = encoded(
        bundle["crosslink_source_columns"]
    )
    bodies["author_filter_evidence.json"] = encoded(policy)

    prepared = {
        "rows": rows,
        "audit": audit,
        "provenance": deepcopy(bundle["provenance"]),
        "fdr_policy": {
            "kind": "author_reported_source_table",
            "level": "unspecified",
            "reported_cutoff": policy["reported_fdr_cutoff"],
            "verification_status": "author_reported_only",
            "per_record_fdr_applied": False,
            "application_verified": False,
            "target_decoy_validation_performed": False,
            "parent_confidence_propagated": False,
        },
        "author_filter_policy": policy,
        "author_filter_evidence_snapshot": "author_filter_evidence.json",
        "extra_source_bodies": bodies,
    }
    return write_xifdr_outputs(
        bundle, mappings, prepared, output_dir,
        pipeline_policy={
            "allowed_crosslinkers": sorted(selections["allowed_crosslinkers"]),
            "allowed_residues": sorted(selections["allowed_residues"]),
            "identity_scope": "exact_sequence",
            "source_csm_peptide_start_basis": 0,
            "output_coordinate_basis": 1,
            "requires_parent_crosslink_match": False,
            "modification_masses_assessed": False,
        },
    )