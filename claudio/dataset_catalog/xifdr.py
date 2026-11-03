"""Import xiFDR PSMs, verify source coordinates, and prepare existing AF analysis.

Coordinates are 1-based. Identity resolution requires full sequence equality.
FDR levels and source classifications are preserved rather than inferred.
"""

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
        match = re.match(r"[a-z][a-z0-9]*", raw[position:])
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


def project_xifdr_candidates(bundle, exact_mappings, candidate_reports, *, max_optimal_alignments=64):
    """Project intact peptides onto each reference hypothesis; retain every PSM.

    A peptide projection must have one distinct, contiguous, sequence-identical
    position among the compatible mappings. Multiple protein identities stay
    separate hypotheses. This stage neither filters FDR nor accepts identities.
    The alignment limit bounds computation; exceeding it leaves an unresolved
    audit record rather than selecting the first alignment.
    Candidate reports may be a dictionary or snapshots in retrieval order;
    the last completed result is reused while failed attempts stay in evidence.
    """
    if type(max_optimal_alignments) is not int or max_optimal_alignments < 1:
        raise ValueError("max_optimal_alignments must be a positive integer")
    candidate_reports, report_indices, report_history = _selected_candidate_reports(candidate_reports)
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
        if source_record["record_status"] != "valid" or source_record["decoy_status"] != "target_target":
            row["reason"] = source_record["record_reason"] or "invalid_decoy_flags"
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
            "projection_rule_version": "intact-peptide-projection-v1", "max_optimal_alignments": max_optimal_alignments,
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


def prepare_af_input(bundle, mappings, *, fdr_level, fdr_cutoff):
    """Project validated target PSMs onto exact UniProt identities, retaining all audit rows."""
    policy = _fdr_policy(fdr_level, fdr_cutoff)
    rows, audit = [], []
    for record in bundle["records"]:
        entry, reasons = deepcopy(record), []
        if record["decoy_status"] == "unknown":
            reasons.append("invalid_decoy_flags")
        elif record["decoy_status"] != "target_target":
            reasons.append("decoy_record")
        if record["record_status"] != "valid":
            reasons.append(record["record_reason"])
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
        fields = record["source_fields"]
        row = {"source_record_id": record["source_record_id"], "dataset_accession": bundle["provenance"]["dataset_accession"],
               "source_psm_sha256": bundle["provenance"]["psm_sha256"], "source_fasta_sha256": bundle["provenance"]["fasta_sha256"],
               "source_record_number": record["source_record_number"], "source_psm_id": fields.get("PSMID", ""),
               "source_run": fields.get("run", ""), "source_scan": fields.get("scan", ""),
               "source_crosslinker": fields.get("Crosslinker", ""), "source_fdr_group": fields.get("fdrGroup", ""),
               "source_fields_json": json.dumps(fields, ensure_ascii=False),
               "selected_fdr_level": fdr_level, "selected_fdr_value": confidence["value"], "selected_fdr_cutoff": float(fdr_cutoff)}
        for site in ("a", "b"):
            endpoint, identity = record["endpoints"][site], selected[site]
            row.update({f"unip_id_{site}": identity["accession"], f"seq_{site}": identity["sequence"],
                        f"pep_{site}": endpoint["peptide"], f"pos_{site}": endpoint["protein_position"],
                        f"res_pos_{site}": endpoint["link_position"], f"position_status_{site}": "valid",
                        f"source_protein_id_{site}": endpoint["source_protein_id"],
                        f"source_protein_position_{site}": endpoint["protein_position"],
                        f"source_modified_peptide_{site}": endpoint["source_peptide"]})
        rows.append(row)
    return {"rows": rows, "audit": audit, "fdr_policy": policy, "provenance": deepcopy(bundle["provenance"])}


def write_claudio_input(bundle, mappings, output_dir, *, fdr_level, fdr_cutoff,
                        allowed_crosslinkers, allowed_residues, source_protein_ids=None):
    """Export exact-identity PSMs for the existing M1 -> M3 -> M2 -> M4 CLI.

    The CSV uses M1's default input column names.
    Selection is explicit and every original PSM remains in the audit. The
    original PSM key is named xlms_source_record_id because M1 creates its own
    source_record_id. Use the AF pipeline's source-preserving preprocessing.
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

    prepared = prepare_af_input(bundle, mappings, fdr_level=fdr_level, fdr_cutoff=fdr_cutoff)
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
    policy = {key: sorted(values) if values is not None else None for key, values in selections.items()}
    return write_xifdr_outputs(bundle, mappings, prepared, output_dir, pipeline_policy=policy)


def write_xifdr_outputs(bundle, mappings, prepared, output_dir, *, pipeline_policy=None):
    """Write source/audit snapshots and the eligible input table to a new directory."""
    output_dir = Path(output_dir)
    copied = deepcopy(mappings)
    source_bodies = {}
    for identity in copied.values():
        for source in identity["sources"]:
            raw = Path(source["absolute_path"]).read_bytes()
            if hashlib.sha256(raw).hexdigest() != source["sha256"]:
                raise ValueError("A protein-identity source changed after lookup")
            name = "identity_sources/" + source["sha256"] + Path(source["absolute_path"]).suffix
            source_bodies[name] = raw
            source["original_path"], source["path"] = source.pop("absolute_path"), name
    output_dir.mkdir(parents=True, exist_ok=False)
    if source_bodies:
        (output_dir / "identity_sources").mkdir()
        for name, raw in source_bodies.items():
            (output_dir / name).write_bytes(raw)
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
    names = ["psm_audit.json", "search_proteins.json", "protein_mapping.json", "source_columns.json", input_name, *source_bodies]
    manifest = {"schema_version": 1, "kind": "xifdr_claudio_input" if pipeline else "xifdr_prepared_input", "status": "complete",
                "provenance": bundle["provenance"], "fdr_policy": prepared["fdr_policy"],
                "counts": {"source_records": len(bundle["records"]),
                           "pipeline_input_rows" if pipeline else "af_input_rows": len(rows),
                           "protein_mapping_statuses": dict(Counter(m["status"] for m in mappings.values()))},
                "assessment_performed": False,
                "output_sha256": {name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest() for name in names}}
    if "parser" in bundle:
        manifest["parser"] = deepcopy(bundle["parser"])
    if pipeline:
        manifest["selection_policy"] = pipeline_policy
        manifest["pipeline_input"] = {
            "path": input_name, "projections": ",".join(core_columns),
            "source_key": "xlms_source_record_id", "requires_source_preservation": True,
        }
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
