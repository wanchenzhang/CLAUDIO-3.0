"""Command-line tools for discovering public XL-MS dataset candidates."""

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import click
import requests
import zipfile
from .inventory import inspect_catalog, summarize_inventories, inspect_content_plan
from .proteomecentral import fetch_dataset_pages, write_discovery_outputs
from .xifdr import FDR_COLUMNS, read_xifdr_psms, read_plink_psms, write_claudio_input, write_import_report, read_merox_csms, write_merox_claudio_input, read_scout_csms

DEFAULT_SEARCH_TERMS = (
    "XL-MS",
    "crosslinking",
    "cross-linking",
    "cross linking",
    "crosslinked",
)


@click.group()
def main():
    """Manage the public dataset candidate catalog."""


@main.command()
@click.option(
    "--search-term", "search_terms",
    multiple=True,
    default=DEFAULT_SEARCH_TERMS,
    show_default=True,
    help="Override the built-in search terms. Repeat for multiple terms.",
)
@click.option(
    "--output", "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="New directory for the discovery snapshot.",
)
@click.option(
    "--page-size",
    type=click.IntRange(min=1),
    default=100,
    show_default=True,
    help="Number of records requested per page.",
)
def discover(search_terms, output_dir, page_size):
    """Discover candidate datasets and preserve their search provenance."""
    search_terms = tuple(term.strip() for term in search_terms)
    if any(not term for term in search_terms):
        raise click.BadParameter(
            "Search terms must not be empty.", param_hint="--search-term"
        )
    search_terms = tuple(dict.fromkeys(search_terms))

    if output_dir.exists():
        raise click.BadParameter(
            "Choose a new output directory.", param_hint="--output"
        )

    searches = {}
    try:
        for term in search_terms:
            pages = fetch_dataset_pages(term, page_size=page_size)
            searches[term] = pages
            count = sum(len(page.data["datasets"]) for page in pages)
            click.echo(f"{term}: {count} records")

        destination = write_discovery_outputs(searches, output_dir)
        manifest = json.loads(
            (destination / "manifest.json").read_text(encoding="utf-8")
        )
    except (requests.exceptions.RequestException, ValueError, OSError) as error:
        raise click.ClickException(str(error)) from error

    click.echo(f"Candidate datasets: {manifest['candidate_count']}")
    click.echo(f"Output: {destination.resolve()}")


def _inventory_progress(completed, total, result):
    if result is None:
        click.echo(f"Saved progress: {completed}/{total}")
        return

    status = result["inspection_status"]
    message = (
        f"[{completed}/{total}] {result['accession']} "
        f"({result.get('repository', 'unknown')}): {status}"
    )
    if status == "collected":
        message += f"; {result['file_count']} file records"
    elif status == "failed":
        # Keep each progress record on one log line.
        error = " ".join(result["error_message"].splitlines())
        message += f"; {result['error_stage']}: {error}"
    click.echo(message)


def _content_progress(completed, total, result):
    if result is None:
        click.echo(f"Saved progress: {completed}/{total}")
        return
    message = (f"[{completed}/{total}] {result['pilot_file_id']} {result['project_accession']}: "
               f"{result['download_status']}; content={result['content_probe']['status']}")
    if result.get("error_message"):
        message += "; " + " ".join(result["error_message"].splitlines())
    click.echo(message)


@main.command(name="inspect")
@click.option(
    "--catalog", "catalog_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Discovery candidates.json for metadata collection; exclusive with --plan.",
)
@click.option(
    "--plan", "plan_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Validated content pilot JSON; downloads only its selected files within its limits.",
)
@click.option(
    "--output", "output_dir",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="New inventory directory, or the existing run with --resume.",
)
@click.option(
    "--resume", is_flag=True,
    help="Validate saved inventories and continue after the last recorded result. "
         "Recorded failures are retained, not retried.",
)
@click.option(
    "--delay", type=float, default=1.0, show_default=True,
    help="Seconds between projects or planned files; must be finite and non-negative.",
)
def inspect_command(catalog_path, plan_path, output_dir, resume, delay):
    """Inspect catalog metadata or download and probe a size-screened content plan.

    Catalog mode collects metadata only. Plan mode downloads its selected files;
    content findings remain separate from scientific eligibility decisions.
    """
    if (catalog_path is None) == (plan_path is None):
        raise click.UsageError("Provide exactly one of --catalog or --plan.")
    if not math.isfinite(delay) or delay < 0:
        raise click.BadParameter(
            "Must be finite and non-negative.", param_hint="--delay"
        )
    checkpoint_name = "content_checkpoint.json" if plan_path is not None else "checkpoint.json"
    if resume:
        if not (output_dir / checkpoint_name).is_file():
            raise click.BadParameter(
                f"No {checkpoint_name} found for --resume.", param_hint="--output"
            )
    elif output_dir.exists():
        raise click.BadParameter(
            "Choose a new directory, or use --resume for an existing run.",
            param_hint="--output",
        )

    if plan_path is not None:
        click.echo("Validating content plan and saved sources...")
        try:
            results = inspect_content_plan(
                plan_path, output_dir, resume=resume, delay=delay,
                on_progress=_content_progress,
            )
        except (ValueError, OSError, requests.exceptions.RequestException) as error:
            raise click.ClickException(str(error)) from error
        counts = Counter(row["download_status"] for row in results)
        click.echo(f"Planned files: {len(results)}")
        click.echo(f"Downloaded: {counts['downloaded']}")
        click.echo(f"Failed: {counts['failed']}")
        click.echo(f"Budget exhausted: {counts['budget_exhausted']}")
        click.echo(f"Skipped unknown size: {counts['skipped_unknown_size']}")
        click.echo(f"Skipped file limit: {counts['skipped_file_limit']}")
        backends = sorted({row.get("download_backend", "legacy") for row in results})
        click.echo(f"Download backend: {', '.join(backends)}")
        if "pridepy" in backends:
            click.echo("Limits apply to planned file sizes, not network traffic")
        for status, count in Counter(row["content_probe"]["status"] for row in results).items():
            click.echo(f"Content {status}: {count}")
        click.echo("Scientific review: unchanged")
        click.echo(f"Output: {output_dir.resolve()}")
        return

    click.echo("Validating saved inventory..." if resume else "Starting inventory...")
    try:
        results = inspect_catalog(
            catalog_path, output_dir, resume=resume, delay=delay,
            on_progress=_inventory_progress,
        )
    except (requests.exceptions.RequestException, ValueError, OSError) as error:
        raise click.ClickException(str(error)) from error

    counts = Counter(result["inspection_status"] for result in results)
    click.echo(f"Recorded outcomes: {len(results)}")
    click.echo(f"Collected: {counts['collected']}")
    click.echo(f"Failed: {counts['failed']}")
    click.echo(f"Adapter pending: {counts['adapter_pending']}")
    click.echo("Scientific review: unchanged")
    click.echo(f"Output: {output_dir.resolve()}")


@main.command(name="summarize")
@click.option(
    "--catalog", "catalog_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path), required=True,
    help="Full discovery candidates.json defining the summary population.",
)
@click.option(
    "--run", "run_dirs", multiple=True, required=True,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Finished inventory run. Repeat in chronological order; the last "
         "successful result is preferred over failures.",
)
@click.option(
    "--output", "output_dir", required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="New directory for the summary and complete attempt history.",
)
def summarize_command(catalog_path, run_dirs, output_dir):
    """Summarize saved runs without downloading data or changing source files."""
    click.echo("Validating saved runs and raw response checksums...")
    try:
        destination = summarize_inventories(catalog_path, run_dirs, output_dir)
        manifest = json.loads((destination / "manifest.json").read_bytes())
    except (ValueError, OSError) as error:
        raise click.ClickException(str(error)) from error
    counts = manifest["status_counts"]
    click.echo(f"Candidate datasets: {manifest['candidate_count']}")
    click.echo(f"Recorded attempts: {manifest['attempt_count']}")
    click.echo(f"Collected: {counts['collected']}")
    click.echo(f"Failed: {counts['failed']}")
    click.echo(f"Adapter pending: {counts['adapter_pending']}")
    screening = manifest["file_screening"]
    click.echo(f"File records screened: {screening['file_count']}")
    click.echo(f"Priority download candidates: {screening['download_candidate_count']}")
    click.echo("Scientific review: unchanged")
    click.echo(f"Output: {destination}")


def _read_protein_mapping(mapping_path):
    """Load live lookup results or an exported mapping with copied evidence."""
    try:
        mappings = json.loads(mapping_path.read_text(encoding="utf-8-sig"))
    except ValueError as error:
        raise ValueError(f"Invalid protein mapping JSON: {error}") from error
    if not isinstance(mappings, dict):
        raise ValueError("Protein mapping must be a JSON object keyed by search FASTA identifier.")
    for identifier, identity in mappings.items():
        if (not isinstance(identity, dict) or not isinstance(identity.get("status"), str)
                or not isinstance(identity.get("sources"), list)):
            raise ValueError(f"Invalid protein mapping for {identifier}.")
        for source in identity["sources"]:
            if not isinstance(source, dict) or not isinstance(source.get("sha256"), str):
                raise ValueError(f"Invalid evidence source in protein mapping for {identifier}.")
            saved_path = source.get("path")
            saved = mapping_path.parent / saved_path if isinstance(saved_path, str) and saved_path else None
            original_path = source.get("absolute_path")
            if saved is not None and saved.is_file():
                path = saved
            elif isinstance(original_path, str) and original_path:
                path = Path(original_path)
            else:
                raise ValueError(f"Protein mapping evidence file is missing for {identifier}.")
            if not path.is_file():
                raise ValueError(f"Protein mapping evidence file is missing: {path}")
            source["absolute_path"] = str(path.resolve())
    return mappings


def _read_reference_projection(path):
    """Load a saved projection and resolve live or copied evidence paths."""
    path = path.resolve()
    projection = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(projection, dict):
        raise ValueError("Reference projection must be a JSON object")
    manifest_path = path.parent / "manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        expected = manifest.get("output_sha256", {}).get(path.name)
        if expected is None or hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError("Reference projection checksum does not match its saved manifest")

    def resolve_sources(value):
        if isinstance(value, list):
            for item in value:
                resolve_sources(item)
        elif isinstance(value, dict):
            if "kind" in value and "sha256" in value and "path" in value:
                saved = path.parent / value["path"]
                original = value.get("absolute_path")
                source = saved if saved.is_file() else Path(original) if original else None
                if source is None or not source.is_file():
                    raise ValueError("Reference projection evidence file is missing")
                value["absolute_path"] = str(source.resolve())
            for child in value.values():
                resolve_sources(child)
    resolve_sources(projection)
    return projection


def _detect_psm_format(results_path):
    """Recognize supported tables and native MeroX containers."""
    if zipfile.is_zipfile(results_path):
        with zipfile.ZipFile(results_path) as archive:
            names = archive.namelist()
        if (
            names.count("Result.csv") == 1
            and names.count("properties.mxf") == 1
        ):
            return "merox_csm"
        raise ValueError(
            "ZIP container is not a supported native MeroX result"
        )

    with results_path.open(encoding="utf-8-sig", newline="") as handle:
        columns = set(next(csv.reader(handle), []))

    signatures = {
        "xifdr_psm": {
            "PSMID", "PepSeq1", "PepSeq2",
            "ProteinLinkPos1", "ProteinLinkPos2",
        },
        "plink_psm": {
            "Title", "Peptide", "Proteins", "Peptide_Type", "Linker",
        },
        "scout_csm": {
            "Scan", "Link-Type", "Alpha peptide", "Beta peptide",
            "Alpha protein(s) position(s)",
            "Beta protein(s) position(s)",
            "IsDecoy", "File",
        },
    }
    matches = [
        name for name, required in signatures.items()
        if required <= columns
    ]
    if len(matches) != 1:
        raise ValueError(
            "Cannot identify a unique results format; supported CSV layouts "
            "are xiFDR PSM, pLink PSM and Scout CSM tables."
        )
    return matches[0]


@main.command(name="standardize")
@click.option(
    "--results", "results_path", required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="xiFDR/pLink/Scout CSM CSV or a supported native MeroX container.",
)
@click.option(
    "--fasta", "fasta_path", required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="FASTA used for the original search.",
)
@click.option(
    "--mapping", "mapping_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Saved protein_mapping.json with its identity evidence files; required for input export.",
)
@click.option(
    "--reference-projection", "reference_projection_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
        help=(
        "Saved reference projection report, or a MeroX reference "
        "mapping with its saved UniProt evidence."
    ),
)
@click.option("--dataset-accession", required=True, help="ProteomeXchange PXD accession.")
@click.option(
    "--fdr-level", type=click.Choice(tuple(FDR_COLUMNS)),
    help="Per-record FDR level; required for export unless using --author-filter-evidence.",
)
@click.option(
    "--fdr-cutoff", type=click.FloatRange(min=0.0, max=1.0),
    help="Inclusive per-record FDR cutoff as a fraction. 0.01 means 1%; omit with --author-filter-evidence.",
)
@click.option(
    "--fdr-column",
    help="Manual pLink column corresponding to --fdr-level; use with a declared --decoy-prefix.",
)
@click.option(
    "--decoy-prefix",
    help="Manual pLink decoy prefix; use with a reported --fdr-column.",
)
@click.option(
    "--plink-unfiltered", "plink_unfiltered_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Native pLink 2 unfiltered PSM CSV with Q-value and Target_Decoy; pair with --plink-parameters.",
)
@click.option(
    "--plink-parameters", "plink_parameters_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Original .plink search parameters; pair with --plink-unfiltered.",
)
@click.option(
    "--author-filter-evidence", "author_filter_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Reviewed declaration for a target-only author-filtered pLink PSM table; preserves missing per-record FDR.",
)
@click.option(
    "--crosslinker", "crosslinkers", multiple=True,
    help="Crosslinker label to include; required for export. Repeat for multiple labels.",
)
@click.option(
    "--residue", "residues", multiple=True,
    type=click.Choice(tuple("ABCDEFGHIJKLMNOPQRSTUVWXYZ")),
    help="Residue letter to include at both sites; required for export. Repeat for multiple letters.",
)
@click.option(
    "--audit-only", is_flag=True,
    help="Save every imported record without identity mapping, filtering, or input export.",
)
@click.option(
    "--source-protein", "source_proteins", multiple=True,
    help="Search FASTA identifier to include; omit to include all resolved proteins.",
)
@click.option(
    "--output", "output_dir", required=True,
    type=click.Path(file_okay=False, path_type=Path),
    help="New directory for the input export or import audit.",
)
def standardize_command(results_path, fasta_path, mapping_path, dataset_accession,
                        fdr_level, fdr_cutoff, fdr_column, decoy_prefix, crosslinkers,
                        residues, source_proteins, audit_only, output_dir,
                        plink_unfiltered_path, plink_parameters_path, author_filter_path,
                        reference_projection_path=None):
    """Standardize supported XL-MS results or preserve a full import audit."""
    if output_dir.exists():
        raise click.BadParameter("Output directory already exists; choose a new directory.", param_hint="--output")
    if fdr_cutoff is not None and not math.isfinite(fdr_cutoff):
        raise click.BadParameter("Must be a finite fraction between 0 and 1.", param_hint="--fdr-cutoff")
    if any(not label.strip() for label in crosslinkers):
        raise click.BadParameter("Crosslinker labels must not be empty.", param_hint="--crosslinker")
    if any(not identifier.strip() for identifier in source_proteins):
        raise click.BadParameter("Protein identifiers must not be empty.", param_hint="--source-protein")
    try:
        source_format = _detect_psm_format(results_path)
    except (ValueError, OSError) as error:
        raise click.ClickException(str(error)) from error
    if source_format == "scout_csm":
        if not audit_only:
            raise click.UsageError(
                "Scout input currently requires --audit-only."
            )
        if len(crosslinkers) != 1:
            raise click.UsageError(
                "Scout import requires exactly one --crosslinker label."
            )

        unsupported = (
            mapping_path, reference_projection_path,
            fdr_level, fdr_cutoff, fdr_column, decoy_prefix,
            plink_unfiltered_path, plink_parameters_path,
            author_filter_path,
        )
        if (
            any(value is not None for value in unsupported)
            or residues
            or source_proteins
        ):
            raise click.UsageError(
                "Scout --audit-only cannot be combined with "
                "identity, FDR, residue, source-protein "
                "or pLink evidence selection options."
            )

        click.echo("Detected format: Scout CSM")
        try:
            bundle = read_scout_csms(
                results_path,
                fasta_path,
                dataset_accession=dataset_accession,
                crosslinker=crosslinkers[0],
            )
            destination = write_import_report(bundle, output_dir)
            manifest = json.loads(
                (destination / "manifest.json").read_bytes()
            )
        except (
            ValueError, OSError, RuntimeError,
            KeyError, TypeError, csv.Error,
        ) as error:
            raise click.ClickException(str(error)) from error

        counts = manifest["counts"]
        loop_links = sum(
            record["record_kind"] == "loop_link"
            for record in bundle["records"]
        )
        click.echo(f"Source CSM records: {counts['source_records']}")
        click.echo(
            f"Parsed CSM records: {bundle['parser']['parsed_records']}"
        )
        click.echo(
            f"Validated source records: {counts['validated_source_records']}"
        )
        click.echo(f"Loop-link records retained: {loop_links}")
        click.echo(f"Records with parser errors: {counts['parser_errors']}")
        click.echo(
            f"Crosslinker label: {crosslinkers[0]} (caller-declared)"
        )
        click.echo(
            "UniProt mapping, FDR filtering and input export: not performed"
        )
        click.echo(
            f"Audit report: "
            f"{(destination / 'import_report.json').resolve()}"
        )
        click.echo(f"Output: {destination.resolve()}")
        return
    if source_format == "merox_csm":
        if audit_only:
            policies = (
                mapping_path, reference_projection_path,
                fdr_level, fdr_cutoff, fdr_column, decoy_prefix,
                plink_unfiltered_path, plink_parameters_path,
                author_filter_path,
            )
            if (
                any(value is not None for value in policies)
                or crosslinkers
                or residues
                or source_proteins
            ):
                raise click.UsageError(
                    "MeroX --audit-only cannot be combined with "
                    "selection or pLink evidence options."
                )
        else:
            missing = [
                name
                for name, absent in (
                    (
                        "--reference-projection",
                        reference_projection_path is None,
                    ),
                    ("--crosslinker", not crosslinkers),
                    ("--residue", not residues),
                )
                if absent
            ]
            if missing:
                raise click.UsageError(
                    "MeroX export requires "
                    + ", ".join(missing)
                    + ". Use --audit-only to save an import audit."
                )

            unsupported = (
                mapping_path, fdr_level, fdr_cutoff,
                fdr_column, decoy_prefix, plink_unfiltered_path,
                plink_parameters_path, author_filter_path,
            )
            if (
                any(value is not None for value in unsupported)
                or source_proteins
            ):
                raise click.UsageError(
                    "MeroX export uses saved native confidence and "
                    "one explicit reference mapping. Do not combine it "
                    "with --mapping, per-record FDR, pLink evidence "
                    "or --source-protein options."
                )

        click.echo("Detected format: MeroX CSM (native container)")
        try:
            bundle = read_merox_csms(
                results_path,
                fasta_path,
                dataset_accession=dataset_accession,
            )
            if audit_only:
                destination = write_import_report(bundle, output_dir)
            else:
                reference_mapping = json.loads(
                    reference_projection_path.read_text(
                        encoding="utf-8-sig"
                    )
                )
                if not isinstance(reference_mapping, dict):
                    raise ValueError(
                        "MeroX reference mapping must be a JSON object"
                    )
                destination = write_merox_claudio_input(
                    bundle,
                    reference_mapping,
                    output_dir,
                    allowed_crosslinkers=crosslinkers,
                    allowed_residues=residues,
                )

            manifest = json.loads(
                (destination / "manifest.json").read_bytes()
            )
        except (
            ValueError, OSError, RuntimeError, KeyError, TypeError,
            zipfile.BadZipFile, csv.Error,
        ) as error:
            raise click.ClickException(str(error)) from error

        failures = sum(
            record.get("record_reason", "").startswith("pyxlms_")
            for record in bundle["records"]
        )
        click.echo(
            f"Source CSM records: {manifest['counts']['source_records']}"
        )
        click.echo(
            f"Parsed CSM records: {bundle['parser']['parsed_records']}"
        )
        click.echo(f"Records with parser errors: {failures}")

        if audit_only:
            click.echo(
                "Scientific filtering, identity mapping and "
                "input export: not performed"
            )
            click.echo(
                f"Audit report: "
                f"{(destination / 'import_report.json').resolve()}"
            )
        else:
            click.echo(
                f"Standardized input rows: "
                f"{manifest['counts']['pipeline_input_rows']}"
            )
            click.echo(
                "Localization: source-reported best; "
                "all source candidates retained in the audit"
            )
            click.echo("Native FDR scope: unspecified")
            click.echo("Model assessment: not performed")
            click.echo(
                f"Input: {(destination / 'claudio_input.csv').resolve()}"
            )
            click.echo(
                f"Audit: {(destination / 'psm_audit.json').resolve()}"
            )

        click.echo(f"Output: {destination.resolve()}")
        return
    plink = source_format == "plink_psm"
    native = plink_unfiltered_path is not None or plink_parameters_path is not None
    author = author_filter_path is not None
    if not plink and (fdr_column is not None or decoy_prefix is not None or native or author):
        raise click.UsageError("Decoy, FDR-column and native pLink evidence options only support pLink PSM inputs.")
    if author and (native or fdr_column is not None or decoy_prefix is not None
                   or fdr_level is not None or fdr_cutoff is not None):
        raise click.UsageError("--author-filter-evidence cannot be mixed with native or per-record FDR options.")
    if (plink_unfiltered_path is None) != (plink_parameters_path is None):
        raise click.UsageError("Provide --plink-unfiltered and --plink-parameters together.")
    if native and (fdr_column is not None or decoy_prefix is not None):
        raise click.UsageError("Do not mix native pLink evidence with --fdr-column or --decoy-prefix.")
    if audit_only:
        selection = {"--mapping": mapping_path is not None, "--fdr-cutoff": fdr_cutoff is not None,
                     "--crosslinker": bool(crosslinkers), "--residue": bool(residues),
                     "--source-protein": bool(source_proteins),
                     "--reference-projection": reference_projection_path is not None}
        conflicts = [option for option, present in selection.items() if present]
        if conflicts:
            raise click.UsageError("Selection options cannot be used with --audit-only: " + ", ".join(conflicts))
        if native and fdr_level is not None:
            raise click.UsageError("Native pLink audit retains PSM confidence; omit --fdr-level with --audit-only.")
        if plink and not native and not author and (fdr_level is None) != (fdr_column is None):
            raise click.UsageError("Provide --fdr-level and --fdr-column together to label reported pLink FDR values.")
        if not plink and fdr_level is not None:
            raise click.UsageError("--fdr-level is not used with --audit-only for xiFDR; all reported levels are retained.")
    else:
        required = {"--mapping": mapping_path is None, "--crosslinker": not crosslinkers,
                    "--residue": not residues}
        if not author:
            required.update({"--fdr-level": fdr_level is None, "--fdr-cutoff": fdr_cutoff is None})
        if plink and not native and not author:
            required.update({"--fdr-column": fdr_column is None, "--decoy-prefix": decoy_prefix is None})
        missing = [option for option, absent in required.items() if absent]
        if missing:
            raise click.UsageError("Input export requires " + ", ".join(missing) + ".")
        if native and fdr_level != "psm":
            raise click.UsageError("Native pLink Q-value evidence supports --fdr-level psm only.")

    click.echo("Detected format: " + ("pLink PSM" if plink else "xiFDR PSM"))
    try:
        if plink:
            bundle = read_plink_psms(
                results_path, fasta_path, dataset_accession=dataset_accession,
                decoy_prefix=decoy_prefix,
                fdr_columns={fdr_level: fdr_column} if fdr_column is not None else None,
                unfiltered_path=plink_unfiltered_path, parameters_path=plink_parameters_path,
                author_filter_path=author_filter_path,
            )
        else:
            bundle = read_xifdr_psms(results_path, fasta_path, dataset_accession=dataset_accession,
                                     backend="pyxlms")
        if audit_only:
            destination = write_import_report(bundle, output_dir)
            manifest = json.loads((destination / "manifest.json").read_text(encoding="utf-8"))
            click.echo(f"Source PSM records: {manifest['counts']['source_records']}")
            click.echo(f"Unresolved source records: {manifest['counts']['unresolved_source_records']}")
            click.echo(f"Records with parser errors: {manifest['counts']['parser_errors']}")
            if "native_evidence_counts" in manifest:
                click.echo(f"Native confidence verified: {manifest['native_evidence_counts'].get('verified', 0)}")
                click.echo(f"Author-reported PSM FDR cutoff: {manifest['native_filter_policy']['cutoff']}")
            if "author_filter_record_counts" in manifest:
                click.echo(f"Author-filtered PSM records: {manifest['author_filter_record_counts'].get('reported_accepted', 0)}")
                click.echo(f"Author-reported PSM FDR cutoff: {manifest['author_filter_policy']['reported_fdr_cutoff']}")
            click.echo(f"Audit report: {(destination / 'import_report.json').resolve()}")
            click.echo(f"Output: {destination.resolve()}")
            return
        mappings = _read_protein_mapping(mapping_path)
        projection = _read_reference_projection(reference_projection_path) if reference_projection_path is not None else None
        destination = write_claudio_input(
            bundle, mappings, output_dir,
            fdr_level=fdr_level, fdr_cutoff=fdr_cutoff,
            allowed_crosslinkers=crosslinkers, allowed_residues=residues,
            source_protein_ids=source_proteins or None,
            reference_projection=projection,
        )
        manifest = json.loads((destination / "manifest.json").read_text(encoding="utf-8"))
    except (ValueError, OSError, RuntimeError) as error:
        raise click.ClickException(str(error)) from error

    click.echo(f"Source PSM records: {manifest['counts']['source_records']}")
    click.echo(f"Standardized input rows: {manifest['counts']['pipeline_input_rows']}")
    failures = sum(record["record_reason"].startswith("pyxlms_") for record in bundle["records"])
    click.echo(f"Records with parser errors: {failures}")
    click.echo(f"Input: {(destination / 'claudio_input.csv').resolve()}")
    click.echo(f"Audit: {(destination / 'psm_audit.json').resolve()}")
    click.echo(f"Output: {destination.resolve()}")


if __name__ == "__main__":
    main()
