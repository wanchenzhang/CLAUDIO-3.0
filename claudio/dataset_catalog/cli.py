"""Command-line tools for discovering public XL-MS dataset candidates."""

import json
import math
from collections import Counter
from pathlib import Path

import click
import requests

from .inventory import inspect_catalog, summarize_inventories, inspect_content_plan
from .proteomecentral import fetch_dataset_pages, write_discovery_outputs

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


if __name__ == "__main__":
    main()