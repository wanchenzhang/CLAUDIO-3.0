"""Command-line tools for discovering public XL-MS dataset candidates."""

import json
from pathlib import Path

import click
import requests

from .outputs import write_discovery_outputs
from .proteomecentral import fetch_dataset_pages


@click.group()
def main():
    """Manage the public dataset candidate catalog."""


@main.command()
@click.option(
    "--search-term", "search_terms",
    multiple=True,
    required=True,
    help="Search term. Repeat this option to search multiple terms.",
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


if __name__ == "__main__":
    main()