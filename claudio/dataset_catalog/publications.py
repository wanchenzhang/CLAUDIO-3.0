"""Extract explicit paper identifiers without verifying publication status."""

import re
from html import unescape
from html.parser import HTMLParser
from urllib.parse import unquote, urlsplit


_DOI_HOSTS = {"doi.org", "dx.doi.org", "www.doi.org"}
_PUBMED_HOSTS = {"pubmed.ncbi.nlm.nih.gov"}
_LEGACY_PUBMED_HOSTS = {"www.ncbi.nlm.nih.gov", "ncbi.nlm.nih.gov"}


def normalize_doi(value):
    """Return a canonical DOI from an identifier or trusted resolver URL."""
    if not isinstance(value, str):
        return None
    text = unescape(value).strip()
    if text[:4].lower() == "doi:":
        text = text[4:].strip()
    if not text:
        return None
    try:
        uri = urlsplit(text)
        if uri.scheme:
            if (
                uri.scheme.lower() not in {"http", "https"}
                or uri.netloc.lower() not in _DOI_HOSTS
            ):
                return None
            text = unquote(
                uri.path.lstrip("/"), encoding="utf-8", errors="strict"
            )
    except (ValueError, UnicodeError):
        return None
    if re.fullmatch(r"10\.[0-9]{4,9}/[^\s<>\"'?#]+", text) is None:
        return None
    return text.lower()


def normalize_pmid(value):
    """Return a positive PMID from a number or trusted PubMed URL."""
    if type(value) is int:
        return str(value) if value > 0 else None
    if not isinstance(value, str):
        return None
    text = unescape(value).strip()
    if text[:5].lower() == "pmid:":
        text = text[5:].strip()
    try:
        uri = urlsplit(text)
        if uri.scheme:
            if uri.scheme.lower() not in {"http", "https"}:
                return None
            if uri.netloc.lower() in _PUBMED_HOSTS:
                match = re.fullmatch(r"/([0-9]+)/?", uri.path)
            elif uri.netloc.lower() in _LEGACY_PUBMED_HOSTS:
                match = re.fullmatch(r"/pubmed/([0-9]+)/?", uri.path)
            else:
                return None
            if match is None:
                return None
            text = match.group(1)
        if re.fullmatch(r"[0-9]+", text) is None:
            return None
        number = int(text)
    except ValueError:
        return None
    return str(number) if number > 0 else None


def _canonical_accession(value):
    if not isinstance(value, str):
        raise ValueError("Dataset accession must be a nonempty string")
    accession = value.strip().upper()
    if re.fullmatch(r"[A-Z][A-Z0-9._-]*", accession) is None:
        raise ValueError("Dataset accession is missing or invalid")
    return accession


def _repository_accession(metadata):
    value = metadata.get("accession")
    if isinstance(value, str):
        return _canonical_accession(value)
    if isinstance(value, dict):
        value = [value]
    if isinstance(value, list):
        bindings = [
            entry.get("value")
            for entry in value
            if isinstance(entry, dict)
            and entry.get("accession") == "MS:1001919"
        ]
        if len(bindings) == 1:
            return _canonical_accession(bindings[0])
    raise ValueError("Repository metadata has no unique dataset accession")


class _PublicationHTML(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.text = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() == "a":
            for key, value in attrs:
                if key.lower() == "href" and isinstance(value, str):
                    self.links.append(value)

    def handle_data(self, data):
        self.text.append(data)


def extract_publication_links(candidate, *, repository_metadata=None):
    """Extract associated identifiers while retaining unverified status."""
    if not isinstance(candidate, dict):
        raise ValueError("Candidate must contain a dataset accession")
    accession = _canonical_accession(candidate.get("accession"))
    if repository_metadata is not None:
        if not isinstance(repository_metadata, dict):
            raise ValueError("Repository dataset metadata must be an object")
        if _repository_accession(repository_metadata) != accession:
            raise ValueError("Repository accession differs from the dataset")

    result = {
        "dataset_accession": accession,
        "publication_status": "unverified",
        "references": [],
        "indicators": [],
    }
    reference_index = {}
    indicator_index = {}

    def add_reference(doi, pmid, paths):
        doi = normalize_doi(doi)
        pmid = normalize_pmid(pmid)
        if doi is None and pmid is None:
            return
        key = (doi, pmid)
        if key not in reference_index:
            reference = {
                "doi": doi, "pmid": pmid, "source_fields": [],
            }
            reference_index[key] = reference
            result["references"].append(reference)
        for path in paths:
            if path not in reference_index[key]["source_fields"]:
                reference_index[key]["source_fields"].append(path)

    def add_indicator(kind, path, raw):
        key = (kind, raw)
        if key not in indicator_index:
            indicator = {
                "kind": kind, "raw_value": raw, "source_fields": [],
            }
            indicator_index[key] = indicator
            result["indicators"].append(indicator)
        if path not in indicator_index[key]["source_fields"]:
            indicator_index[key]["source_fields"].append(path)

    def scan_text(value, path):
        parser = _PublicationHTML()
        parser.feed(value)
        parser.close()
        for link in parser.links:
            add_reference(link, link, [path])

        if "<" not in value and ">" not in value:
            stripped = value.strip()
            pmid = (
                None if re.fullmatch(r"[0-9]+", stripped)
                else normalize_pmid(stripped)
            )
            add_reference(normalize_doi(stripped), pmid, [path])

        visible = " ".join(" ".join(parser.text).lower().split())
        if re.search(
            r"\b(?:publication|manuscript)\b.{0,40}\bpending\b", visible
        ):
            add_indicator("pending", path, value)
        if re.search(
            r"\bno\s+(?:(?:associated|published)\s+)*"
            r"(?:publication|manuscript|paper)\b",
            visible,
        ):
            add_indicator("no_publication", path, value)

    def scan_value(value, path):
        if isinstance(value, str):
            scan_text(value, path)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                scan_value(item, f"{path}[{index}]")
        elif isinstance(value, dict):
            term = value.get("accession")
            if term == "MS:1002853":
                add_indicator(
                    "no_publication", path,
                    str(value.get("name", "MS:1002853")),
                )
            elif term == "MS:1000879":
                add_reference(None, value.get("value"), [path])
            elif term == "MS:1001922":
                add_reference(value.get("value"), None, [path])
            else:
                doi = normalize_doi(value.get("doi"))
                pmid_key = "pubmedID" if "pubmedID" in value else "pmid"
                pmid = normalize_pmid(value.get(pmid_key))
                paths = []
                if doi is not None:
                    paths.append(path + ".doi")
                if pmid is not None:
                    paths.append(path + "." + pmid_key)
                add_reference(doi, pmid, paths)
                for key in ("cvParams", "publications", "references"):
                    if key in value:
                        scan_value(value[key], path + "." + key)

    def scan_fields(metadata, prefix):
        for key in ("publications", "references"):
            if key in metadata:
                scan_value(metadata[key], prefix + "." + key)

    scan_fields(candidate, "candidate")
    variants = candidate.get("metadata_variants", [])
    if isinstance(variants, list):
        for index, variant in enumerate(variants):
            if (
                isinstance(variant, dict)
                and isinstance(variant.get("metadata"), dict)
            ):
                scan_fields(
                    variant["metadata"],
                    f"candidate.metadata_variants[{index}].metadata",
                )
    if repository_metadata is not None:
        scan_fields(repository_metadata, "repository")
    return result
