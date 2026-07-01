"""Pfam annotation helpers for ProteinHunter."""

from __future__ import annotations

import csv
import re
from io import StringIO

import requests

from core.cache import JsonCache
from core.exceptions import PfamAnnotationError
from core.models import DomainHit


PFAM_SEARCH_URL = "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"


def parse_pfam_response(text: str) -> list[DomainHit]:
    """Parse a tolerant Pfam or HMMER-like response into domain hits."""
    hits: list[DomainHit] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if not line or line.startswith("#") or line.startswith(">"):
            continue

        hit = _parse_pfam_line(line)
        if hit is not None:
            hits.append(hit)

    return hits


def search_pfam_by_sequence(
    protein_id: str,
    sequence: str,
    cache: JsonCache | None = None,
    timeout: int = 60,
) -> list[DomainHit]:
    """Search Pfam/HMMER for one protein sequence and return domain hits."""
    if cache is not None and cache.has("pfam", protein_id):
        cached = cache.get("pfam", protein_id)
        if isinstance(cached, list):
            return [
                domain_hit_from_dict(item)
                for item in cached
                if isinstance(item, dict)
            ]

    if not sequence.strip():
        return []

    data = {
        "seq": sequence.strip(),
        "seqdb": "pfam",
        "output": "text",
    }

    try:
        response = requests.post(PFAM_SEARCH_URL, data=data, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise PfamAnnotationError(
            f"Pfam search failed for '{protein_id}'. Please check the network connection."
        ) from exc

    try:
        hits = parse_pfam_response(response.text)
    except Exception as exc:
        raise PfamAnnotationError(
            f"Pfam returned an unexpected response for '{protein_id}'."
        ) from exc

    if cache is not None:
        cache.set("pfam", protein_id, [domain_hit_to_dict(hit) for hit in hits])

    return hits


def domain_hit_to_dict(hit: DomainHit) -> dict[str, str | int | float | None]:
    """Convert a DomainHit to a JSON-serializable dictionary."""
    return {
        "source": hit.source,
        "accession": hit.accession,
        "name": hit.name,
        "description": hit.description,
        "evalue": hit.evalue,
        "bitscore": hit.bitscore,
        "start": hit.start,
        "end": hit.end,
    }


def domain_hit_from_dict(data: dict[str, object]) -> DomainHit:
    """Convert cached domain data back into a DomainHit."""
    return DomainHit(
        source=_string_or_default(data.get("source"), "Pfam"),
        accession=_string_or_default(data.get("accession"), ""),
        name=_string_or_default(data.get("name"), ""),
        description=_string_or_default(data.get("description"), ""),
        evalue=_optional_float(data.get("evalue")),
        bitscore=_optional_float(data.get("bitscore")),
        start=_optional_int(data.get("start")),
        end=_optional_int(data.get("end")),
    )


def _parse_pfam_line(line: str) -> DomainHit | None:
    """Parse one non-comment Pfam/HMMER-like line."""
    parts = _split_line(line)
    if len(parts) < 2:
        return None

    accession_index = _find_accession_index(parts)
    if accession_index is None:
        return None

    accession = _clean_accession(parts[accession_index])
    name = _find_name(parts, accession_index, accession)
    evalue = _find_labeled_float(parts, ("evalue", "e-value", "eval", "i-evalue"))
    bitscore = _find_labeled_float(parts, ("bitscore", "bit_score", "score"))
    start = _find_labeled_int(parts, ("start", "from", "ali_from", "env_from"))
    end = _find_labeled_int(parts, ("end", "to", "ali_to", "env_to"))

    if evalue is None:
        evalue = _first_float_after(parts, accession_index)

    if bitscore is None:
        bitscore = _second_float_after(parts, accession_index)

    if start is None or end is None:
        coordinates = _find_coordinate_pair(parts)
        if coordinates is not None:
            start, end = coordinates

    description = _description_from_parts(parts, accession_index)

    return DomainHit(
        source="Pfam",
        accession=accession,
        name=name,
        description=description,
        evalue=evalue,
        bitscore=bitscore,
        start=start,
        end=end,
    )


def _split_line(line: str) -> list[str]:
    """Split Pfam-like text using tabs, commas, or whitespace."""
    if "\t" in line:
        return [part.strip() for part in line.split("\t") if part.strip()]

    try:
        csv_parts = next(csv.reader(StringIO(line)))
    except csv.Error:
        csv_parts = []

    if len(csv_parts) > 1:
        return [part.strip() for part in csv_parts if part.strip()]

    return [part.strip() for part in re.split(r"\s{2,}|\s+", line) if part.strip()]


def _find_accession_index(parts: list[str]) -> int | None:
    """Find the first Pfam-like accession field."""
    for index, part in enumerate(parts):
        if re.search(r"\bPF\d{5}(?:\.\d+)?\b", part, re.IGNORECASE):
            return index

    return None


def _clean_accession(value: str) -> str:
    """Return the Pfam accession without surrounding punctuation."""
    match = re.search(r"\bPF\d{5}(?:\.\d+)?\b", value, re.IGNORECASE)
    if match:
        return match.group(0)

    return value.strip()


def _find_name(parts: list[str], accession_index: int, default: str) -> str:
    """Find a readable Pfam name near the accession."""
    if accession_index > 0 and not _looks_like_value(parts[accession_index - 1]):
        return parts[accession_index - 1]

    if accession_index + 1 < len(parts):
        return parts[accession_index + 1]

    return default


def _description_from_parts(parts: list[str], accession_index: int) -> str:
    """Return a readable description when extra fields are available."""
    start_index = accession_index + 2
    description_parts = [
        part
        for part in parts[start_index:]
        if not _looks_like_labeled_value(part)
        and _optional_float(part) is None
        and _optional_int(part) is None
        and _find_coordinate_pair([part]) is None
    ]
    return " ".join(description_parts)


def _find_labeled_float(parts: list[str], labels: tuple[str, ...]) -> float | None:
    """Find a float from fields like evalue=1e-5 or score:42.0."""
    for part in parts:
        key, value = _split_labeled_value(part)
        if key in labels:
            return _optional_float(value)

    return None


def _find_labeled_int(parts: list[str], labels: tuple[str, ...]) -> int | None:
    """Find an int from fields like start=10 or end:80."""
    for part in parts:
        key, value = _split_labeled_value(part)
        if key in labels:
            return _optional_int(value)

    return None


def _split_labeled_value(part: str) -> tuple[str, str]:
    """Split a key-value field into lowercase key and raw value."""
    if "=" in part:
        key, value = part.split("=", 1)
    elif ":" in part:
        key, value = part.split(":", 1)
    else:
        return "", part

    return key.strip().lower(), value.strip()


def _looks_like_labeled_value(part: str) -> bool:
    """Return True for simple key=value or key:value fields."""
    return "=" in part or ":" in part


def _looks_like_value(part: str) -> bool:
    """Return True for fields that look numeric, labeled, or coordinate-like."""
    return (
        _looks_like_labeled_value(part)
        or _optional_float(part) is not None
        or _find_coordinate_pair([part]) is not None
    )


def _first_float_after(parts: list[str], index: int) -> float | None:
    """Return the first float after an index."""
    for part in parts[index + 1 :]:
        value = _optional_float(part)
        if value is not None:
            return value

    return None


def _second_float_after(parts: list[str], index: int) -> float | None:
    """Return the second float after an index."""
    seen = 0

    for part in parts[index + 1 :]:
        value = _optional_float(part)
        if value is None:
            continue

        seen += 1
        if seen == 2:
            return value

    return None


def _find_coordinate_pair(parts: list[str]) -> tuple[int, int] | None:
    """Find a simple coordinate field such as 10-80."""
    for part in parts:
        match = re.fullmatch(r"(\d+)\s*[-.]+\s*(\d+)", part)
        if match:
            return int(match.group(1)), int(match.group(2))

    return None


def _string_or_default(value: object, default: str) -> str:
    """Return a string value or a fallback."""
    if isinstance(value, str):
        return value

    return default


def _optional_float(value: object) -> float | None:
    """Return a float when the value can be parsed."""
    if value is None or isinstance(value, bool):
        return None

    try:
        return float(str(value).strip())
    except ValueError:
        return None


def _optional_int(value: object) -> int | None:
    """Return an int when the value can be parsed."""
    if value is None or isinstance(value, bool):
        return None

    try:
        return int(str(value).strip())
    except ValueError:
        return None


__all__: tuple[str, ...] = (
    "PFAM_SEARCH_URL",
    "domain_hit_from_dict",
    "domain_hit_to_dict",
    "parse_pfam_response",
    "search_pfam_by_sequence",
)
