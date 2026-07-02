"""GFF helpers for old/locus tag annotation."""

from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import unquote

from core.models import ProteinRecord


OLD_LOCUS_TAG_PATTERN = re.compile(r"\bMA_\d{4}\b")


def parse_gff_attributes(attribute_text: str) -> dict[str, list[str]]:
    """Parse a GFF attribute column into key/value lists."""
    attributes: dict[str, list[str]] = {}

    for part in attribute_text.split(";"):
        if not part.strip():
            continue

        if "=" in part:
            key, value = part.split("=", 1)
        else:
            key, value = part, ""

        values = [unquote(item.strip()) for item in value.split(",") if item.strip()]
        attributes.setdefault(key.strip(), []).extend(values)

    return attributes


def normalize_protein_id(value: str) -> str:
    """Normalize common GFF protein identifier forms for matching."""
    normalized = value.strip()

    for prefix in ("Genbank:", "RefSeq:", "gb:", "ref:"):
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix) :]

    for prefix in ("cds-", "gene-"):
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix) :]

    return normalized


def load_gff_locus_map(path: str | Path) -> dict[str, str]:
    """Load a protein_id to old/locus tag mapping from a GFF file."""
    gff_path = Path(path).expanduser().resolve()
    mapping: dict[str, str] = {}

    with gff_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue

            feature_type = columns[2].lower()
            if feature_type not in {"cds", "gene"}:
                continue

            attributes = parse_gff_attributes(columns[8])
            locus_tag = _locus_tag_from_attributes(attributes)
            if locus_tag is None:
                continue

            for protein_id in _protein_ids_from_attributes(attributes):
                mapping.setdefault(protein_id, locus_tag)

    return mapping


def annotate_records_with_gff_locus_tags(
    records: dict[str, ProteinRecord],
    mapping: dict[str, str],
) -> int:
    """Apply GFF old/locus tags to records and return updated record count."""
    updated = 0

    for record in records.values():
        locus_tag = mapping.get(normalize_protein_id(record.protein_id))
        if locus_tag is None:
            continue

        record.old_locus_tag = locus_tag
        updated += 1

    return updated


def _locus_tag_from_attributes(attributes: dict[str, list[str]]) -> str | None:
    """Return the preferred MA_#### locus tag from GFF attributes."""
    for key in ("old_locus_tag", "locus_tag"):
        tag = _first_ma_tag(attributes.get(key, []))
        if tag is not None:
            return tag

    tag = _first_ma_tag(attributes.get("gene", []))
    if tag is not None:
        return tag

    for values in attributes.values():
        tag = _first_ma_tag(values)
        if tag is not None:
            return tag

    return None


def _protein_ids_from_attributes(attributes: dict[str, list[str]]) -> list[str]:
    """Return normalized protein IDs found in GFF attributes."""
    candidates: list[str] = []

    for key in ("protein_id", "Name", "ID"):
        candidates.extend(attributes.get(key, []))

    for value in attributes.get("Dbxref", []):
        candidates.append(value)

    normalized_ids: list[str] = []
    seen: set[str] = set()

    for candidate in candidates:
        normalized = normalize_protein_id(candidate)
        if not normalized or normalized in seen:
            continue

        seen.add(normalized)
        normalized_ids.append(normalized)

    return normalized_ids


def _first_ma_tag(values: list[str]) -> str | None:
    """Return the first MA_#### tag found in a list of strings."""
    for value in values:
        match = OLD_LOCUS_TAG_PATTERN.search(value)
        if match:
            return match.group(0)

    return None


__all__: tuple[str, ...] = (
    "OLD_LOCUS_TAG_PATTERN",
    "annotate_records_with_gff_locus_tags",
    "load_gff_locus_map",
    "normalize_protein_id",
    "parse_gff_attributes",
)
