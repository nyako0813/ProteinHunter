"""GFF helpers for old/locus tag annotation."""

from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path
from urllib.parse import unquote

from core.models import ProteinRecord


OLD_LOCUS_TAG_PATTERN = re.compile(r"\bMA_\d{4}\b")
COMPACT_OLD_LOCUS_TAG_PATTERN = re.compile(r"\bMA(\d{4})\b")


@dataclass(frozen=True)
class GffFeatureLocation:
    """Genomic coordinates for one GFF feature."""

    contig: str
    start: int
    end: int
    strand: str | None


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

        decoded_value = unquote(value.strip())
        values = [item.strip() for item in decoded_value.split(",") if item.strip()]
        attributes.setdefault(key.strip(), []).extend(values)

    return attributes


def normalize_protein_id(value: str) -> str:
    """Normalize common GFF protein identifier forms for matching."""
    normalized = value.strip()

    normalized = normalized.split()[0] if normalized else ""

    for prefix in ("GenBank:", "Genbank:", "RefSeq:", "NCBI_GP:", "gb:", "ref:"):
        if normalized.lower().startswith(prefix.lower()):
            normalized = normalized[len(prefix) :]

    for prefix in ("cds-", "gene-"):
        if normalized.startswith(prefix):
            normalized = normalized[len(prefix) :]

    return normalized


def load_gff_locus_map(path: str | Path) -> dict[str, str]:
    """Load a protein_id to old/locus tag mapping from a GFF file."""
    gff_path = Path(path).expanduser().resolve()
    mapping: dict[str, str] = {}
    gene_old_locus_tags: dict[str, str] = {}
    cds_attributes: list[dict[str, list[str]]] = []

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
            if feature_type == "gene":
                gene_id = _first_value(attributes, "ID")
                old_locus_tag = _old_locus_tag_from_attributes(attributes)
                if gene_id is not None and old_locus_tag is not None:
                    gene_old_locus_tags[gene_id] = old_locus_tag
                continue

            cds_attributes.append(attributes)

    for attributes in cds_attributes:
        locus_tag = _cds_locus_tag_from_attributes(attributes, gene_old_locus_tags)
        if locus_tag is None:
            continue

        for protein_id in _protein_ids_from_attributes(attributes):
            for key in _protein_id_lookup_keys(protein_id):
                mapping.setdefault(key, locus_tag)

    return mapping


def load_gff_feature_map(path: str | Path) -> dict[str, GffFeatureLocation]:
    """Load protein/locus identifiers to genomic coordinates from a GFF file."""
    gff_path = Path(path).expanduser().resolve()
    feature_map: dict[str, GffFeatureLocation] = {}
    gene_features: dict[str, GffFeatureLocation] = {}
    gene_old_locus_tags: dict[str, str] = {}
    cds_rows: list[tuple[GffFeatureLocation, dict[str, list[str]]]] = []

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

            location = GffFeatureLocation(
                contig=columns[0],
                start=int(columns[3]),
                end=int(columns[4]),
                strand=columns[6] if columns[6] not in {"", "."} else None,
            )
            attributes = parse_gff_attributes(columns[8])

            if feature_type == "gene":
                gene_id = _first_value(attributes, "ID")
                old_locus_tag = _old_locus_tag_from_attributes(attributes)
                if gene_id is not None:
                    gene_features[gene_id] = location
                    feature_map.setdefault(gene_id, location)
                if gene_id is not None and old_locus_tag is not None:
                    gene_old_locus_tags[gene_id] = old_locus_tag
                    feature_map.setdefault(old_locus_tag, location)
                locus_tag = _locus_tag_from_attributes(attributes)
                if locus_tag is not None:
                    feature_map.setdefault(locus_tag, location)
                continue

            cds_rows.append((location, attributes))

    for location, attributes in cds_rows:
        locus_tag = _cds_locus_tag_from_attributes(attributes, gene_old_locus_tags)
        if locus_tag is not None:
            feature_map.setdefault(locus_tag, location)

        for parent_id in attributes.get("Parent", []):
            parent_location = gene_features.get(parent_id)
            if parent_location is not None and locus_tag is not None:
                feature_map.setdefault(locus_tag, parent_location)

        for protein_id in _protein_ids_from_attributes(attributes):
            for key in _protein_id_lookup_keys(protein_id):
                feature_map.setdefault(key, location)

    return feature_map


def annotate_records_with_gff_locus_tags(
    records: dict[str, ProteinRecord],
    mapping: dict[str, str],
) -> int:
    """Apply GFF old/locus tags to records and return updated record count."""
    updated = 0

    for record in records.values():
        protein_id = normalize_protein_id(record.protein_id)
        locus_tag = mapping.get(protein_id)
        if locus_tag is None:
            locus_tag = mapping.get(_without_version(protein_id))
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


def _cds_locus_tag_from_attributes(
    attributes: dict[str, list[str]],
    gene_old_locus_tags: dict[str, str],
) -> str | None:
    """Return the best old/locus tag for a CDS feature."""
    direct_old_locus_tag = _old_locus_tag_from_attributes(attributes)
    if direct_old_locus_tag is not None:
        return direct_old_locus_tag

    for parent_id in attributes.get("Parent", []):
        parent_old_locus_tag = gene_old_locus_tags.get(parent_id)
        if parent_old_locus_tag is not None:
            return parent_old_locus_tag

    return _locus_tag_from_attributes(attributes)


def _old_locus_tag_from_attributes(attributes: dict[str, list[str]]) -> str | None:
    """Return the preferred old_locus_tag from GFF attributes."""
    return _first_ma_tag(attributes.get("old_locus_tag", []))


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


def _protein_id_lookup_keys(protein_id: str) -> list[str]:
    """Return versioned and unversioned lookup keys for a protein ID."""
    keys = [protein_id]
    unversioned = _without_version(protein_id)
    if unversioned != protein_id:
        keys.append(unversioned)

    return keys


def _without_version(protein_id: str) -> str:
    """Return a protein ID without a trailing numeric version."""
    return re.sub(r"\.\d+$", "", protein_id)


def _first_value(attributes: dict[str, list[str]], key: str) -> str | None:
    """Return the first attribute value for a key."""
    values = attributes.get(key, [])
    if not values:
        return None

    return values[0]


def _first_ma_tag(values: list[str]) -> str | None:
    """Return the first MA_#### tag found in a list of strings."""
    for value in values:
        match = OLD_LOCUS_TAG_PATTERN.search(value)
        if match:
            return match.group(0)

    for value in values:
        match = COMPACT_OLD_LOCUS_TAG_PATTERN.search(value)
        if match:
            return f"MA_{match.group(1)}"

    return None


__all__: tuple[str, ...] = (
    "GffFeatureLocation",
    "OLD_LOCUS_TAG_PATTERN",
    "annotate_records_with_gff_locus_tags",
    "load_gff_feature_map",
    "load_gff_locus_map",
    "normalize_protein_id",
    "parse_gff_attributes",
)
