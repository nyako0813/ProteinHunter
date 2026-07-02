"""Small UniProt REST API helpers for ProteinHunter annotations."""

from __future__ import annotations

import re
from typing import Any

import requests

from core.cache import JsonCache
from core.exceptions import UniProtAnnotationError


Metadata = dict[str, str | int | float | bool | None]
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
OLD_LOCUS_TAG_PATTERN = re.compile(r"\bMA_\d{4}\b")


def search_uniprot_by_protein_id(
    protein_id: str,
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> Metadata:
    """Search UniProt for a protein ID and return compact metadata."""
    if cache is not None and cache.has("uniprot", protein_id):
        cached = cache.get("uniprot", protein_id)
        if isinstance(cached, dict):
            return _coerce_metadata(cached)

    params = {
        "query": protein_id,
        "format": "json",
        "size": "1",
        "fields": "accession,id,protein_name,organism_name,reviewed,gene_names",
    }

    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
    except requests.RequestException as exc:
        raise UniProtAnnotationError(
            f"UniProt search failed for '{protein_id}'. Please check the network connection."
        ) from exc
    except ValueError as exc:
        raise UniProtAnnotationError(
            f"UniProt returned invalid JSON for '{protein_id}'."
        ) from exc

    try:
        metadata = _parse_uniprot_payload(protein_id, payload)
    except (KeyError, TypeError) as exc:
        raise UniProtAnnotationError(
            f"UniProt returned an unexpected response for '{protein_id}'."
        ) from exc

    if cache is not None:
        cache.set("uniprot", protein_id, metadata)

    return metadata


def extract_uniprot_accession(metadata: dict[str, object]) -> str | None:
    """Return a UniProt accession string when present and usable."""
    accession = metadata.get("accession")

    if isinstance(accession, str) and accession.strip():
        return accession

    return None


def extract_uniprot_old_locus_tag(metadata: dict[str, object]) -> str | None:
    """Return a UniProt old/locus tag string when present and usable."""
    old_locus_tag = metadata.get("old_locus_tag")

    if isinstance(old_locus_tag, str) and old_locus_tag.strip():
        return old_locus_tag

    return _find_old_locus_tag_in_text(metadata)


def _parse_uniprot_payload(protein_id: str, payload: Any) -> Metadata:
    """Convert a UniProt search response into compact metadata."""
    if not isinstance(payload, dict):
        raise TypeError("UniProt payload must be a JSON object.")

    results = payload.get("results", [])
    if not isinstance(results, list):
        raise TypeError("UniProt results must be a list.")

    if not results:
        return {
            "query": protein_id,
            "accession": None,
            "id": None,
            "protein_name": None,
            "organism": None,
            "reviewed": False,
            "old_locus_tag": None,
            "old_locus_tag_note": _old_locus_tag_note(None),
        }

    first_result = results[0]
    if not isinstance(first_result, dict):
        raise TypeError("UniProt result must be a JSON object.")

    old_locus_tag = _old_locus_tag(first_result)

    return {
        "query": protein_id,
        "accession": _optional_string(first_result.get("primaryAccession")),
        "id": _optional_string(first_result.get("uniProtkbId")),
        "protein_name": _protein_name(first_result),
        "organism": _organism_name(first_result),
        "reviewed": _is_reviewed(first_result),
        "old_locus_tag": old_locus_tag,
        "old_locus_tag_note": _old_locus_tag_note(old_locus_tag),
    }


def _coerce_metadata(value: dict[str, Any]) -> Metadata:
    """Return cached metadata with the expected value type."""
    old_locus_tag = _optional_string(value.get("old_locus_tag"))
    metadata: Metadata = {
        "query": _optional_string(value.get("query")),
        "accession": _optional_string(value.get("accession")),
        "id": _optional_string(value.get("id")),
        "protein_name": _optional_string(value.get("protein_name")),
        "organism": _optional_string(value.get("organism")),
        "reviewed": bool(value.get("reviewed")),
        "old_locus_tag": old_locus_tag,
    }

    inferred_tag = extract_uniprot_old_locus_tag(metadata)
    metadata["old_locus_tag"] = inferred_tag
    metadata["old_locus_tag_note"] = _old_locus_tag_note(inferred_tag)
    return metadata


def _optional_string(value: object) -> str | None:
    """Return a non-empty string, or None."""
    if isinstance(value, str) and value.strip():
        return value

    return None


def _protein_name(result: dict[str, Any]) -> str | None:
    """Extract the recommended protein name when available."""
    description = result.get("proteinDescription")
    if not isinstance(description, dict):
        return None

    recommended = description.get("recommendedName")
    if not isinstance(recommended, dict):
        return None

    full_name = recommended.get("fullName")
    if not isinstance(full_name, dict):
        return None

    return _optional_string(full_name.get("value"))


def _organism_name(result: dict[str, Any]) -> str | None:
    """Extract the organism scientific name when available."""
    organism = result.get("organism")
    if not isinstance(organism, dict):
        return None

    return _optional_string(organism.get("scientificName"))


def _is_reviewed(result: dict[str, Any]) -> bool:
    """Return True when UniProt marks the entry as reviewed."""
    entry_type = result.get("entryType")

    if isinstance(entry_type, str):
        return "reviewed" in entry_type.lower()

    reviewed = result.get("reviewed")
    return bool(reviewed)


def _old_locus_tag(result: dict[str, Any]) -> str | None:
    """Extract the first ordered locus name, falling back to the first ORF name."""
    genes = result.get("genes")
    if isinstance(genes, list):
        for gene in genes:
            if not isinstance(gene, dict):
                continue

            tag = _first_gene_name(gene.get("orderedLocusNames"))
            if tag is not None:
                return tag

        for gene in genes:
            if not isinstance(gene, dict):
                continue

            tag = _first_gene_name(gene.get("orfNames"))
            if tag is not None:
                return tag

    return _find_old_locus_tag_in_text(result)


def _first_gene_name(value: object) -> str | None:
    """Return the first gene name from UniProt gene-name structures."""
    if not isinstance(value, list):
        return None

    for item in value:
        if isinstance(item, str):
            name = _old_locus_tag_from_string(item)
            if name is not None:
                return name

            continue

        if not isinstance(item, dict):
            continue

        name = _optional_string(item.get("value"))
        if name is not None:
            tag = _old_locus_tag_from_string(name)
            if tag is not None:
                return tag

    return None


def _find_old_locus_tag_in_text(value: object) -> str | None:
    """Search nested UniProt metadata text for an MA_#### locus tag."""
    if isinstance(value, str):
        return _old_locus_tag_from_string(value)

    if isinstance(value, dict):
        for item in value.values():
            tag = _find_old_locus_tag_in_text(item)
            if tag is not None:
                return tag

    if isinstance(value, list):
        for item in value:
            tag = _find_old_locus_tag_in_text(item)
            if tag is not None:
                return tag

    return None


def _old_locus_tag_from_string(value: str) -> str | None:
    """Return the first MA_#### tag in one text value."""
    match = OLD_LOCUS_TAG_PATTERN.search(value)

    if match:
        return match.group(0)

    return None


def _old_locus_tag_note(old_locus_tag: str | None) -> str | None:
    """Return a diagnostic note when UniProt did not provide an MA tag."""
    if old_locus_tag is not None:
        return None

    return "UniProt metadata did not contain an MA_#### old locus tag."


__all__: tuple[str, ...] = (
    "Metadata",
    "OLD_LOCUS_TAG_PATTERN",
    "UNIPROT_SEARCH_URL",
    "extract_uniprot_accession",
    "extract_uniprot_old_locus_tag",
    "search_uniprot_by_protein_id",
)
