"""Small UniProt REST API helpers for ProteinHunter annotations."""

from __future__ import annotations

from typing import Any

import requests

from core.cache import JsonCache
from core.exceptions import UniProtAnnotationError


Metadata = dict[str, str | int | float | bool | None]
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


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
        "fields": "accession,id,protein_name,organism_name,reviewed",
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
        }

    first_result = results[0]
    if not isinstance(first_result, dict):
        raise TypeError("UniProt result must be a JSON object.")

    return {
        "query": protein_id,
        "accession": _optional_string(first_result.get("primaryAccession")),
        "id": _optional_string(first_result.get("uniProtkbId")),
        "protein_name": _protein_name(first_result),
        "organism": _organism_name(first_result),
        "reviewed": _is_reviewed(first_result),
    }


def _coerce_metadata(value: dict[str, Any]) -> Metadata:
    """Return cached metadata with the expected value type."""
    return {
        "query": _optional_string(value.get("query")),
        "accession": _optional_string(value.get("accession")),
        "id": _optional_string(value.get("id")),
        "protein_name": _optional_string(value.get("protein_name")),
        "organism": _optional_string(value.get("organism")),
        "reviewed": bool(value.get("reviewed")),
    }


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


__all__: tuple[str, ...] = (
    "Metadata",
    "UNIPROT_SEARCH_URL",
    "extract_uniprot_accession",
    "search_uniprot_by_protein_id",
)
