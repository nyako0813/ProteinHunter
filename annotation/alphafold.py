"""Small AlphaFold DB helpers for ProteinHunter annotations."""

from __future__ import annotations

import requests

from core.cache import JsonCache
from core.exceptions import AlphaFoldAnnotationError


ALPHAFOLD_ENTRY_URL = "https://alphafold.ebi.ac.uk/entry/{accession}"


def build_alphafold_url(accession: str | None) -> str | None:
    """Return the AlphaFold DB entry URL for a UniProt accession."""
    if accession is None or not accession.strip():
        return None

    return ALPHAFOLD_ENTRY_URL.format(accession=accession.strip())


def check_alphafold_exists(
    accession: str,
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> bool:
    """Return True when AlphaFold DB has a prediction for an accession."""
    cache_key = accession.strip()
    if cache is not None and cache.has("alphafold", cache_key):
        return bool(cache.get("alphafold", cache_key))

    url = build_alphafold_url(cache_key)
    if url is None:
        return False

    try:
        response = requests.head(url, allow_redirects=True, timeout=timeout)
    except requests.RequestException as exc:
        raise AlphaFoldAnnotationError(
            f"AlphaFold DB check failed for '{accession}'. Please check the network connection."
        ) from exc

    if 200 <= response.status_code < 400:
        exists = True
    elif response.status_code == 404:
        exists = False
    else:
        raise AlphaFoldAnnotationError(
            f"AlphaFold DB returned status {response.status_code} for '{accession}'."
        )

    if cache is not None:
        cache.set("alphafold", cache_key, exists)

    return exists


def get_alphafold_url_if_exists(
    accession: str | None,
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> str | None:
    """Return an AlphaFold DB URL only when the entry exists."""
    url = build_alphafold_url(accession)

    if url is None or accession is None:
        return None

    if check_alphafold_exists(accession, cache=cache, timeout=timeout):
        return url

    return None


__all__: tuple[str, ...] = (
    "ALPHAFOLD_ENTRY_URL",
    "build_alphafold_url",
    "check_alphafold_exists",
    "get_alphafold_url_if_exists",
)
