"""Loader for a UniProtKB REST bulk-export JSON (static file, no live API calls).

This is a static-file counterpart to annotation/cdd.py: instead of querying
NCBI's CDD service per-protein, it reads one UniProtKB search-result export
(``{"results": [...]}``) covering an entire organism/strain, and indexes it
by ``old_locus_tag`` (e.g. ``"MA_4115"``) so ``domain_family_map.py`` can
classify query/candidate proteins into coarse domain-family categories using
UniProt's own Pfam/InterPro/SUPFAM cross-references -- see
claude_code_instructions_domain_complementarity_v3.md section 3.2.

Unlike ``core/cache.py``'s ``JsonCache`` (broken cache = hard error), a
missing or malformed bulk export here means "the new domain-family feature
is disabled" -- callers get an empty mapping rather than an exception, so
existing environments that never set up this optional file keep working
unchanged.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

#: UniProtKB cross-reference databases relevant to domain-family classification.
_RELEVANT_DATABASES: frozenset[str] = frozenset({"Pfam", "InterPro", "SUPFAM"})


@dataclass(frozen=True, slots=True)
class UniProtDomainInfo:
    """One protein's domain-family cross-references from a UniProtKB entry."""

    pfam: frozenset[str] = frozenset()
    interpro: frozenset[str] = frozenset()
    supfam: frozenset[str] = frozenset()
    protein_existence: str = ""


def load_uniprot_bulk_domain_map(path: str | Path) -> dict[str, UniProtDomainInfo]:
    """Load a UniProtKB bulk-export JSON into an ``old_locus_tag -> UniProtDomainInfo`` map.

    ``genes[].orderedLocusNames[].value`` entries are used as old_locus_tag
    keys (an entry with multiple such names is registered under all of
    them). Only ``uniProtKBCrossReferences`` entries whose ``database`` is
    Pfam, InterPro, or SUPFAM are collected, into their respective
    frozensets.

    Returns an empty dict -- never raises -- when ``path`` does not exist,
    is not readable, or is not valid JSON matching the expected
    ``{"results": [...]}`` shape. This is a deliberate "missing bulk data
    disables the feature" policy, not the JsonCache "corrupt cache is a
    hard error" one.
    """
    return _load_uniprot_bulk_domain_map_cached(str(path))


@lru_cache(maxsize=4)
def _load_uniprot_bulk_domain_map_cached(path_str: str) -> dict[str, UniProtDomainInfo]:
    """Process-lifetime memoized parse -- re-parsing a ~36MB export per call is wasteful."""
    path = Path(path_str)
    try:
        raw_text = path.read_text(encoding="utf-8")
        payload = json.loads(raw_text)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError):
        return {}

    if not isinstance(payload, dict):
        return {}
    results = payload.get("results")
    if not isinstance(results, list):
        return {}

    domain_map: dict[str, UniProtDomainInfo] = {}
    for entry in results:
        if not isinstance(entry, dict):
            continue
        old_locus_tags = _old_locus_tags(entry)
        if not old_locus_tags:
            continue

        info = UniProtDomainInfo(
            pfam=frozenset(_cross_reference_ids(entry, "Pfam")),
            interpro=frozenset(_cross_reference_ids(entry, "InterPro")),
            supfam=frozenset(_cross_reference_ids(entry, "SUPFAM")),
            protein_existence=str(entry.get("proteinExistence") or ""),
        )
        for tag in old_locus_tags:
            domain_map[tag] = info

    return domain_map


def _old_locus_tags(entry: dict[str, object]) -> list[str]:
    """Return every ``genes[].orderedLocusNames[].value`` string on one entry."""
    genes = entry.get("genes")
    if not isinstance(genes, list):
        return []

    tags: list[str] = []
    for gene in genes:
        if not isinstance(gene, dict):
            continue
        ordered_locus_names = gene.get("orderedLocusNames")
        if not isinstance(ordered_locus_names, list):
            continue
        for name in ordered_locus_names:
            if not isinstance(name, dict):
                continue
            value = name.get("value")
            if isinstance(value, str) and value:
                tags.append(value)
    return tags


def _cross_reference_ids(entry: dict[str, object], database: str) -> list[str]:
    """Return every ``id`` from ``uniProtKBCrossReferences`` matching ``database``."""
    cross_references = entry.get("uniProtKBCrossReferences")
    if not isinstance(cross_references, list):
        return []

    ids: list[str] = []
    for cross_reference in cross_references:
        if not isinstance(cross_reference, dict):
            continue
        if cross_reference.get("database") != database:
            continue
        cross_reference_id = cross_reference.get("id")
        if isinstance(cross_reference_id, str) and cross_reference_id:
            ids.append(cross_reference_id)
    return ids


__all__: tuple[str, ...] = (
    "UniProtDomainInfo",
    "load_uniprot_bulk_domain_map",
)
