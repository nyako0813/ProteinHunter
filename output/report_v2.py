"""Phase 6-8 Stage 1: unified 12-sheet workbook consolidation.

Builds the data rows for the redesigned Excel workbook (see
claude/phase678_excel_word_redesign_investigation.md and the Stage 1
implementation directive it followed) by merging what used to be up to 10
base classification sheets and up to 11 ``Interaction_*`` bucket sheets
into a single ``candidate_source`` column per sheet. This module only
builds rows/DataFrames; ``output/excel.py`` still owns actual worksheet
writing (formatting, hyperlinks, the Index sheet).

Two independent axes are consolidated:

- **Base classification** (``analysis/blast_pipeline.py::BlastClassificationResult``):
  every candidate source bucket (Candidates, Candidates_relaxed, ...,
  Negative_hit) is always computed, regardless of whether
  interaction_scoring is enabled. Used for 03_Candidate_Overview, and as
  the fallback source for 02_Final_Score when interaction_scoring produced
  no rows for a candidate (or is disabled entirely).
- **Interaction scoring** (``analysis/interaction_scoring.py::InteractionScoringResult.source_rows``):
  one row per (query, candidate, candidate_source bucket) -- when more than
  one enabled bucket contains the same candidate for the same query (e.g.
  Candidates is a subset of Candidates_relaxed, or negative_hit overlaps
  its own strong/medium/weak sub-buckets, see PR #11), this module keeps
  exactly one row per (query, candidate) pair.

Both axes use the same bucket-priority order for deduplication/consolidation,
derived from ``analysis.interaction_scoring.CANDIDATE_PRIORITY_BASE`` (the
existing implicit "how positive is this bucket" ordering interaction_scoring
already uses inside its own source_classification scoring component).
"""

from __future__ import annotations

from typing import Any

from analysis.interaction_scoring import (
    CANDIDATE_PRIORITY_BASE,
    compute_protein_hunter_only_final_score,
)
from analysis.scoring_engine_config import ScoringEngineConfig
from core.models import ProteinRecord

#: Negative_strong/medium/weak_hit are strict subsets of Negative_hit (PR
#: #11's duplicate-scoring finding) -- for candidate_source consolidation
#: they are all the same bucket; their strength lives in the
#: negative_hit_strength column instead (see analysis/interaction_scoring.py
#: M1: both _score_pair and _score_pair_v2 now expose it directly).
_NEGATIVE_SUBBUCKET_NORMALIZE: dict[str, str] = {
    "Negative_strong_hit": "Negative_hit",
    "Negative_medium_hit": "Negative_hit",
    "Negative_weak_hit": "Negative_hit",
}

#: Highest-to-lowest dedup priority for the 6 buckets a consolidated
#: candidate_source can resolve to. Matches CANDIDATE_PRIORITY_BASE's own
#: score ordering (30/30/25/20/15/0); Candidates is preferred over
#: Positive_all_sources on their 30/30 tie only because it is listed first
#: here -- an arbitrary but deterministic tiebreak (both are already
#: strict-positive, negative-free buckets, so which one "wins" the display
#: label does not change any score).
CANDIDATE_SOURCE_DEDUP_ORDER: tuple[str, ...] = (
    "Candidates",
    "Positive_all_sources",
    "Candidates_relaxed",
    "No_hit",
    "Negative_unmatched",
    "Negative_hit",
)

#: blast_classification attribute holding each base bucket's records, in the
#: same priority order as CANDIDATE_SOURCE_DEDUP_ORDER.
_BASE_BUCKET_ATTRS: tuple[str, ...] = (
    "positive_only_records",
    "positive_all_sources_records",
    "candidates_relaxed_records",
    "no_hit_records",
    "negative_unmatched_records",
    "negative_hit_records",
)


def normalize_candidate_source(raw_source: str) -> str:
    """Collapse the negative_hit strength sub-buckets down to "Negative_hit"."""
    return _NEGATIVE_SUBBUCKET_NORMALIZE.get(raw_source, raw_source)


def _dedup_priority(raw_source: str) -> tuple[float, int]:
    """Sort key for picking the single winning bucket for one candidate.

    Higher is better: primary key is CANDIDATE_PRIORITY_BASE's own score
    (already the existing "how positive is this bucket" ordering); the
    secondary key breaks the Candidates/Positive_all_sources 30/30 tie
    using CANDIDATE_SOURCE_DEDUP_ORDER's fixed position (earlier wins).
    """
    normalized = normalize_candidate_source(raw_source)
    base = CANDIDATE_PRIORITY_BASE.get(normalized, 10.0)
    try:
        tie_rank = CANDIDATE_SOURCE_DEDUP_ORDER.index(normalized)
    except ValueError:
        tie_rank = len(CANDIDATE_SOURCE_DEDUP_ORDER)
    return (base, -tie_rank)


def candidate_source_for_protein(protein_id: str, blast_classification: Any) -> str:
    """Return the single consolidated candidate_source label for one protein.

    Iterates the base classification buckets in dedup-priority order and
    returns the first (highest-priority) bucket the protein belongs to.
    Every protein in ``blast_classification.all_records`` falls into at
    least one bucket (No_hit is the catch-all for "no positive, no
    negative"), so "Unclassified" should not normally occur -- it only
    guards against a caller passing a protein_id outside ``all_records``.
    """
    for attr, label in zip(_BASE_BUCKET_ATTRS, CANDIDATE_SOURCE_DEDUP_ORDER):
        bucket = getattr(blast_classification, attr, None) or {}
        if protein_id in bucket:
            return label
    return "Unclassified"


def consolidate_interaction_rows(
    source_rows: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    """Merge every Interaction_* bucket's rows into one row per (query, candidate).

    When the same candidate was scored under more than one enabled
    candidate_sources bucket for the same query, keeps the row from the
    highest dedup-priority bucket and normalizes its displayed
    candidate_source (see normalize_candidate_source). Row order is not
    meaningful here -- callers that need a ranked order should sort/rerank
    the result themselves (see rerank_final_score_rows).
    """
    best: dict[tuple[str, str], dict[str, Any]] = {}
    best_priority: dict[tuple[str, str], tuple[float, int]] = {}

    for rows in source_rows.values():
        for row in rows:
            key = (str(row["query_id"]), str(row["candidate_protein_id"]))
            priority = _dedup_priority(str(row["candidate_source"]))
            if key not in best or priority > best_priority[key]:
                merged = dict(row)
                merged["candidate_source"] = normalize_candidate_source(str(row["candidate_source"]))
                best[key] = merged
                best_priority[key] = priority

    return list(best.values())


def rerank_final_score_rows(rows: list[dict[str, Any]]) -> None:
    """Sort ``rows`` in place into 02_Final_Score's display order and reassign candidate_rank.

    Reuses analysis.interaction_scoring._rerank_by_final_score for the
    actual rank assignment (identical tie-break rules: eligible rows
    1..N descending by final_score, AlphaFold-recommended first, then
    candidate_protein_id; ineligible rows get candidate_rank=0) -- this
    sheet is explicitly "final-score-first" by design regardless of
    whichever ranking_metric the run itself used for the original
    Interaction_* sheets' own candidate_rank.
    """
    from analysis.interaction_scoring import _rerank_by_final_score  # local: intentionally private, tightly coupled

    _rerank_by_final_score(rows)
    rows.sort(
        key=lambda row: (
            str(row["query_id"]),
            row.get("candidate_rank") in (0, None),
            row.get("candidate_rank") or 0,
            str(row["candidate_protein_id"]),
        )
    )


def build_base_overview_rows(blast_classification: Any) -> list[dict[str, Any]]:
    """Return one 03_Candidate_Overview-shaped row per protein in all_records.

    Global identity + consolidated candidate_source + negative_hit_strength,
    independent of any query -- see candidate_source_for_protein.
    """
    rows: list[dict[str, Any]] = []
    for protein_id, record in blast_classification.all_records.items():
        source = candidate_source_for_protein(protein_id, blast_classification)
        rows.append(_overview_row(record, source))
    return rows


def _overview_row(record: ProteinRecord, candidate_source: str) -> dict[str, Any]:
    from output.excel import _best_hit  # local: avoids a hard import-time cycle with output/excel.py

    best_positive = _best_hit(record.positive_hits)
    return {
        "protein_id": record.protein_id,
        "old_locus_tag": record.old_locus_tag or "",
        "description": record.description,
        "candidate_source": candidate_source,
        "negative_hit_strength": record.negative_hit_strength or "none",
        "protein_hunter_score": record.score.total_score if record.score else None,
        "sequence_length": record.length,
        "blast_status": _blast_status(record),
        "positive_hit_count": len(record.positive_hits),
        "negative_hit_count": len(record.negative_hits),
        "best_positive_hit": best_positive.subject_id if best_positive else None,
        "domain_count": len(record.domains),
        "uniprot_accession": record.uniprot_accession,
        "alphafold_url": record.alphafold_url,
        "positive_source_count": record.positive_source_count,
        "positive_sources_hit": "; ".join(record.positive_sources_hit),
        "positive_sources_missing": "; ".join(record.positive_sources_missing),
        "notes": "; ".join(record.notes),
    }


def _blast_status(record: ProteinRecord) -> str:
    has_positive = bool(record.positive_hits)
    has_negative = bool(record.negative_hits)
    if has_positive and has_negative:
        return "positive_and_negative"
    if has_positive:
        return "positive_only"
    if has_negative:
        return "negative_only"
    return "no_hits"


def apply_wider_protein_hunter_scores(
    overview_rows: list[dict[str, Any]],
    protein_hunter_scores: dict[str, Any],
) -> None:
    """Overwrite protein_hunter_score in place with interaction_scoring's wider-scope values.

    ``record.score`` (used by _overview_row's fallback) is only ever
    computed for positive_only_records ("Candidates") -- see
    resolve_protein_hunter_scores's own docstring. When interaction_scoring
    is enabled, its resolve_protein_hunter_scores(...) result covers every
    candidate_sources bucket instead, so 03_Candidate_Overview should
    prefer that wider value whenever it exists for a given protein.
    """
    for row in overview_rows:
        score = protein_hunter_scores.get(row["protein_id"])
        if score is not None:
            row["protein_hunter_score"] = score.total_score


def build_no_query_final_score_rows(
    overview_rows: list[dict[str, Any]],
    engine_config: ScoringEngineConfig,
) -> list[dict[str, Any]]:
    """02_Final_Score fallback rows when interaction_scoring produced no data at all.

    Every candidate gets Final Score computed from protein_hunter_score
    alone (compute_protein_hunter_only_final_score), the same fallback
    every in-run candidate with a MISSING interaction_score already gets
    from _attach_final_score_columns -- just exposed here for the
    "no query configured / interaction_scoring disabled" case, which never
    reaches that per-query code path at all.
    """
    rows: list[dict[str, Any]] = []
    for overview_row in overview_rows:
        final_score, tier = compute_protein_hunter_only_final_score(
            overview_row.get("protein_hunter_score"), engine_config
        )
        rows.append(
            {
                "query_id": "",
                "query_protein_id": "",
                "candidate_rank": 0,
                "candidate_protein_id": overview_row["protein_id"],
                "candidate_old_locus_tag": overview_row["old_locus_tag"],
                "candidate_source": overview_row["candidate_source"],
                "candidate_description": overview_row["description"],
                "negative_hit_strength": overview_row["negative_hit_strength"],
                "protein_hunter_score": overview_row.get("protein_hunter_score"),
                "interaction_score": None,
                "final_score": final_score,
                "final_score_tier": tier,
                "candidate_priority_score": None,
                "functional_domain_score": None,
                "evolutionary_score": None,
                "same_gene_neighborhood_score": None,
                "interaction_evidence_score": None,
                "evidence_category_count": None,
                "evidence_component_count": None,
                "available_weight_total": None,
                "alphafold_recommended": None,
            }
        )
    rerank_final_score_rows(rows)
    return rows


__all__: tuple[str, ...] = (
    "CANDIDATE_SOURCE_DEDUP_ORDER",
    "apply_wider_protein_hunter_scores",
    "build_base_overview_rows",
    "build_no_query_final_score_rows",
    "candidate_source_for_protein",
    "consolidate_interaction_rows",
    "normalize_candidate_source",
    "rerank_final_score_rows",
)
