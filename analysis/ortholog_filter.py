"""Ortholog-aware negative BLAST hit classification helpers."""

from __future__ import annotations

from typing import Any

from core.models import BlastHit, ProteinRecord


NEGATIVE_STRENGTH_ORDER: dict[str, int] = {
    "none": 0,
    "weak": 1,
    "medium": 2,
    "strong": 3,
}


def classify_negative_hit_strength(
    hit: BlastHit,
    ortholog_filter: Any,
    query_length: int | None = None,
) -> str:
    """Classify one negative BLAST hit as strong, medium, weak, or none."""
    coverage = hit.query_coverage
    if coverage is None and query_length is not None and query_length > 0:
        coverage = hit.alignment_length / query_length * 100
    if coverage is None:
        return "none"

    if _passes_threshold(hit, coverage, ortholog_filter.strong):
        return "strong"
    if _passes_threshold(hit, coverage, ortholog_filter.medium):
        return "medium"
    if _passes_threshold(hit, coverage, ortholog_filter.weak):
        return "weak"

    return "none"


def populate_negative_hit_evidence(
    records: dict[str, ProteinRecord],
    ortholog_filter: Any,
) -> None:
    """Populate negative hit strength and exclusion evidence on each record."""
    for record in records.values():
        _populate_record_negative_evidence(record, ortholog_filter)


def is_excluded_by_negative_mode(record: ProteinRecord, mode: str) -> bool:
    """Return whether a record is excluded under the configured negative mode."""
    if mode == "none":
        return False
    if mode == "any_hit":
        return bool(record.negative_hits)
    if mode == "strong_only":
        return record.negative_strong_hit_count > 0
    if mode == "strong_or_medium":
        return (
            record.negative_strong_hit_count > 0
            or record.negative_medium_hit_count > 0
        )

    return bool(record.negative_hits)


def negative_exclusion_reason(record: ProteinRecord, mode: str) -> str:
    """Return a compact explanation for retention or exclusion."""
    if not record.negative_hits:
        return "no negative hit"
    if mode == "none":
        return "retained: negative exclusion mode is none"
    if mode == "any_hit":
        return "excluded: negative hit present; exclusion mode is any_hit"
    if record.negative_strong_hit_count > 0:
        if mode in {"strong_only", "strong_or_medium"}:
            return "excluded: strong negative hit"
        return "retained: strong negative hit; exclusion mode is none"
    if record.negative_medium_hit_count > 0:
        if mode == "strong_or_medium":
            return "excluded: medium negative hit; exclusion mode is strong_or_medium"
        return "retained: medium negative hit; exclusion mode is strong_only"
    if record.negative_weak_hit_count > 0:
        return "retained: weak negative hit only"

    return "retained: negative hit below weak threshold"


def _populate_record_negative_evidence(
    record: ProteinRecord,
    ortholog_filter: Any,
) -> None:
    counts = {"strong": 0, "medium": 0, "weak": 0}
    representative_hit: BlastHit | None = None
    representative_strength = "none"
    query_length = record.length or None

    for hit in record.negative_hits:
        if hit.query_length is None and query_length is not None:
            hit.query_length = query_length
        strength = classify_negative_hit_strength(
            hit,
            ortholog_filter,
            query_length=query_length,
        )
        if strength in counts:
            counts[strength] += 1
        if _is_better_representative(
            hit,
            strength,
            representative_hit,
            representative_strength,
        ):
            representative_hit = hit
            representative_strength = strength

    record.negative_strong_hit_count = counts["strong"]
    record.negative_medium_hit_count = counts["medium"]
    record.negative_weak_hit_count = counts["weak"]
    record.negative_hit_strength = representative_strength
    if representative_hit is None:
        record.negative_best_identity = None
        record.negative_best_query_coverage = None
        record.negative_best_evalue = None
        record.negative_best_source = None
    else:
        record.negative_best_identity = representative_hit.percent_identity
        record.negative_best_query_coverage = representative_hit.query_coverage
        record.negative_best_evalue = representative_hit.evalue
        record.negative_best_source = representative_hit.source

    record.negative_exclusion_reason = negative_exclusion_reason(
        record,
        ortholog_filter.negative_exclusion_mode,
    )


def _passes_threshold(hit: BlastHit, coverage: float, threshold: Any) -> bool:
    return (
        hit.percent_identity >= threshold.min_identity
        and coverage >= threshold.min_query_coverage
        and hit.evalue <= threshold.max_evalue
    )


def _is_better_representative(
    hit: BlastHit,
    strength: str,
    current_hit: BlastHit | None,
    current_strength: str,
) -> bool:
    if current_hit is None:
        return True

    strength_rank = NEGATIVE_STRENGTH_ORDER[strength]
    current_rank = NEGATIVE_STRENGTH_ORDER[current_strength]
    if strength_rank != current_rank:
        return strength_rank > current_rank

    return (hit.bitscore, -hit.evalue) > (current_hit.bitscore, -current_hit.evalue)


__all__: tuple[str, ...] = (
    "NEGATIVE_STRENGTH_ORDER",
    "classify_negative_hit_strength",
    "is_excluded_by_negative_mode",
    "negative_exclusion_reason",
    "populate_negative_hit_evidence",
)
