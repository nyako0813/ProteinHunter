"""Candidate record assembly helpers based on BLAST results."""

from __future__ import annotations

from core.models import BlastHit, ProteinRecord


def group_hits_by_query(hits: list[BlastHit]) -> dict[str, list[BlastHit]]:
    """Group BLAST hits by query protein while preserving input order."""
    grouped_hits: dict[str, list[BlastHit]] = {}

    for hit in hits:
        grouped_hits.setdefault(hit.query_id, []).append(hit)

    return grouped_hits


def get_best_hit(hits: list[BlastHit]) -> BlastHit | None:
    """Return the best hit by highest bitscore, then lowest e-value."""
    if not hits:
        return None

    return max(hits, key=lambda hit: (hit.bitscore, -hit.evalue))


def build_candidate_records(
    protein_ids: list[str],
    descriptions: dict[str, str] | None = None,
    sequences: dict[str, str] | None = None,
    positive_hits: list[BlastHit] | None = None,
    negative_hits: list[BlastHit] | None = None,
) -> dict[str, ProteinRecord]:
    """Build ProteinRecord objects and attach grouped BLAST hits."""
    descriptions = descriptions or {}
    sequences = sequences or {}
    grouped_positive = group_hits_by_query(positive_hits or [])
    grouped_negative = group_hits_by_query(negative_hits or [])

    records: dict[str, ProteinRecord] = {}

    for protein_id in protein_ids:
        records[protein_id] = ProteinRecord(
            protein_id=protein_id,
            description=descriptions.get(protein_id, ""),
            sequence=sequences.get(protein_id, ""),
            positive_hits=list(grouped_positive.get(protein_id, [])),
            negative_hits=list(grouped_negative.get(protein_id, [])),
        )

    return records


def filter_positive_without_negative(
    records: dict[str, ProteinRecord],
) -> dict[str, ProteinRecord]:
    """Return records with positive BLAST hits and no negative BLAST hits."""
    return {
        protein_id: record
        for protein_id, record in records.items()
        if record.positive_hits and not record.negative_hits
    }


def summarize_blast_status(record: ProteinRecord) -> str:
    """Return a simple status label for positive and negative BLAST hits."""
    has_positive = bool(record.positive_hits)
    has_negative = bool(record.negative_hits)

    if has_positive and has_negative:
        return "positive_and_negative"

    if has_positive:
        return "positive_only"

    if has_negative:
        return "negative_only"

    return "no_hits"


__all__: tuple[str, ...] = (
    "build_candidate_records",
    "filter_positive_without_negative",
    "get_best_hit",
    "group_hits_by_query",
    "summarize_blast_status",
)
