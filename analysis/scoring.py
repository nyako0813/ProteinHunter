"""Simple candidate scoring helpers for ProteinHunter."""

from __future__ import annotations

from core.models import CandidateScore, ProteinRecord


DEFAULT_WEIGHTS: dict[str, float] = {
    "positive_hit": 5.0,
    "no_negative_hit": 5.0,
    "domain_hit": 4.0,
    "uniprot_accession": 2.0,
    "alphafold_url": 2.0,
    "annotation_warning": -1.0,
}


def score_record(
    record: ProteinRecord,
    weights: dict[str, float] | None = None,
) -> ProteinRecord:
    """Score one protein record in place and return it."""
    active_weights = DEFAULT_WEIGHTS | (weights or {})
    score = CandidateScore(protein_id=record.protein_id)

    if record.positive_hits:
        score.add_component(
            "positive_hit",
            active_weights["positive_hit"],
            "This protein has at least one positive BLAST hit.",
        )

    if not record.negative_hits:
        score.add_component(
            "no_negative_hit",
            active_weights["no_negative_hit"],
            "This protein has no negative BLAST hits.",
        )

    if record.domains:
        score.add_component(
            "domain_hit",
            active_weights["domain_hit"],
            "This protein has at least one domain annotation.",
        )

    if record.uniprot_accession:
        score.add_component(
            "uniprot_accession",
            active_weights["uniprot_accession"],
            "A UniProt accession was found for this protein.",
        )

    if record.alphafold_url:
        score.add_component(
            "alphafold_url",
            active_weights["alphafold_url"],
            "An AlphaFold structure link was found for this protein.",
        )

    if record.notes:
        score.add_component(
            "annotation_warning",
            active_weights["annotation_warning"],
            "This protein has annotation notes that may need review.",
        )

    record.score = score
    return record


def score_records(
    records: dict[str, ProteinRecord],
    weights: dict[str, float] | None = None,
) -> dict[str, ProteinRecord]:
    """Score all records in place and return the same dictionary."""
    for record in records.values():
        score_record(record, weights=weights)

    return records


def get_sorted_records(
    records: dict[str, ProteinRecord],
    descending: bool = True,
) -> list[ProteinRecord]:
    """Return records sorted by total score."""
    return sorted(
        records.values(),
        key=lambda record: record.score.total_score if record.score else 0.0,
        reverse=descending,
    )


__all__: tuple[str, ...] = (
    "DEFAULT_WEIGHTS",
    "get_sorted_records",
    "score_record",
    "score_records",
)
