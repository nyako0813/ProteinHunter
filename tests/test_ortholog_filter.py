"""Tests for ortholog-aware negative hit classification."""

from __future__ import annotations

from config import ORTHOLOG_FILTER_DEFAULT
from core.models import BlastHit, ProteinRecord
from analysis.ortholog_filter import (
    classify_negative_hit_strength,
    is_excluded_by_negative_mode,
    populate_negative_hit_evidence,
)


def hit(
    identity: float,
    coverage: float,
    evalue: float,
    query_id: str = "T",
) -> BlastHit:
    """Create a negative hit with a chosen query coverage."""
    return BlastHit(
        query_id=query_id,
        subject_id="negative_subject",
        percent_identity=identity,
        alignment_length=int(coverage),
        evalue=evalue,
        bitscore=50.0,
        source="negative_source",
        query_length=100,
    )


def test_hit_strength_classification_thresholds() -> None:
    """Negative hits should be classified by identity, query coverage, and e-value."""
    assert classify_negative_hit_strength(
        hit(45.0, 80.0, 1e-20),
        ORTHOLOG_FILTER_DEFAULT,
    ) == "strong"
    assert classify_negative_hit_strength(
        hit(35.0, 80.0, 1e-20),
        ORTHOLOG_FILTER_DEFAULT,
    ) == "medium"
    assert classify_negative_hit_strength(
        hit(27.0, 60.0, 1e-4),
        ORTHOLOG_FILTER_DEFAULT,
    ) == "weak"
    assert classify_negative_hit_strength(
        hit(24.0, 80.0, 1e-20),
        ORTHOLOG_FILTER_DEFAULT,
    ) == "none"
    assert classify_negative_hit_strength(
        hit(45.0, 30.0, 1e-20),
        ORTHOLOG_FILTER_DEFAULT,
    ) != "strong"


def test_candidates_relaxed_exclusion_in_strong_only_mode() -> None:
    """Strong-only mode should retain medium/weak/no-negative positive hits."""
    records = {
        "strong": ProteinRecord(
            protein_id="strong",
            sequence="M" * 100,
            positive_hits=[hit(90.0, 100.0, 1e-30, query_id="strong")],
            negative_hits=[hit(45.0, 80.0, 1e-20, query_id="strong")],
        ),
        "medium": ProteinRecord(
            protein_id="medium",
            sequence="M" * 100,
            positive_hits=[hit(90.0, 100.0, 1e-30, query_id="medium")],
            negative_hits=[hit(35.0, 80.0, 1e-20, query_id="medium")],
        ),
        "weak": ProteinRecord(
            protein_id="weak",
            sequence="M" * 100,
            positive_hits=[hit(90.0, 100.0, 1e-30, query_id="weak")],
            negative_hits=[hit(27.0, 60.0, 1e-4, query_id="weak")],
        ),
        "none": ProteinRecord(
            protein_id="none",
            sequence="M" * 100,
            positive_hits=[hit(90.0, 100.0, 1e-30, query_id="none")],
        ),
    }
    populate_negative_hit_evidence(records, ORTHOLOG_FILTER_DEFAULT)

    retained = {
        protein_id
        for protein_id, record in records.items()
        if record.positive_hits
        and not is_excluded_by_negative_mode(record, "strong_only")
    }

    assert retained == {"medium", "weak", "none"}
    assert records["strong"].negative_hit_strength == "strong"
    assert records["medium"].negative_hit_strength == "medium"
    assert records["weak"].negative_hit_strength == "weak"
    assert records["none"].negative_hit_strength == "none"
