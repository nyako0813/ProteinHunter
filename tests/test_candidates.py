"""Tests for candidate assembly helpers."""

from __future__ import annotations

from analysis.candidates import (
    build_candidate_records,
    filter_positive_without_negative,
    get_best_hit,
    group_hits_by_query,
    summarize_blast_status,
)
from core.models import BlastHit, ProteinRecord


def make_hit(
    query_id: str,
    subject_id: str,
    bitscore: float,
    evalue: float = 1e-5,
    source: str = "blast",
) -> BlastHit:
    """Create a BLAST hit for candidate tests."""
    return BlastHit(
        query_id=query_id,
        subject_id=subject_id,
        percent_identity=80.0,
        alignment_length=100,
        evalue=evalue,
        bitscore=bitscore,
        source=source,
    )


def test_group_hits_by_query_preserves_group_order() -> None:
    """Hits should be grouped by query_id without reordering within groups."""
    first = make_hit("protein_1", "subject_a", 20.0)
    second = make_hit("protein_2", "subject_b", 30.0)
    third = make_hit("protein_1", "subject_c", 40.0)

    grouped = group_hits_by_query([first, second, third])

    assert grouped == {
        "protein_1": [first, third],
        "protein_2": [second],
    }


def test_get_best_hit_uses_highest_bitscore() -> None:
    """The best hit should be the one with the highest bitscore."""
    low = make_hit("protein_1", "low", 20.0, 1e-20)
    high = make_hit("protein_1", "high", 40.0, 1e-2)

    assert get_best_hit([low, high]) == high


def test_get_best_hit_tie_breaks_by_lower_evalue() -> None:
    """When bitscores tie, the lower e-value should win."""
    weaker = make_hit("protein_1", "weaker", 40.0, 1e-3)
    stronger = make_hit("protein_1", "stronger", 40.0, 1e-20)

    assert get_best_hit([weaker, stronger]) == stronger


def test_get_best_hit_returns_none_for_empty_list() -> None:
    """An empty hit list should return None."""
    assert get_best_hit([]) is None


def test_build_candidate_records_attaches_available_data() -> None:
    """Candidate records should receive descriptions, sequences, and hits."""
    positive = make_hit("protein_1", "positive_subject", 50.0, source="positive")
    negative = make_hit("protein_1", "negative_subject", 45.0, source="negative")
    unused = make_hit("protein_3", "unused_subject", 70.0)

    records = build_candidate_records(
        protein_ids=["protein_1", "protein_2"],
        descriptions={"protein_1": "Known-like candidate"},
        sequences={"protein_1": "MSTN", "protein_2": "AAAA"},
        positive_hits=[positive, unused],
        negative_hits=[negative],
    )

    assert set(records) == {"protein_1", "protein_2"}
    assert records["protein_1"].description == "Known-like candidate"
    assert records["protein_1"].sequence == "MSTN"
    assert records["protein_1"].positive_hits == [positive]
    assert records["protein_1"].negative_hits == [negative]
    assert records["protein_2"].description == ""
    assert records["protein_2"].sequence == "AAAA"
    assert records["protein_2"].positive_hits == []
    assert records["protein_2"].negative_hits == []


def test_filter_positive_without_negative_returns_positive_only_records() -> None:
    """Only records with positive hits and no negative hits should be returned."""
    positive_only = ProteinRecord(
        protein_id="positive_only",
        positive_hits=[make_hit("positive_only", "subject", 10.0)],
    )
    mixed = ProteinRecord(
        protein_id="mixed",
        positive_hits=[make_hit("mixed", "subject", 10.0)],
        negative_hits=[make_hit("mixed", "negative", 9.0)],
    )
    negative_only = ProteinRecord(
        protein_id="negative_only",
        negative_hits=[make_hit("negative_only", "negative", 9.0)],
    )
    no_hits = ProteinRecord(protein_id="no_hits")
    records = {
        "positive_only": positive_only,
        "mixed": mixed,
        "negative_only": negative_only,
        "no_hits": no_hits,
    }

    filtered = filter_positive_without_negative(records)

    assert filtered == {"positive_only": positive_only}
    assert set(records) == {"positive_only", "mixed", "negative_only", "no_hits"}


def test_summarize_blast_status_all_statuses() -> None:
    """BLAST status labels should cover the four expected hit combinations."""
    positive_hit = make_hit("protein", "positive", 10.0)
    negative_hit = make_hit("protein", "negative", 9.0)

    assert (
        summarize_blast_status(
            ProteinRecord(protein_id="protein", positive_hits=[positive_hit])
        )
        == "positive_only"
    )
    assert (
        summarize_blast_status(
            ProteinRecord(
                protein_id="protein",
                positive_hits=[positive_hit],
                negative_hits=[negative_hit],
            )
        )
        == "positive_and_negative"
    )
    assert (
        summarize_blast_status(
            ProteinRecord(protein_id="protein", negative_hits=[negative_hit])
        )
        == "negative_only"
    )
    assert summarize_blast_status(ProteinRecord(protein_id="protein")) == "no_hits"
