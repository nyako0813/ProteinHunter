"""Tests for candidate scoring helpers."""

from __future__ import annotations

from analysis.scoring import get_sorted_records, score_record, score_records
from core.models import BlastHit, DomainHit, ProteinRecord


def make_hit(
    query_id: str = "protein_1",
    subject_id: str = "subject_1",
) -> BlastHit:
    """Create a small BLAST hit for scoring tests."""
    return BlastHit(
        query_id=query_id,
        subject_id=subject_id,
        percent_identity=80.0,
        alignment_length=100,
        evalue=1e-20,
        bitscore=50.0,
    )


def make_domain() -> DomainHit:
    """Create a small domain hit for scoring tests."""
    return DomainHit(source="CDD", accession="cd12345", name="Domain")


def test_positive_hit_adds_positive_hit_score() -> None:
    """A positive BLAST hit should add the positive_hit component."""
    record = ProteinRecord(
        protein_id="protein_1",
        positive_hits=[make_hit()],
        negative_hits=[make_hit(subject_id="negative")],
    )

    score_record(record)

    assert record.score is not None
    assert record.score.components["positive_hit"] == 5.0
    assert record.score.total_score == 5.0


def test_no_negative_hit_adds_no_negative_hit_score() -> None:
    """No negative BLAST hits should add the no_negative_hit component."""
    record = ProteinRecord(protein_id="protein_1")

    score_record(record)

    assert record.score is not None
    assert record.score.components["no_negative_hit"] == 5.0
    assert record.score.total_score == 5.0


def test_domain_hit_adds_domain_hit_score() -> None:
    """A domain annotation should add the domain_hit component."""
    record = ProteinRecord(
        protein_id="protein_1",
        negative_hits=[make_hit(subject_id="negative")],
        domains=[make_domain()],
    )

    score_record(record)

    assert record.score is not None
    assert record.score.components["domain_hit"] == 4.0
    assert record.score.total_score == 4.0


def test_uniprot_and_alphafold_add_scores() -> None:
    """UniProt and AlphaFold annotations should add their components."""
    record = ProteinRecord(
        protein_id="protein_1",
        negative_hits=[make_hit(subject_id="negative")],
        uniprot_accession="P12345",
        alphafold_url="https://alphafold.ebi.ac.uk/entry/P12345",
    )

    score_record(record)

    assert record.score is not None
    assert record.score.components["uniprot_accession"] == 2.0
    assert record.score.components["alphafold_url"] == 2.0
    assert record.score.total_score == 4.0


def test_notes_add_annotation_warning_penalty() -> None:
    """Annotation notes should add the annotation warning penalty."""
    record = ProteinRecord(
        protein_id="protein_1",
        negative_hits=[make_hit(subject_id="negative")],
        notes=["CDD annotation failed"],
    )

    score_record(record)

    assert record.score is not None
    assert record.score.components["annotation_warning"] == -1.0
    assert record.score.total_score == -1.0


def test_custom_weights_override_defaults() -> None:
    """Custom weights should replace default component scores."""
    record = ProteinRecord(
        protein_id="protein_1",
        positive_hits=[make_hit()],
        domains=[make_domain()],
    )

    score_record(
        record,
        weights={
            "positive_hit": 10.0,
            "no_negative_hit": 3.0,
            "domain_hit": 8.0,
        },
    )

    assert record.score is not None
    assert record.score.components["positive_hit"] == 10.0
    assert record.score.components["no_negative_hit"] == 3.0
    assert record.score.components["domain_hit"] == 8.0
    assert record.score.total_score == 21.0


def test_score_records_scores_multiple_records() -> None:
    """score_records should score each record and return the same dict."""
    records = {
        "protein_1": ProteinRecord("protein_1", positive_hits=[make_hit()]),
        "protein_2": ProteinRecord(
            "protein_2",
            negative_hits=[make_hit(subject_id="negative")],
        ),
    }

    result = score_records(records)

    assert result is records
    assert records["protein_1"].score is not None
    assert records["protein_2"].score is not None
    assert records["protein_1"].score.total_score == 10.0
    assert records["protein_2"].score.total_score == 0.0


def test_get_sorted_records_sorts_by_total_score() -> None:
    """Records should sort by total_score, with unscored records treated as zero."""
    high = ProteinRecord("high", positive_hits=[make_hit("high")])
    middle = ProteinRecord("middle", domains=[make_domain()])
    low = ProteinRecord("low")
    unscored = ProteinRecord("unscored")
    score_record(high)
    score_record(middle)
    score_record(low)
    records = {
        "low": low,
        "unscored": unscored,
        "high": high,
        "middle": middle,
    }

    sorted_records = get_sorted_records(records)

    assert [record.protein_id for record in sorted_records] == [
        "high",
        "middle",
        "low",
        "unscored",
    ]

    ascending_records = get_sorted_records(records, descending=False)
    assert [record.protein_id for record in ascending_records] == [
        "unscored",
        "low",
        "middle",
        "high",
    ]
