"""Tests for ProteinHunter data models."""

from __future__ import annotations

from core.models import (
    AnnotationResult,
    BlastHit,
    CandidateScore,
    DomainHit,
    ProteinRecord,
)


def make_blast_hit() -> BlastHit:
    """Create a small BLAST hit for tests."""
    return BlastHit(
        query_id="query_1",
        subject_id="subject_1",
        percent_identity=91.5,
        alignment_length=120,
        evalue=1e-20,
        bitscore=88.0,
    )


def test_annotation_result_defaults_are_independent() -> None:
    """AnnotationResult list and dict defaults should not be shared."""
    first = AnnotationResult(protein_id="protein_1", source="pfam", success=True)
    second = AnnotationResult(protein_id="protein_2", source="pfam", success=True)

    first.domains.append(DomainHit(source="pfam", accession="PF00001", name="Domain"))
    first.motifs.append("CXXC")
    first.metadata["hits"] = 1

    assert second.domains == []
    assert second.motifs == []
    assert second.metadata == {}


def test_candidate_score_defaults_are_independent() -> None:
    """CandidateScore list and dict defaults should not be shared."""
    first = CandidateScore(protein_id="protein_1")
    second = CandidateScore(protein_id="protein_2")

    first.components["blast"] = 3.0
    first.reasons.append("Strong BLAST support")

    assert second.components == {}
    assert second.reasons == []


def test_protein_record_defaults_are_independent() -> None:
    """ProteinRecord list and dict defaults should not be shared."""
    first = ProteinRecord(protein_id="protein_1")
    second = ProteinRecord(protein_id="protein_2")

    first.positive_hits.append(make_blast_hit())
    first.negative_hits.append(make_blast_hit())
    first.domains.append(DomainHit(source="cdd", accession="CDD:1", name="Domain"))
    first.motifs.append("HXXH")
    first.annotations["pfam"] = AnnotationResult(
        protein_id="protein_1",
        source="pfam",
        success=True,
    )
    first.notes.append("Reviewed")

    assert second.positive_hits == []
    assert second.negative_hits == []
    assert second.domains == []
    assert second.motifs == []
    assert second.annotations == {}
    assert second.notes == []


def test_candidate_score_add_component_updates_total_score() -> None:
    """Adding score components should keep total_score in sync."""
    score = CandidateScore(protein_id="protein_1")

    score.add_component("blast", 2.5, "Positive BLAST hit")
    score.add_component("domain", 4.0)
    score.add_component("blast", 3.0, "Updated BLAST evidence")

    assert score.components == {"blast": 3.0, "domain": 4.0}
    assert score.total_score == 7.0
    assert score.reasons == ["Positive BLAST hit", "Updated BLAST evidence"]


def test_protein_record_length_and_has_negative_hit() -> None:
    """ProteinRecord properties should reflect sequence and negative hits."""
    record = ProteinRecord(protein_id="protein_1", sequence="MSTNPKPQR")

    assert record.length == 9
    assert record.has_negative_hit is False

    record.negative_hits.append(make_blast_hit())

    assert record.has_negative_hit is True
