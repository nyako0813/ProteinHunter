"""Tests for lightweight interaction candidate ranking."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from config import (
    INTERACTION_ALPHAFOLD_DEFAULT,
    INTERACTION_CANDIDATE_SOURCE_DEFAULTS,
    INTERACTION_SCORING_WEIGHTS_DEFAULT,
    InteractionQueryConfig,
    InteractionScoringConfig,
)
from core.models import DomainHit, ProteinRecord
from analysis.interaction_scoring import (
    interaction_pair_columns,
    run_interaction_scoring,
)


def interaction_config(
    *,
    enabled: bool = True,
    query_proteins: tuple[InteractionQueryConfig, ...] = (),
    query_fasta: Path | None = None,
    candidate_sources: dict[str, bool] | None = None,
    max_candidates_per_query: int = 200,
    include_sequences_in_excel: bool = False,
) -> SimpleNamespace:
    """Build a minimal app config for interaction scoring tests."""
    return SimpleNamespace(
        interaction_scoring=InteractionScoringConfig(
            enabled=enabled,
            query_proteins=query_proteins,
            query_fasta=query_fasta,
            candidate_sources=(
                dict(INTERACTION_CANDIDATE_SOURCE_DEFAULTS)
                if candidate_sources is None
                else candidate_sources
            ),
            max_candidates_per_query=max_candidates_per_query,
            include_sequences_in_excel=include_sequences_in_excel,
            scoring_weights=INTERACTION_SCORING_WEIGHTS_DEFAULT,
            alphafold=INTERACTION_ALPHAFOLD_DEFAULT,
        )
    )


def record(
    protein_id: str,
    *,
    old_locus_tag: str = "",
    sequence: str = "MSTNPKPQR",
    description: str = "enzyme protein",
    positive_sources_hit: list[str] | None = None,
) -> ProteinRecord:
    """Create a target/candidate protein record."""
    return ProteinRecord(
        protein_id=protein_id,
        old_locus_tag=old_locus_tag or None,
        sequence=sequence,
        description=description,
        positive_sources_hit=positive_sources_hit or [],
    )


def classification(records: dict[str, ProteinRecord]) -> SimpleNamespace:
    """Build a classification object with candidate source groups."""
    return SimpleNamespace(
        all_records=records,
        positive_only_records={"candidate": records["candidate"]},
        candidates_relaxed_records={
            "candidate": records["candidate"],
            "relaxed": records["relaxed"],
        },
        positive_all_sources_records={},
        negative_unmatched_records={},
        no_hit_records={"novel": records["novel"]},
        negative_hit_records={},
        negative_strong_hit_records={},
        negative_medium_hit_records={},
        negative_weak_hit_records={},
    )


def test_interaction_scoring_disabled_returns_none() -> None:
    """Disabled interaction scoring should not create output."""
    cfg = interaction_config(enabled=False)
    result = run_interaction_scoring(cfg, SimpleNamespace(all_records={}))

    assert result is None


def test_query_resolution_by_protein_id_old_locus_and_sequence() -> None:
    """Queries should resolve by protein_id, old_locus_tag, or explicit sequence."""
    records = {
        "query": record("query", old_locus_tag="MA_0001"),
        "locus_query": record("locus_query", old_locus_tag="MA_0002"),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(
            InteractionQueryConfig("query", "", ""),
            InteractionQueryConfig("", "MA_0002", ""),
            InteractionQueryConfig("external", "", "MPEPTIDE"),
            InteractionQueryConfig("missing", "", ""),
        ),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    statuses = {row["query_id"]: row["resolution_status"] for row in result.query_rows}
    assert statuses["query"] == "resolved"
    assert statuses["MA_0002"] == "resolved"
    assert statuses["external"] == "resolved"
    assert statuses["missing"] == "unresolved"
    external = next(row for row in result.query_rows if row["query_id"] == "external")
    assert external["sequence_length"] == len("MPEPTIDE")


def test_query_fasta_is_loaded(tmp_path: Path) -> None:
    """query_fasta records should become interaction queries."""
    query_fasta = tmp_path / "queries.faa"
    query_fasta.write_text(">query_fasta_1 description\nMSTNPK\n", encoding="utf-8")
    records = {
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_fasta=query_fasta,
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert result.query_rows[0]["query_id"] == "query_fasta_1"
    assert result.query_rows[0]["resolution_status"] == "resolved"


def test_candidate_sources_true_false_and_ranking_limit_are_respected() -> None:
    """Only enabled source sheets should produce interaction rows."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed", description="carrier protein"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={
            "candidates": True,
            "candidates_relaxed": False,
            "no_hit": True,
        },
        max_candidates_per_query=1,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert set(result.source_rows) == {"Interaction_Candidates", "Interaction_No_hit"}
    assert len(result.source_rows["Interaction_Candidates"]) == 1
    assert result.source_rows["Interaction_Candidates"][0]["candidate_rank"] == 1
    assert "Interaction_Candidates_relaxed" not in result.source_rows


def test_missing_enabled_source_is_skipped_with_warning() -> None:
    """Enabled but empty candidate sources should not crash or create a sheet."""
    records = {
        "query": record("query"),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"positive_all_sources": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert result.source_rows == {}
    assert result.warnings == [
        "interaction candidate source has no records: positive_all_sources"
    ]


def test_include_sequences_controls_pair_columns() -> None:
    """Sequence columns should only be present when explicitly requested."""
    assert "query_sequence" not in interaction_pair_columns(False)
    assert "candidate_sequence" not in interaction_pair_columns(False)
    assert "query_sequence" in interaction_pair_columns(True)
    assert "candidate_sequence" in interaction_pair_columns(True)


def test_scoring_adds_reasons_and_alphafold_readiness_without_running_alphafold() -> None:
    """Scoring should use lightweight evidence only."""
    query = record(
        "query",
        description="enzyme protein",
        positive_sources_hit=["A"],
    )
    candidate = record(
        "candidate",
        description="carrier enzyme protein",
        positive_sources_hit=["A"],
    )
    candidate.domains.append(
        DomainHit(source="Pfam", accession="PF00001", name="carrier")
    )
    records = {
        "query": query,
        "candidate": candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_priority_score"] > 0
    assert row["alphafold_recommended"] is True
    assert "compatible for manual AlphaFold" in row["interaction_score_reasons"]
