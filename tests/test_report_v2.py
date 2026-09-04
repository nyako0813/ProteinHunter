"""Tests for the Phase 6-8 Stage 1 unified 12-sheet consolidation layer."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from config import (
    INTERACTION_ALPHAFOLD_DEFAULT,
    INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    INTERACTION_NEIGHBORHOOD_DEFAULT,
    INTERACTION_SCORING_WEIGHTS_DEFAULT,
    InteractionScoringConfig,
)
from core.models import ProteinRecord
from analysis.scoring_engine_config import load_scoring_engine_config
from output.report_v2 import (
    TIER_SAFETY_NET,
    apply_wider_protein_hunter_scores,
    bookmark_name,
    build_base_overview_rows,
    build_no_query_final_score_rows,
    build_workbook_sheets,
    candidate_source_for_protein,
    consolidate_interaction_rows,
    normalize_candidate_source,
    rerank_final_score_rows,
    select_top_candidates_per_query,
)


def record(protein_id: str, **kwargs) -> ProteinRecord:
    defaults: dict = dict(protein_id=protein_id, old_locus_tag=None, sequence="MSTNPK", description="protein")
    defaults.update(kwargs)
    return ProteinRecord(**defaults)


def blast_classification(**buckets: dict[str, ProteinRecord]) -> SimpleNamespace:
    defaults: dict[str, dict[str, ProteinRecord]] = {
        "all_records": {},
        "positive_only_records": {},
        "candidates_relaxed_records": {},
        "positive_all_sources_records": {},
        "negative_unmatched_records": {},
        "no_hit_records": {},
        "negative_hit_records": {},
        "negative_strong_hit_records": {},
        "negative_medium_hit_records": {},
        "negative_weak_hit_records": {},
    }
    defaults.update(buckets)
    return SimpleNamespace(**defaults)


# ---------------------------------------------------------------------------
# normalize_candidate_source
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("Negative_strong_hit", "Negative_hit"),
        ("Negative_medium_hit", "Negative_hit"),
        ("Negative_weak_hit", "Negative_hit"),
        ("Negative_hit", "Negative_hit"),
        ("Candidates", "Candidates"),
        ("Positive_all_sources", "Positive_all_sources"),
    ],
)
def test_normalize_candidate_source(raw: str, expected: str) -> None:
    assert normalize_candidate_source(raw) == expected


# ---------------------------------------------------------------------------
# candidate_source_for_protein
# ---------------------------------------------------------------------------


def test_candidate_source_for_protein_prefers_strict_candidates_over_relaxed() -> None:
    """A protein in both Candidates and Candidates_relaxed should resolve to Candidates."""
    p = record("p1")
    bc = blast_classification(
        all_records={"p1": p},
        positive_only_records={"p1": p},
        candidates_relaxed_records={"p1": p},
    )
    assert candidate_source_for_protein("p1", bc) == "Candidates"


def test_candidate_source_for_protein_negative_hit_when_only_in_negative_hit() -> None:
    p = record("p1")
    bc = blast_classification(all_records={"p1": p}, negative_hit_records={"p1": p})
    assert candidate_source_for_protein("p1", bc) == "Negative_hit"


def test_candidate_source_for_protein_unclassified_when_in_no_bucket() -> None:
    p = record("p1")
    bc = blast_classification(all_records={"p1": p})
    assert candidate_source_for_protein("p1", bc) == "Unclassified"


# ---------------------------------------------------------------------------
# consolidate_interaction_rows
# ---------------------------------------------------------------------------


def _pair_row(query_id: str, candidate_id: str, candidate_source: str, **extra) -> dict:
    base = {
        "query_id": query_id,
        "candidate_protein_id": candidate_id,
        "candidate_source": candidate_source,
        "final_score": 10.0,
        "alphafold_recommended": False,
    }
    base.update(extra)
    return base


def test_consolidate_interaction_rows_dedups_across_overlapping_buckets() -> None:
    """A candidate scored under both Candidates and Candidates_relaxed for the same query -> one row, Candidates wins."""
    source_rows = {
        "Interaction_Candidates": [_pair_row("q1", "c1", "Candidates", final_score=50.0)],
        "Interaction_Candidates_relaxed": [_pair_row("q1", "c1", "Candidates_relaxed", final_score=40.0)],
    }
    merged = consolidate_interaction_rows(source_rows)
    assert len(merged) == 1
    assert merged[0]["candidate_source"] == "Candidates"
    assert merged[0]["final_score"] == 50.0


def test_consolidate_interaction_rows_normalizes_negative_subbuckets() -> None:
    """A candidate scored under both negative_hit and negative_strong_hit (PR #11 duplicate) collapses to one Negative_hit row."""
    source_rows = {
        "Interaction_Neg_hit": [_pair_row("q1", "c1", "Negative_hit")],
        "Interaction_Neg_strong": [_pair_row("q1", "c1", "Negative_strong_hit")],
    }
    merged = consolidate_interaction_rows(source_rows)
    assert len(merged) == 1
    assert merged[0]["candidate_source"] == "Negative_hit"


def test_consolidate_interaction_rows_keeps_distinct_pairs_separate() -> None:
    source_rows = {
        "Interaction_Candidates": [
            _pair_row("q1", "c1", "Candidates"),
            _pair_row("q2", "c1", "Candidates"),
            _pair_row("q1", "c2", "Candidates"),
        ],
    }
    merged = consolidate_interaction_rows(source_rows)
    keys = {(r["query_id"], r["candidate_protein_id"]) for r in merged}
    assert keys == {("q1", "c1"), ("q2", "c1"), ("q1", "c2")}


# ---------------------------------------------------------------------------
# rerank_final_score_rows
# ---------------------------------------------------------------------------


def test_rerank_final_score_rows_orders_by_final_score_descending_within_query() -> None:
    rows = [
        _pair_row("q1", "low", "Candidates", final_score=10.0),
        _pair_row("q1", "high", "Candidates", final_score=90.0),
        _pair_row("q2", "only", "Candidates", final_score=5.0),
    ]
    rerank_final_score_rows(rows)
    q1_rows = [r for r in rows if r["query_id"] == "q1"]
    assert q1_rows[0]["candidate_protein_id"] == "high"
    assert q1_rows[0]["candidate_rank"] == 1
    assert q1_rows[1]["candidate_protein_id"] == "low"
    assert q1_rows[1]["candidate_rank"] == 2


def test_rerank_final_score_rows_puts_ineligible_last() -> None:
    rows = [
        _pair_row("q1", "eligible", "Candidates", final_score=10.0),
        _pair_row("q1", "ineligible", "Candidates", final_score=None),
    ]
    rerank_final_score_rows(rows)
    assert rows[0]["candidate_protein_id"] == "eligible"
    assert rows[1]["candidate_protein_id"] == "ineligible"
    assert rows[1]["candidate_rank"] == 0


# ---------------------------------------------------------------------------
# build_base_overview_rows / apply_wider_protein_hunter_scores
# ---------------------------------------------------------------------------


def test_build_base_overview_rows_one_row_per_protein_with_consolidated_source() -> None:
    p1 = record("p1", positive_hits=[])
    p1.negative_hit_strength = "weak"
    bc = blast_classification(
        all_records={"p1": p1},
        negative_hit_records={"p1": p1},
        negative_weak_hit_records={"p1": p1},
    )
    rows = build_base_overview_rows(bc)
    assert len(rows) == 1
    assert rows[0]["protein_id"] == "p1"
    assert rows[0]["candidate_source"] == "Negative_hit"
    assert rows[0]["negative_hit_strength"] == "weak"


def test_apply_wider_protein_hunter_scores_overrides_in_place() -> None:
    rows = [{"protein_id": "p1", "protein_hunter_score": None}]
    apply_wider_protein_hunter_scores(rows, {"p1": SimpleNamespace(total_score=14)})
    assert rows[0]["protein_hunter_score"] == 14


def test_apply_wider_protein_hunter_scores_leaves_unmatched_rows_untouched() -> None:
    rows = [{"protein_id": "p1", "protein_hunter_score": 5}]
    apply_wider_protein_hunter_scores(rows, {})
    assert rows[0]["protein_hunter_score"] == 5


# ---------------------------------------------------------------------------
# build_no_query_final_score_rows
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# build_workbook_sheets
# ---------------------------------------------------------------------------


def test_build_workbook_sheets_falls_back_when_interaction_scoring_disabled() -> None:
    p1 = record("p1")
    p1.score = SimpleNamespace(total_score=9, components={}, reasons=[])
    bc = blast_classification(all_records={"p1": p1}, positive_only_records={"p1": p1})
    config = SimpleNamespace(interaction_scoring=SimpleNamespace(enabled=False, scoring_engine_config=None))

    sheets = build_workbook_sheets(config, bc, None)

    assert len(sheets["overview_rows"]) == 1
    assert len(sheets["final_score_rows"]) == 1
    assert sheets["final_score_rows"][0]["candidate_protein_id"] == "p1"
    assert sheets["evidence_detail_rows"] == []
    assert sheets["query_rows"] == []
    assert sheets["neighborhood_rows"] == []


def test_build_workbook_sheets_uses_consolidated_interaction_rows_when_available() -> None:
    p1 = record("p1")
    bc = blast_classification(all_records={"p1": p1}, positive_only_records={"p1": p1})
    config = SimpleNamespace(
        interaction_scoring=InteractionScoringConfig(
            enabled=True,
            query_proteins=(),
            query_fasta=None,
            candidate_sources={"candidates": True},
            max_candidates_per_query=200,
            include_sequences_in_excel=False,
            scoring_weights=INTERACTION_SCORING_WEIGHTS_DEFAULT,
            alphafold=INTERACTION_ALPHAFOLD_DEFAULT,
            neighborhood=INTERACTION_NEIGHBORHOOD_DEFAULT,
            scoring_model="legacy_additive",
            scoring_engine_config=None,
            functional_complementarity_ruleset=None,
            pih_evidence_bundle=None,
            evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
        )
    )
    interaction_result = SimpleNamespace(
        source_rows={"Interaction_Candidates": [_pair_row("q1", "p1", "Candidates", final_score=42.0)]},
        evidence_detail_rows=[{"category": "genomic_context"}],
        evidence_detail_scoring_model="v2_evidence_based",
        query_rows=[{"query_id": "q1"}],
        neighborhood_rows=[],
    )

    sheets = build_workbook_sheets(config, bc, interaction_result)

    assert len(sheets["final_score_rows"]) == 1
    assert sheets["final_score_rows"][0]["final_score"] == 42.0
    assert sheets["final_score_rows"][0]["candidate_rank"] == 1
    assert sheets["evidence_detail_scoring_model"] == "v2_evidence_based"
    assert sheets["query_rows"] == [{"query_id": "q1"}]


def test_build_no_query_final_score_rows_uses_protein_hunter_alone() -> None:
    overview_rows = [
        {
            "protein_id": "p1",
            "old_locus_tag": "MA_0001",
            "description": "d",
            "candidate_source": "Candidates",
            "negative_hit_strength": "none",
            "protein_hunter_score": 14,
        },
        {
            "protein_id": "p2",
            "old_locus_tag": "MA_0002",
            "description": "d2",
            "candidate_source": "No_hit",
            "negative_hit_strength": "none",
            "protein_hunter_score": None,
        },
    ]
    engine_config = load_scoring_engine_config(None)
    rows = build_no_query_final_score_rows(overview_rows, engine_config)

    assert len(rows) == 2
    p1_row = next(r for r in rows if r["candidate_protein_id"] == "p1")
    assert p1_row["final_score"] is not None
    assert p1_row["interaction_score"] is None
    p2_row = next(r for r in rows if r["candidate_protein_id"] == "p2")
    assert p2_row["final_score"] is None
    assert p2_row["candidate_rank"] == 0
    # p1 (protein_hunter_score=14) should outrank p2 (no score at all).
    assert p1_row["candidate_rank"] == 1


# ---------------------------------------------------------------------------
# Phase 6-8 Stage 2: Word report candidate selection / bookmark naming
# ---------------------------------------------------------------------------


def _ranked_row(query_id: str, rank: int, tier: str = "Tier4_Weak", candidate: str | None = None) -> dict:
    return {
        "query_id": query_id,
        "candidate_protein_id": candidate or f"c{rank}",
        "candidate_rank": rank,
        "final_score_tier": tier,
    }


def test_bookmark_name_is_deterministic() -> None:
    assert bookmark_name("q1", "cand_1") == bookmark_name("q1", "cand_1")


def test_bookmark_name_differs_for_different_inputs() -> None:
    assert bookmark_name("q1", "cand_1") != bookmark_name("q1", "cand_2")
    assert bookmark_name("q1", "cand_1") != bookmark_name("q2", "cand_1")


def test_bookmark_name_is_word_legal() -> None:
    name = bookmark_name("MA_4115", "WP_012345678.1")
    assert name[0].isalpha() or name[0] == "_"
    assert all(ch.isalnum() or ch == "_" for ch in name)
    assert len(name) <= 40


def test_bookmark_name_long_ids_do_not_collide() -> None:
    long_id_a = "protein_" + "a" * 60
    long_id_b = "protein_" + "a" * 59 + "b"
    assert bookmark_name("q1", long_id_a) != bookmark_name("q1", long_id_b)


def test_select_top_candidates_per_query_keeps_only_top_n() -> None:
    rows = [_ranked_row("q1", rank) for rank in range(1, 6)]
    selected = select_top_candidates_per_query(rows, max_per_query=3)
    assert [row["candidate_protein_id"] for row in selected] == ["c1", "c2", "c3"]


def test_select_top_candidates_per_query_includes_safety_net_beyond_n() -> None:
    rows = [_ranked_row("q1", rank) for rank in range(1, 6)]
    rows[4]["final_score_tier"] = "Tier1_VeryStrong"  # rank 5, beyond max_per_query=3
    selected = select_top_candidates_per_query(rows, max_per_query=3)
    assert [row["candidate_protein_id"] for row in selected] == ["c1", "c2", "c3", "c5"]


def test_select_top_candidates_per_query_excludes_no_query_rows() -> None:
    rows = [_ranked_row("", 1)]
    assert select_top_candidates_per_query(rows, max_per_query=15) == []


def test_select_top_candidates_per_query_excludes_ineligible_rows() -> None:
    rows = [_ranked_row("q1", 0), _ranked_row("q1", 1)]
    selected = select_top_candidates_per_query(rows, max_per_query=15)
    assert [row["candidate_rank"] for row in selected] == [1]


def test_select_top_candidates_per_query_is_per_query_independent() -> None:
    rows = [_ranked_row("q1", 1), _ranked_row("q1", 2), _ranked_row("q2", 1), _ranked_row("q2", 2)]
    selected = select_top_candidates_per_query(rows, max_per_query=1)
    assert [(row["query_id"], row["candidate_rank"]) for row in selected] == [("q1", 1), ("q2", 1)]


def test_tier_safety_net_contains_only_top_two_tiers() -> None:
    assert TIER_SAFETY_NET == {"Tier1_VeryStrong", "Tier2_Strong"}
