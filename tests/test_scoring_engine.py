"""Tests for analysis/scoring_engine.py."""

from __future__ import annotations

import pytest

from core.evidence import EvidenceComponent
from analysis.scoring_engine import rank_candidates, score_candidate, UNCLASSIFIED_TIER
from analysis.scoring_engine_config import (
    MinimumEvidenceConfig,
    ScoringEngineConfig,
    TierThresholds,
)
from core.exceptions import ConfigError


def make_engine_config(**overrides) -> ScoringEngineConfig:
    defaults = dict(
        output_scale=100.0,
        category_caps={
            "source_classification": 30.0,
            "genomic_context": 25.0,
            "functional_annotation": 20.0,
        },
        negative_penalty_cap=30.0,
        minimum_evidence=MinimumEvidenceConfig(min_categories=1, min_available_weight=0.0),
        tiers=TierThresholds(),
        tie_precision=3,
    )
    defaults.update(overrides)
    return ScoringEngineConfig(**defaults)


def test_full_evidence_scores_near_output_scale() -> None:
    components = [
        EvidenceComponent.available("source", "source_classification", 1.0, 30.0),
        EvidenceComponent.available("distance", "genomic_context", 1.0, 25.0),
        EvidenceComponent.available("domain", "functional_annotation", 1.0, 20.0),
    ]
    breakdown = score_candidate(components, make_engine_config())
    assert breakdown.eligible is True
    assert breakdown.final_score == pytest.approx(100.0)
    assert breakdown.evidence_category_count == 3
    assert breakdown.tier == "Tier1_VeryStrong"


def test_zero_evidence_scores_zero_not_none_when_evaluated() -> None:
    components = [
        EvidenceComponent.available("source", "source_classification", 0.0, 30.0),
    ]
    breakdown = score_candidate(components, make_engine_config())
    assert breakdown.eligible is True
    assert breakdown.final_score == pytest.approx(0.0)
    assert breakdown.tier == "Tier4_Weak"


def test_missing_evidence_is_excluded_from_denominator_not_zeroed() -> None:
    # Only genomic_context is evaluated; functional_annotation was never run.
    from core.evidence import EvidenceStatus

    components = [
        EvidenceComponent.available("distance", "genomic_context", 0.5, 25.0),
        EvidenceComponent.unavailable("domain", "functional_annotation", EvidenceStatus.MISSING),
    ]
    breakdown = score_candidate(components, make_engine_config())
    # Only genomic_context (cap 25) is active, so total_cap == 25, not 45+.
    assert breakdown.total_cap == pytest.approx(25.0)
    assert breakdown.evidence_category_count == 1
    assert breakdown.final_score == pytest.approx(50.0)  # 0.5 normalized * 100


def test_category_cap_prevents_double_counting() -> None:
    # Two components share functional_annotation; even both maxed out must
    # not exceed that category's cap.
    components = [
        EvidenceComponent.available("co_occurrence", "functional_annotation", 1.0, 10.0),
        EvidenceComponent.available("domain_complementarity", "functional_annotation", 1.0, 10.0),
    ]
    breakdown = score_candidate(components, make_engine_config())
    assert breakdown.category_scores["functional_annotation"].capped_score == pytest.approx(20.0)
    assert breakdown.final_score == pytest.approx(100.0)  # only category active, capped at 20/20


def test_negative_evidence_applies_capped_penalty() -> None:
    components = [
        EvidenceComponent.available("source", "source_classification", 1.0, 30.0),
        EvidenceComponent.available(
            "incompatible_localization",
            "cellular_compatibility",
            1.0,
            60.0,
            is_negative=True,
        ),
    ]
    breakdown = score_candidate(components, make_engine_config(negative_penalty_cap=15.0))
    assert breakdown.negative_penalty_points == pytest.approx(15.0)
    assert breakdown.final_score == pytest.approx(100.0 - 15.0)


def test_insufficient_evidence_yields_no_formal_score() -> None:
    components = [
        EvidenceComponent.available("distance", "genomic_context", 0.9, 25.0),
    ]
    engine_config = make_engine_config(
        minimum_evidence=MinimumEvidenceConfig(min_categories=2, min_available_weight=0.0)
    )
    breakdown = score_candidate(components, engine_config)
    assert breakdown.eligible is False
    assert breakdown.final_score is None
    assert breakdown.tier == UNCLASSIFIED_TIER


def test_unknown_category_without_cap_raises_config_error() -> None:
    components = [EvidenceComponent.available("mystery", "not_configured", 1.0, 10.0)]
    with pytest.raises(ConfigError):
        score_candidate(components, make_engine_config())


def test_ranking_orders_by_score_then_categories_then_weight_then_id() -> None:
    high = score_candidate(
        [EvidenceComponent.available("s", "source_classification", 1.0, 30.0)],
        make_engine_config(),
    )
    tie_a = score_candidate(
        [
            EvidenceComponent.available("s", "source_classification", 0.5, 30.0),
            EvidenceComponent.available("g", "genomic_context", 0.5, 25.0),
        ],
        make_engine_config(),
    )
    tie_b_same_score_fewer_categories = score_candidate(
        [EvidenceComponent.available("s", "source_classification", 0.5, 30.0)],
        make_engine_config(),
    )
    ranked = rank_candidates(
        [
            ("candidate_high", high),
            ("candidate_tie_a", tie_a),
            ("candidate_tie_b", tie_b_same_score_fewer_categories),
        ]
    )
    order = [row.candidate_id for row in ranked]
    # tie_a is listed before tie_b (more evidence categories orders display),
    # but both share the same dense rank because secondary keys stabilize
    # display order without breaking a score tie (design spec, section 21).
    assert order == ["candidate_high", "candidate_tie_a", "candidate_tie_b"]
    assert ranked[0].rank == 1
    assert ranked[1].rank == 2
    assert ranked[2].rank == 2


def test_equal_score_and_categories_use_dense_rank() -> None:
    a = score_candidate(
        [EvidenceComponent.available("s", "source_classification", 1.0, 30.0)],
        make_engine_config(),
    )
    b = score_candidate(
        [EvidenceComponent.available("s", "source_classification", 1.0, 30.0)],
        make_engine_config(),
    )
    c = score_candidate(
        [EvidenceComponent.available("s", "source_classification", 0.0, 30.0)],
        make_engine_config(),
    )
    ranked = rank_candidates([("b_candidate", b), ("a_candidate", a), ("c_candidate", c)])
    ranks = {row.candidate_id: row.rank for row in ranked}
    # a_candidate and b_candidate tie for rank 1 (alphabetical tie-break for
    # display order only); c_candidate gets dense rank 2, not 3.
    assert ranks["a_candidate"] == 1
    assert ranks["b_candidate"] == 1
    assert ranks["c_candidate"] == 2


def test_ineligible_candidates_get_no_numeric_rank_and_sort_last() -> None:
    eligible = score_candidate(
        [EvidenceComponent.available("s", "source_classification", 0.1, 30.0)],
        make_engine_config(),
    )
    ineligible = score_candidate(
        [EvidenceComponent.available("g", "genomic_context", 0.9, 25.0)],
        make_engine_config(
            minimum_evidence=MinimumEvidenceConfig(min_categories=5, min_available_weight=0.0)
        ),
    )
    ranked = rank_candidates([("weak_but_eligible", eligible), ("strong_but_ineligible", ineligible)])
    assert [row.candidate_id for row in ranked] == ["weak_but_eligible", "strong_but_ineligible"]
    assert ranked[0].rank == 1
    assert ranked[1].rank is None
