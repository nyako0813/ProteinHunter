"""Evidence-based scoring engine (scoring model v2).

Replaces plain "add up fixed points" scoring with the model described in
the project design specification (sections 11-24): components are grouped
into categories, each category is normalized against only the evidence that
was actually evaluated, category contributions are capped so correlated
evidence cannot be double counted, negative (contradicting) evidence is
applied as a separate, capped penalty, and every candidate keeps a full,
auditable breakdown instead of a single opaque number.

This module has no knowledge of BLAST, GFF, or Excel. It only knows how to
turn a list of :class:`~core.evidence.EvidenceComponent` into a
:class:`ScoreBreakdown`. ``analysis/interaction_scoring.py`` is responsible
for building the components from ProteinHunter's actual evidence.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from core.evidence import EvidenceComponent, EvidenceStatus
from analysis.scoring_engine_config import ScoringEngineConfig
from core.exceptions import ConfigError


UNCLASSIFIED_TIER = "Unclassified"
TIER_LABELS: tuple[str, ...] = (
    "Tier1_VeryStrong",
    "Tier2_Strong",
    "Tier3_Moderate",
    "Tier4_Weak",
    UNCLASSIFIED_TIER,
)


@dataclass(slots=True, frozen=True)
class CategoryScore:
    """Aggregated, capped score for one evidence category."""

    category: str
    available_weight: float
    raw_weighted_sum: float
    normalized_score: float  # 0.0-1.0, within-category
    cap: float
    capped_score: float  # 0.0-cap, in output_scale units
    component_count: int


@dataclass(slots=True, frozen=True)
class ScoreBreakdown:
    """Full, auditable result of scoring one candidate pair."""

    components: tuple[EvidenceComponent, ...]
    category_scores: dict[str, CategoryScore]
    positive_raw_total: float
    total_cap: float
    negative_penalty_points: float
    final_score: float | None  # None when evidence is insufficient
    evidence_category_count: int
    evidence_component_count: int
    available_weight_total: float
    tier: str
    eligible: bool


def score_candidate(
    components: list[EvidenceComponent], engine_config: ScoringEngineConfig
) -> ScoreBreakdown:
    """Score one candidate pair from its evidence components.

    Only components with ``status == AVAILABLE`` contribute. Positive
    components are aggregated per category (capped); negative-flagged
    components are summed separately into a capped penalty. Evidence
    completeness (category count, available weight) gates whether a formal
    ``final_score`` is produced at all -- an ineligible candidate keeps its
    breakdown (for audit) but sorts after every eligible one.
    """
    positive = [c for c in components if c.status is EvidenceStatus.AVAILABLE and not c.is_negative]
    negative = [c for c in components if c.status is EvidenceStatus.AVAILABLE and c.is_negative]

    category_scores = _score_categories(positive, engine_config.category_caps)
    active_categories = {
        category: score for category, score in category_scores.items() if score.available_weight > 0
    }

    total_cap = sum(score.cap for score in active_categories.values())
    positive_raw_total = sum(score.capped_score for score in active_categories.values())

    available_weight_total = sum(score.available_weight for score in active_categories.values())
    evidence_category_count = len(active_categories)
    evidence_component_count = len(positive) + len(negative)

    minimum = engine_config.minimum_evidence
    eligible = (
        evidence_category_count >= minimum.min_categories
        and available_weight_total >= minimum.min_available_weight
        and total_cap > 0
    )

    negative_penalty_points = _negative_penalty(negative, engine_config.negative_penalty_cap)

    final_score: float | None
    if not eligible:
        final_score = None
    else:
        raw_score = positive_raw_total / total_cap * engine_config.output_scale
        final_score = _clamp(raw_score - negative_penalty_points, 0.0, engine_config.output_scale)

    tier = _classify_tier(final_score, evidence_category_count, engine_config.tiers)

    return ScoreBreakdown(
        components=tuple(components),
        category_scores=category_scores,
        positive_raw_total=positive_raw_total,
        total_cap=total_cap,
        negative_penalty_points=negative_penalty_points,
        final_score=final_score,
        evidence_category_count=evidence_category_count,
        evidence_component_count=evidence_component_count,
        available_weight_total=available_weight_total,
        tier=tier,
        eligible=eligible,
    )


def _score_categories(
    positive_components: list[EvidenceComponent], category_caps: dict[str, float]
) -> dict[str, CategoryScore]:
    by_category: dict[str, list[EvidenceComponent]] = {}
    for component in positive_components:
        by_category.setdefault(component.category, []).append(component)

    scores: dict[str, CategoryScore] = {}
    for category, category_components in by_category.items():
        if category not in category_caps:
            raise ConfigError(
                f"evidence category '{category}' has no configured cap in "
                "category_caps. Add it to the scoring engine config, or fix "
                "the category name used when building evidence components."
            )
        available_weight = sum(c.effective_weight for c in category_components)
        raw_weighted_sum = sum(c.contribution for c in category_components)
        normalized_score = (
            _clamp(raw_weighted_sum / available_weight, 0.0, 1.0) if available_weight > 0 else 0.0
        )
        cap = category_caps[category]
        scores[category] = CategoryScore(
            category=category,
            available_weight=available_weight,
            raw_weighted_sum=raw_weighted_sum,
            normalized_score=normalized_score,
            cap=cap,
            capped_score=normalized_score * cap,
            component_count=len(category_components),
        )

    # Categories with a configured cap but zero available evidence this run
    # are still reported (available_weight == 0) so the breakdown always
    # shows every configured category, not just the ones that fired.
    for category, cap in category_caps.items():
        if category not in scores:
            scores[category] = CategoryScore(
                category=category,
                available_weight=0.0,
                raw_weighted_sum=0.0,
                normalized_score=0.0,
                cap=cap,
                capped_score=0.0,
                component_count=0,
            )

    return scores


def _negative_penalty(
    negative_components: list[EvidenceComponent], penalty_cap: float | None
) -> float:
    if not negative_components:
        return 0.0
    raw_penalty = sum(c.contribution for c in negative_components)
    if penalty_cap is None:
        return raw_penalty
    return min(raw_penalty, penalty_cap)


def _classify_tier(
    final_score: float | None, category_count: int, thresholds
) -> str:
    if final_score is None:
        return UNCLASSIFIED_TIER
    if final_score >= thresholds.tier1_min_score and category_count >= thresholds.tier1_min_categories:
        return "Tier1_VeryStrong"
    if final_score >= thresholds.tier2_min_score and category_count >= thresholds.tier2_min_categories:
        return "Tier2_Strong"
    if final_score >= thresholds.tier3_min_score and category_count >= thresholds.tier3_min_categories:
        return "Tier3_Moderate"
    return "Tier4_Weak"


def _clamp(value: float, low: float, high: float) -> float:
    if value < low:
        return low
    if value > high:
        return high
    return value


@dataclass(slots=True, frozen=True)
class RankedCandidate:
    """One row of a deterministic per-query ranking.

    ``rank`` is ``None`` for candidates whose evidence was insufficient for
    a formal score (``ScoreBreakdown.eligible is False``); they are still
    returned, sorted after every ranked candidate, so nothing is silently
    dropped from the audit trail.
    """

    candidate_id: str
    breakdown: ScoreBreakdown
    rank: int | None


def rank_candidates(
    candidates: list[tuple[str, ScoreBreakdown]], tie_precision: int = 3
) -> list[RankedCandidate]:
    """Rank candidates deterministically for one query.

    Sort order: eligible before ineligible, final score descending,
    evidence category count descending, available evidence weight
    descending, candidate ID ascending. Equal scores (after rounding to
    ``tie_precision``) receive the same dense rank (1, 1, 2, ...), matching
    the tie-break rule in the design specification (section 21). Ineligible
    candidates (no formal score) are never assigned a numeric rank.
    """

    def sort_key(item: tuple[str, ScoreBreakdown]) -> tuple[bool, float, float, float, str]:
        candidate_id, breakdown = item
        score = breakdown.final_score
        score_key = -round(score, tie_precision) if score is not None else 0.0
        return (
            not breakdown.eligible,
            score_key,
            -breakdown.evidence_category_count,
            -breakdown.available_weight_total,
            candidate_id,
        )

    ordered = sorted(candidates, key=sort_key)

    ranked: list[RankedCandidate] = []
    previous_score_key: float | None = None
    current_rank = 0
    for candidate_id, breakdown in ordered:
        if not breakdown.eligible:
            ranked.append(RankedCandidate(candidate_id=candidate_id, breakdown=breakdown, rank=None))
            continue
        score_key = round(breakdown.final_score, tie_precision)
        if previous_score_key is None or score_key != previous_score_key:
            current_rank += 1
        previous_score_key = score_key
        ranked.append(RankedCandidate(candidate_id=candidate_id, breakdown=breakdown, rank=current_rank))

    return ranked


__all__: tuple[str, ...] = (
    "CategoryScore",
    "RankedCandidate",
    "ScoreBreakdown",
    "TIER_LABELS",
    "UNCLASSIFIED_TIER",
    "rank_candidates",
    "score_candidate",
)
