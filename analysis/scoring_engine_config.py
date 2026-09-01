"""Configuration for the evidence-based scoring engine (scoring model v2).

This is intentionally a small, standalone YAML file rather than another
section inside the already large ``config.py`` validator. Weights, caps, and
penalties are exactly the kind of numbers that need to be retuned once real
calibration data (known interacting / non-interacting pairs) is available,
so they must never be hard-coded in Python (see the project design
specification, sections 14 and 36).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

from core.exceptions import ConfigError


@dataclass(slots=True, frozen=True)
class TierThresholds:
    """Score/category-count thresholds for the four confidence tiers.

    A candidate that does not reach ``tier3`` is still reported (as
    ``Tier4_Weak``); a candidate with no formal score at all (insufficient
    evidence) is reported as ``Unclassified``, never silently dropped.
    """

    tier1_min_score: float = 70.0
    tier1_min_categories: int = 3
    tier2_min_score: float = 50.0
    tier2_min_categories: int = 2
    tier3_min_score: float = 25.0
    tier3_min_categories: int = 1


@dataclass(slots=True, frozen=True)
class MinimumEvidenceConfig:
    """Eligibility gate for producing a formal (ranked) score at all."""

    min_categories: int = 1
    min_available_weight: float = 0.0


@dataclass(slots=True, frozen=True)
class ScoringEngineConfig:
    """Everything the scoring engine needs that is not biology-specific."""

    output_scale: float = 100.0
    category_caps: dict[str, float] = field(default_factory=dict)
    negative_penalty_cap: float | None = 30.0
    minimum_evidence: MinimumEvidenceConfig = field(default_factory=MinimumEvidenceConfig)
    tiers: TierThresholds = field(default_factory=TierThresholds)
    tie_precision: int = 3


#: Category caps that reproduce the point budget of the legacy additive
#: scorer (candidate_priority=30, gene_neighborhood=25,
#: co_occurrence+domain_complementarity=20 combined instead of 20+15=35, to
#: reduce double counting -- see the design document, section 3.1/3.2).
#:
#: The three "pih_*" caps are always present but only ever active for a
#: pair when the optional ProteinInteractionHunter evidence bridge
#: (analysis/pih_evidence_bridge.py) is configured and actually finds
#: matching evidence; a run that never uses the bridge is unaffected,
#: because a category with no available evidence contributes nothing and
#: is excluded from the score denominator. Values must match
#: analysis.pih_evidence_bridge.BRIDGED_PIH_CATEGORY_CAPS.
DEFAULT_CATEGORY_CAPS: dict[str, float] = {
    "source_classification": 30.0,
    "genomic_context": 25.0,
    "functional_annotation": 20.0,
    "pih_cellular_compatibility": 5.0,
    "pih_evolutionary": 10.0,
    "pih_direct_interaction": 20.0,
}

DEFAULT_SCORING_ENGINE_CONFIG = ScoringEngineConfig(
    category_caps=dict(DEFAULT_CATEGORY_CAPS)
)


def load_scoring_engine_config(path: Path | None) -> ScoringEngineConfig:
    """Load a :class:`ScoringEngineConfig` from ``path``, or return defaults."""
    if path is None:
        return DEFAULT_SCORING_ENGINE_CONFIG

    if not path.exists():
        raise ConfigError(
            f"scoring engine config file was not found: {path}. "
            "Either create it or remove 'scoring_engine_config' from "
            "config.yaml to use the built-in defaults."
        )

    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ConfigError(
            f"scoring engine config file is not valid YAML: {path} ({exc})"
        ) from exc

    if not isinstance(raw, dict):
        raise ConfigError(
            f"scoring engine config file must contain a YAML mapping: {path}"
        )

    return _parse_scoring_engine_config(raw, path)


def _parse_scoring_engine_config(raw: dict[object, object], path: Path) -> ScoringEngineConfig:
    output_scale = _positive_float(raw.get("output_scale", 100.0), "output_scale", path)

    raw_caps = raw.get("category_caps", DEFAULT_CATEGORY_CAPS)
    if not isinstance(raw_caps, dict) or not raw_caps:
        raise ConfigError(
            f"'category_caps' in {path} must be a non-empty mapping of "
            "category name to a positive number."
        )
    category_caps: dict[str, float] = {}
    for category, cap in raw_caps.items():
        if not isinstance(category, str) or not category:
            raise ConfigError(f"'category_caps' keys in {path} must be non-empty strings.")
        category_caps[category] = _positive_float(cap, f"category_caps.{category}", path)

    negative_penalty_cap_raw = raw.get("negative_penalty_cap", 30.0)
    negative_penalty_cap: float | None
    if negative_penalty_cap_raw is None:
        negative_penalty_cap = None
    else:
        negative_penalty_cap = _positive_float(
            negative_penalty_cap_raw, "negative_penalty_cap", path, allow_zero=True
        )

    minimum_evidence_raw = raw.get("minimum_evidence", {}) or {}
    if not isinstance(minimum_evidence_raw, dict):
        raise ConfigError(f"'minimum_evidence' in {path} must be a mapping.")
    minimum_evidence = MinimumEvidenceConfig(
        min_categories=_positive_int(
            minimum_evidence_raw.get("min_categories", 1), "minimum_evidence.min_categories", path
        ),
        min_available_weight=_positive_float(
            minimum_evidence_raw.get("min_available_weight", 0.0),
            "minimum_evidence.min_available_weight",
            path,
            allow_zero=True,
        ),
    )

    tiers_raw = raw.get("tiers", {}) or {}
    if not isinstance(tiers_raw, dict):
        raise ConfigError(f"'tiers' in {path} must be a mapping.")
    tiers = TierThresholds(
        tier1_min_score=_positive_float(
            tiers_raw.get("tier1_min_score", 70.0), "tiers.tier1_min_score", path
        ),
        tier1_min_categories=_positive_int(
            tiers_raw.get("tier1_min_categories", 3), "tiers.tier1_min_categories", path
        ),
        tier2_min_score=_positive_float(
            tiers_raw.get("tier2_min_score", 50.0), "tiers.tier2_min_score", path
        ),
        tier2_min_categories=_positive_int(
            tiers_raw.get("tier2_min_categories", 2), "tiers.tier2_min_categories", path
        ),
        tier3_min_score=_positive_float(
            tiers_raw.get("tier3_min_score", 25.0), "tiers.tier3_min_score", path
        ),
        tier3_min_categories=_positive_int(
            tiers_raw.get("tier3_min_categories", 1), "tiers.tier3_min_categories", path
        ),
    )

    tie_precision = _positive_int(raw.get("tie_precision", 3), "tie_precision", path, allow_zero=True)

    return ScoringEngineConfig(
        output_scale=output_scale,
        category_caps=category_caps,
        negative_penalty_cap=negative_penalty_cap,
        minimum_evidence=minimum_evidence,
        tiers=tiers,
        tie_precision=tie_precision,
    )


def _positive_float(
    value: object, field_name: str, path: Path, *, allow_zero: bool = False
) -> float:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"'{field_name}' in {path} must be a number.") from exc
    if allow_zero and number < 0:
        raise ConfigError(f"'{field_name}' in {path} must be >= 0.")
    if not allow_zero and number <= 0:
        raise ConfigError(f"'{field_name}' in {path} must be a positive number.")
    return number


def _positive_int(
    value: object, field_name: str, path: Path, *, allow_zero: bool = False
) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ConfigError(f"'{field_name}' in {path} must be a whole number.")
    if allow_zero and value < 0:
        raise ConfigError(f"'{field_name}' in {path} must be >= 0.")
    if not allow_zero and value <= 0:
        raise ConfigError(f"'{field_name}' in {path} must be a positive whole number.")
    return value


__all__: tuple[str, ...] = (
    "DEFAULT_CATEGORY_CAPS",
    "DEFAULT_SCORING_ENGINE_CONFIG",
    "MinimumEvidenceConfig",
    "ScoringEngineConfig",
    "TierThresholds",
    "load_scoring_engine_config",
)
