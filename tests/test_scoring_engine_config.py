"""Tests for analysis/scoring_engine_config.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from core.exceptions import ConfigError
from analysis.scoring_engine_config import (
    DEFAULT_CATEGORY_CAPS,
    DEFAULT_SCORING_ENGINE_CONFIG,
    load_scoring_engine_config,
)

#: config/scoring_engine.example.yaml is documentation a user is expected to
#: copy and edit; it is never loaded automatically. If it silently drifts
#: out of sync with DEFAULT_CATEGORY_CAPS (e.g. a new evidence category is
#: added to the code but not to the example file), a user who copies the
#: example verbatim as their scoring_engine_config would hit a ConfigError
#: the first time that category actually fires for a real pair -- see
#: analysis/scoring_engine.py::_score_categories. This path matches the
#: layout documented in <repo_root>/config/scoring_engine.example.yaml.
_EXAMPLE_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config" / "scoring_engine.example.yaml"


def test_default_config_used_when_no_path() -> None:
    config = load_scoring_engine_config(None)
    assert config is DEFAULT_SCORING_ENGINE_CONFIG
    assert config.category_caps["source_classification"] == 30.0
    assert config.negative_penalty_cap == 30.0


def test_missing_file_raises_config_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigError):
        load_scoring_engine_config(tmp_path / "does_not_exist.yaml")


def test_invalid_yaml_raises_config_error(tmp_path: Path) -> None:
    path = tmp_path / "bad.yaml"
    path.write_text("category_caps: [unterminated", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_scoring_engine_config(path)


def test_non_mapping_yaml_raises_config_error(tmp_path: Path) -> None:
    path = tmp_path / "list.yaml"
    path.write_text("- a\n- b\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_scoring_engine_config(path)


def test_valid_custom_config_is_parsed(tmp_path: Path) -> None:
    path = tmp_path / "scoring.yaml"
    path.write_text(
        """
output_scale: 100
category_caps:
  source_classification: 40
  genomic_context: 20
  functional_annotation: 15
negative_penalty_cap: 20
minimum_evidence:
  min_categories: 2
  min_available_weight: 5
tiers:
  tier1_min_score: 80
  tier1_min_categories: 3
  tier2_min_score: 55
  tier2_min_categories: 2
  tier3_min_score: 30
  tier3_min_categories: 1
tie_precision: 2
""",
        encoding="utf-8",
    )
    config = load_scoring_engine_config(path)
    assert config.output_scale == 100.0
    assert config.category_caps == {
        "source_classification": 40.0,
        "genomic_context": 20.0,
        "functional_annotation": 15.0,
    }
    assert config.negative_penalty_cap == 20.0
    assert config.minimum_evidence.min_categories == 2
    assert config.minimum_evidence.min_available_weight == 5.0
    assert config.tiers.tier1_min_score == 80.0
    assert config.tie_precision == 2


def test_null_negative_penalty_cap_means_uncapped(tmp_path: Path) -> None:
    path = tmp_path / "scoring.yaml"
    path.write_text(
        "category_caps:\n  source_classification: 30\nnegative_penalty_cap: null\n",
        encoding="utf-8",
    )
    config = load_scoring_engine_config(path)
    assert config.negative_penalty_cap is None


def test_empty_category_caps_rejected(tmp_path: Path) -> None:
    path = tmp_path / "scoring.yaml"
    path.write_text("category_caps: {}\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_scoring_engine_config(path)


def test_non_numeric_cap_rejected(tmp_path: Path) -> None:
    path = tmp_path / "scoring.yaml"
    path.write_text("category_caps:\n  source_classification: not_a_number\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_scoring_engine_config(path)


def test_example_config_matches_defaults() -> None:
    """config/scoring_engine.example.yaml must stay in sync with the code.

    This guards against exactly the kind of drift found during self-review:
    a new evidence category (e.g. the ProteinInteractionHunter bridge's
    pih_* categories) added to DEFAULT_CATEGORY_CAPS in code but forgotten
    in the example file a user is expected to copy.
    """
    assert _EXAMPLE_CONFIG_PATH.exists(), (
        f"expected example scoring engine config at {_EXAMPLE_CONFIG_PATH}"
    )
    config = load_scoring_engine_config(_EXAMPLE_CONFIG_PATH)
    assert config.category_caps == DEFAULT_CATEGORY_CAPS
    assert config.negative_penalty_cap == DEFAULT_SCORING_ENGINE_CONFIG.negative_penalty_cap
    assert config.tiers == DEFAULT_SCORING_ENGINE_CONFIG.tiers
    assert config.minimum_evidence == DEFAULT_SCORING_ENGINE_CONFIG.minimum_evidence
    assert config.tie_precision == DEFAULT_SCORING_ENGINE_CONFIG.tie_precision
