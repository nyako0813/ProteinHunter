"""Tests for analysis/functional_complementarity_rules.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from core.exceptions import ConfigError
from analysis.functional_complementarity_rules import (
    DEFAULT_RULESET_PATH,
    load_functional_complementarity_ruleset,
)


def test_default_ruleset_loads_and_has_known_rule() -> None:
    ruleset = load_functional_complementarity_ruleset()
    assert ruleset.version == "v1"
    assert any(rule.rule_id == "radical_sam_iron_sulfur" for rule in ruleset.rules)
    assert "atpase" in ruleset.meaningful_keywords
    assert "hypothetical" in ruleset.stopwords


def test_default_ruleset_path_exists() -> None:
    assert DEFAULT_RULESET_PATH.exists()


def test_rule_matches_either_orientation() -> None:
    ruleset = load_functional_complementarity_ruleset()
    rule = next(r for r in ruleset.rules if r.rule_id == "radical_sam_iron_sulfur")
    assert rule.matches({"radical sam"}, {"iron-sulfur"}) is True
    assert rule.matches({"iron-sulfur"}, {"radical sam"}) is True
    assert rule.matches({"radical sam"}, {"nad"}) is False


def test_find_match_returns_first_hit() -> None:
    ruleset = load_functional_complementarity_ruleset()
    match = ruleset.find_match({"radical sam"}, {"iron-sulfur"})
    assert match is not None
    assert match.rule_id == "radical_sam_iron_sulfur"


def test_find_match_returns_none_when_nothing_matches() -> None:
    ruleset = load_functional_complementarity_ruleset()
    assert ruleset.find_match({"unrelated"}, {"also-unrelated"}) is None


def test_missing_file_raises_config_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigError):
        load_functional_complementarity_ruleset(tmp_path / "missing.yaml")


def test_duplicate_rule_id_rejected(tmp_path: Path) -> None:
    path = tmp_path / "rules.yaml"
    path.write_text(
        """
version: test
rules:
  - rule_id: dup
    left_terms: [a]
    right_terms: [b]
  - rule_id: dup
    left_terms: [c]
    right_terms: [d]
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigError):
        load_functional_complementarity_ruleset(path)


def test_rule_missing_terms_rejected(tmp_path: Path) -> None:
    path = tmp_path / "rules.yaml"
    path.write_text(
        """
version: test
rules:
  - rule_id: only_left
    left_terms: [a]
    right_terms: []
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigError):
        load_functional_complementarity_ruleset(path)


def test_custom_ruleset_overrides_default(tmp_path: Path) -> None:
    path = tmp_path / "rules.yaml"
    path.write_text(
        """
version: custom
rules:
  - rule_id: custom_pair
    left_terms: [foo]
    right_terms: [bar]
meaningful_keywords: [foo, bar]
stopwords: [protein]
""",
        encoding="utf-8",
    )
    ruleset = load_functional_complementarity_ruleset(path)
    assert ruleset.version == "custom"
    assert len(ruleset.rules) == 1
    assert ruleset.rules[0].rule_id == "custom_pair"
    assert ruleset.meaningful_keywords == frozenset({"foo", "bar"})
