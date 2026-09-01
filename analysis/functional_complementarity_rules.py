"""Loader for the functional complementarity ruleset (scoring model v2).

The rule pairs, keyword list, and stopword list used to be hard-coded
Python constants in ``analysis/interaction_scoring.py``
(``COMPLEMENTARY_TERM_PAIRS``, ``MEANINGFUL_KEYWORDS``,
``DESCRIPTION_STOPWORDS``). Those constants still exist and are unchanged,
so the original ("legacy_additive") scoring path behaves exactly as before.
This module lets the new evidence-based ("v2_evidence_based") path load the
same rules -- or a project/organism-specific override -- from a versioned
YAML file instead, per design specification section 3.4.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from core.exceptions import ConfigError


@dataclass(slots=True, frozen=True)
class ComplementarityRule:
    """One directional-or-symmetric term-pair rule."""

    rule_id: str
    left_terms: frozenset[str]
    right_terms: frozenset[str]
    note: str = ""

    def matches(self, left_present: set[str], right_present: set[str]) -> bool:
        """Return True if this rule fires for either orientation of terms."""
        forward = bool(self.left_terms & left_present) and bool(self.right_terms & right_present)
        backward = bool(self.right_terms & left_present) and bool(self.left_terms & right_present)
        return forward or backward


@dataclass(slots=True, frozen=True)
class FunctionalComplementarityRuleset:
    """A versioned collection of complementarity rules plus vocabulary."""

    version: str
    rules: tuple[ComplementarityRule, ...]
    meaningful_keywords: frozenset[str]
    stopwords: frozenset[str]

    def find_match(
        self, left_terms: set[str], right_terms: set[str]
    ) -> ComplementarityRule | None:
        """Return the first rule whose terms match, or None."""
        for rule in self.rules:
            if rule.matches(left_terms, right_terms):
                return rule
        return None


DEFAULT_RULESET_PATH = (
    Path(__file__).resolve().parent.parent / "config" / "functional_complementarity_rules.v1.yaml"
)


def load_functional_complementarity_ruleset(
    path: Path | None = None,
) -> FunctionalComplementarityRuleset:
    """Load a ruleset from ``path``, defaulting to the shipped v1 ruleset."""
    resolved_path = path or DEFAULT_RULESET_PATH
    if not resolved_path.exists():
        raise ConfigError(
            f"functional complementarity ruleset file was not found: {resolved_path}"
        )

    try:
        raw = yaml.safe_load(resolved_path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ConfigError(
            f"functional complementarity ruleset is not valid YAML: {resolved_path} ({exc})"
        ) from exc

    if not isinstance(raw, dict):
        raise ConfigError(
            f"functional complementarity ruleset must contain a mapping: {resolved_path}"
        )

    version = str(raw.get("version", "unknown"))

    raw_rules = raw.get("rules", [])
    if not isinstance(raw_rules, list):
        raise ConfigError(f"'rules' in {resolved_path} must be a list.")

    rules: list[ComplementarityRule] = []
    seen_ids: set[str] = set()
    for index, raw_rule in enumerate(raw_rules):
        if not isinstance(raw_rule, dict):
            raise ConfigError(f"'rules[{index}]' in {resolved_path} must be a mapping.")
        rule_id = str(raw_rule.get("rule_id") or "").strip()
        if not rule_id:
            raise ConfigError(f"'rules[{index}].rule_id' in {resolved_path} is required.")
        if rule_id in seen_ids:
            raise ConfigError(f"duplicate rule_id '{rule_id}' in {resolved_path}.")
        seen_ids.add(rule_id)

        left_terms = _string_list(raw_rule.get("left_terms"), f"rules[{index}].left_terms", resolved_path)
        right_terms = _string_list(raw_rule.get("right_terms"), f"rules[{index}].right_terms", resolved_path)
        if not left_terms or not right_terms:
            raise ConfigError(
                f"'rules[{index}]' ({rule_id}) in {resolved_path} needs at least one "
                "left_terms entry and one right_terms entry."
            )
        rules.append(
            ComplementarityRule(
                rule_id=rule_id,
                left_terms=frozenset(left_terms),
                right_terms=frozenset(right_terms),
                note=str(raw_rule.get("note", "")),
            )
        )

    meaningful_keywords = frozenset(
        _string_list(raw.get("meaningful_keywords"), "meaningful_keywords", resolved_path)
    )
    stopwords = frozenset(_string_list(raw.get("stopwords"), "stopwords", resolved_path))

    return FunctionalComplementarityRuleset(
        version=version,
        rules=tuple(rules),
        meaningful_keywords=meaningful_keywords,
        stopwords=stopwords,
    )


def _string_list(value: object, field_name: str, path: Path) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ConfigError(f"'{field_name}' in {path} must be a list.")
    result: list[str] = []
    for item in value:
        if not isinstance(item, str) or not item:
            raise ConfigError(f"'{field_name}' in {path} must contain only non-empty strings.")
        result.append(item.lower())
    return result


__all__: tuple[str, ...] = (
    "ComplementarityRule",
    "DEFAULT_RULESET_PATH",
    "FunctionalComplementarityRuleset",
    "load_functional_complementarity_ruleset",
)
