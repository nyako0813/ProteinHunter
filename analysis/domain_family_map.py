"""Loader for the domain family map (scoring model v2, UniProt Pfam/InterPro-based).

Coarse domain-family categories (e.g. "electron_carrier",
"atp_dependent_activator") group Pfam/InterPro/SUPFAM accessions so
``_domain_complementarity_status_and_value`` (analysis/interaction_scoring.py)
can recognize a functional relationship even when the two proteins' NCBI
FASTA header descriptions share no keyword at all -- see
claude_code_instructions_domain_complementarity_v3.md section 3.3 for the
motivating example (MA_4115's misleading "alpha hydrolase" annotation).

Follows the same YAML load -> validate -> frozen dataclass pattern as
analysis/functional_complementarity_rules.py: a configured but
missing/invalid file is a hard ``ConfigError`` (this file is explicitly
opted into by the user, unlike the UniProt bulk export, which silently
disables itself -- see annotation/uniprot_bulk.py).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from annotation.uniprot_bulk import UniProtDomainInfo
from core.exceptions import ConfigError


@dataclass(frozen=True, slots=True)
class DomainFamilyCategory:
    """One domain-family category's Pfam/InterPro/SUPFAM accession set."""

    pfam: frozenset[str] = frozenset()
    interpro: frozenset[str] = frozenset()
    supfam: frozenset[str] = frozenset()


@dataclass(frozen=True, slots=True)
class DomainFamilyMap:
    """A versioned collection of domain-family categories plus complementarity rules."""

    version: str
    categories: dict[str, DomainFamilyCategory]
    category_rules: tuple[tuple[str, str, str], ...]  # (left_category, right_category, note)

    def categories_for(self, info: "UniProtDomainInfo | None") -> set[str]:
        """Return every category ``info`` belongs to (empty set if ``info`` is None)."""
        if info is None:
            return set()

        matched: set[str] = set()
        for category_id, category in self.categories.items():
            if (
                (info.pfam & category.pfam)
                or (info.interpro & category.interpro)
                or (info.supfam & category.supfam)
            ):
                matched.add(category_id)
        return matched

    def find_category_match(
        self, query_categories: set[str], candidate_categories: set[str]
    ) -> tuple[str, str, str] | None:
        """Return the first ``category_rules`` entry that matches in either orientation."""
        for left, right, note in self.category_rules:
            forward = left in query_categories and right in candidate_categories
            backward = right in query_categories and left in candidate_categories
            if forward or backward:
                return (left, right, note)
        return None


def load_domain_family_map(path: str | Path) -> DomainFamilyMap:
    """Load a domain family map from ``path``."""
    resolved_path = Path(path)
    if not resolved_path.exists():
        raise ConfigError(f"domain family map file was not found: {resolved_path}")

    try:
        raw = yaml.safe_load(resolved_path.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ConfigError(
            f"domain family map is not valid YAML: {resolved_path} ({exc})"
        ) from exc

    if not isinstance(raw, dict):
        raise ConfigError(f"domain family map must contain a mapping: {resolved_path}")

    version = str(raw.get("version", "unknown"))

    raw_categories = raw.get("categories", {})
    if not isinstance(raw_categories, dict):
        raise ConfigError(f"'categories' in {resolved_path} must be a mapping.")

    categories: dict[str, DomainFamilyCategory] = {}
    for category_id, raw_category in raw_categories.items():
        if not isinstance(raw_category, dict):
            raise ConfigError(f"'categories.{category_id}' in {resolved_path} must be a mapping.")
        categories[str(category_id)] = DomainFamilyCategory(
            pfam=frozenset(_string_list(raw_category.get("pfam"), f"categories.{category_id}.pfam", resolved_path)),
            interpro=frozenset(
                _string_list(raw_category.get("interpro"), f"categories.{category_id}.interpro", resolved_path)
            ),
            supfam=frozenset(
                _string_list(raw_category.get("supfam"), f"categories.{category_id}.supfam", resolved_path)
            ),
        )

    raw_category_rules = raw.get("category_rules", [])
    if not isinstance(raw_category_rules, list):
        raise ConfigError(f"'category_rules' in {resolved_path} must be a list.")

    category_rules: list[tuple[str, str, str]] = []
    for index, raw_rule in enumerate(raw_category_rules):
        if not isinstance(raw_rule, dict):
            raise ConfigError(f"'category_rules[{index}]' in {resolved_path} must be a mapping.")
        left = str(raw_rule.get("left") or "").strip()
        right = str(raw_rule.get("right") or "").strip()
        if not left or not right:
            raise ConfigError(
                f"'category_rules[{index}]' in {resolved_path} needs both 'left' and 'right'."
            )
        if left not in categories:
            raise ConfigError(
                f"'category_rules[{index}].left' ({left!r}) in {resolved_path} is not a defined category."
            )
        if right not in categories:
            raise ConfigError(
                f"'category_rules[{index}].right' ({right!r}) in {resolved_path} is not a defined category."
            )
        category_rules.append((left, right, str(raw_rule.get("note", ""))))

    return DomainFamilyMap(
        version=version,
        categories=categories,
        category_rules=tuple(category_rules),
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
        result.append(item)
    return result


__all__: tuple[str, ...] = (
    "DomainFamilyCategory",
    "DomainFamilyMap",
    "load_domain_family_map",
)
