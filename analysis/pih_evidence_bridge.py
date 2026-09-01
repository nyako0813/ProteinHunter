"""Optional bridge to ProteinInteractionHunter's evidence output.

ProteinInteractionHunter (PIH) is a separate, independent project and is
never imported here -- see docs/scoring_engine_v2.md and the companion
design specification for why. This module only reads PIH's own
machine-readable output file (`candidate_evidence_bundle.jsonl`, produced
by `protein-interaction-hunter generate-candidates`) as plain JSON lines,
the same way any other external tool's report could be read.

PIH's MVP-1K "integrated scoring" (see its docs/scoring.md) groups
evidence into five categories: genomic_context, functional_annotation,
cellular_compatibility, evolutionary, and direct_interaction. v5 already
computes its own genomic_context (GFF gene distance) and
functional_annotation (BLAST source overlap + domain/description
keywords) independently. Folding PIH's versions of those same two
categories into v5's score would double count the same kind of signal
computed two different ways, so this bridge deliberately only imports the
three categories v5 has no equivalent for: cellular_compatibility
(localization/topology), evolutionary (orthology / phylogenetic profile),
and direct_interaction (gene fusion / known interaction databases). Their
category names are prefixed with "pih_" in v5's scoring breakdown so they
are always visually and structurally distinct from v5's own categories.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

#: PIH category_name values that have no v5 equivalent and are safe to
#: import without double counting. PIH's own "genomic_context" and
#: "functional_annotation" categories are intentionally excluded.
BRIDGED_PIH_CATEGORIES: tuple[str, ...] = (
    "cellular_compatibility",
    "evolutionary",
    "direct_interaction",
)

#: Weight budget (in v5's output_scale points) each bridged PIH category
#: contributes, mirroring analysis/interaction_scoring.py::V2_COMPONENT_WEIGHTS.
#: These are provisional defaults, not a calibrated model -- see
#: docs/scoring_engine_v2.md.
PIH_CATEGORY_WEIGHTS: dict[str, float] = {
    "cellular_compatibility": 5.0,
    "evolutionary": 10.0,
    "direct_interaction": 20.0,
}

#: analysis/scoring_engine_config.py::DEFAULT_CATEGORY_CAPS extension: the
#: cap each bridged category may contribute to v5's total score, active
#: only for pairs where PIH actually produced that category's evidence.
BRIDGED_PIH_CATEGORY_CAPS: dict[str, float] = {
    "pih_cellular_compatibility": 5.0,
    "pih_evolutionary": 10.0,
    "pih_direct_interaction": 20.0,
}


@dataclass(slots=True, frozen=True)
class PihCategoryEvidence:
    """One bridged PIH category score for one query/candidate pair."""

    category_name: str  # PIH's own name, e.g. "direct_interaction"
    normalized_score: float  # 0.0-1.0, already PIH's own within-category normalization
    available_weight: float  # PIH's own available_weight; 0 means "not active for this pair"


@dataclass(slots=True, frozen=True)
class PihEvidenceBundle:
    """A parsed, queryable index of one PIH candidate_evidence_bundle.jsonl."""

    pairs: dict[tuple[str, str], dict[str, PihCategoryEvidence]]
    warnings: tuple[str, ...]

    def lookup(self, query_keys: list[str], candidate_keys: list[str]) -> dict[str, PihCategoryEvidence]:
        """Return bridged category evidence for the first matching ID pair.

        Both PIH and v5 are given several candidate spellings of the same
        protein (protein_id, old_locus_tag, a version-stripped protein_id)
        because the two tools were not designed to share an identifier
        convention. The first (query_key, candidate_key) combination found
        in the bundle wins; returns {} when nothing matches.
        """
        for query_key in query_keys:
            if not query_key:
                continue
            for candidate_key in candidate_keys:
                if not candidate_key:
                    continue
                match = self.pairs.get((query_key, candidate_key))
                if match is not None:
                    return match
        return {}


def load_pih_evidence_bundle(path: Path) -> PihEvidenceBundle:
    """Parse a PIH candidate_evidence_bundle.jsonl file.

    Malformed lines are skipped with a warning rather than aborting the
    whole run -- this bridge is optional, best-effort evidence, exactly
    like any other external adapter in the design specification's
    boundary rules (an unavailable or broken external source must not
    stop a local run).
    """
    warnings: list[str] = []
    pairs: dict[tuple[str, str], dict[str, PihCategoryEvidence]] = {}

    if not path.exists():
        return PihEvidenceBundle(pairs={}, warnings=(f"PIH evidence bundle not found: {path}",))

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                warnings.append(f"{path}:{line_number}: not valid JSON ({exc}); skipped")
                continue
            if not isinstance(record, dict):
                warnings.append(f"{path}:{line_number}: not a JSON object; skipped")
                continue

            query_id = record.get("query_id")
            candidate_id = record.get("candidate_id")
            if not isinstance(query_id, str) or not isinstance(candidate_id, str):
                warnings.append(f"{path}:{line_number}: missing query_id/candidate_id; skipped")
                continue

            categories = _extract_bridged_categories(record, path, line_number, warnings)
            if categories:
                pairs[(query_id, candidate_id)] = categories

    return PihEvidenceBundle(pairs=pairs, warnings=tuple(warnings))


def _extract_bridged_categories(
    record: dict, path: Path, line_number: int, warnings: list[str]
) -> dict[str, PihCategoryEvidence]:
    integrated_scoring = record.get("integrated_scoring") or record.get("score")
    if not isinstance(integrated_scoring, dict):
        return {}

    raw_categories = integrated_scoring.get("category_scores")
    if not isinstance(raw_categories, list):
        return {}

    result: dict[str, PihCategoryEvidence] = {}
    for raw_category in raw_categories:
        if not isinstance(raw_category, dict):
            continue
        category_name = raw_category.get("category_name")
        if category_name not in BRIDGED_PIH_CATEGORIES:
            continue
        try:
            normalized_score = float(raw_category.get("normalized_score", 0.0))
            available_weight = float(raw_category.get("available_weight", 0.0))
        except (TypeError, ValueError):
            warnings.append(
                f"{path}:{line_number}: category '{category_name}' has a non-numeric score; skipped"
            )
            continue
        if available_weight <= 0.0:
            continue  # PIH itself had no active evidence for this category on this pair
        result[category_name] = PihCategoryEvidence(
            category_name=category_name,
            normalized_score=max(0.0, min(1.0, normalized_score)),
            available_weight=available_weight,
        )
    return result


def without_version(protein_id: str) -> str:
    """Strip a trailing '.N' version suffix, mirroring identifier handling
    already used elsewhere in this project (see e.g. _without_version in
    analysis/interaction_scoring.py)."""
    return re.sub(r"\.\d+$", "", protein_id)


__all__: tuple[str, ...] = (
    "BRIDGED_PIH_CATEGORIES",
    "BRIDGED_PIH_CATEGORY_CAPS",
    "PIH_CATEGORY_WEIGHTS",
    "PihCategoryEvidence",
    "PihEvidenceBundle",
    "load_pih_evidence_bundle",
    "without_version",
)
