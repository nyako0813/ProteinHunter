"""Synthetic ProteinInteractionHunter (PIH) candidate_evidence_bundle.jsonl records.

Shaped after PIH's own schema (schemas/candidate_evidence_bundle.schema.json in
the PIH repository) and its writer (application/pipeline.py), not after what
this project assumed PIH looked like: `integrated_scoring` is a LIST holding at
most one IntegratedScore, and each ScoreCategory carries the schema's required
`raw_weighted_sum` / `configured_cap` fields next to `normalized_score`. The
records here were validated against the PIH schema (Draft 2020-12) when they
were written; the PIH repository is not a dependency of this project, so that
check is not part of the test suite.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

#: PIH's five ScoreCategory names (application/scoring.py::COMPONENT_CATEGORIES)
#: mapped to the component names that feed them, used only to build a plausible
#: component_scores list.
_PIH_CATEGORY_COMPONENTS = {
    "genomic_context": "genome_context",
    "functional_annotation": "domain_pair",
    "cellular_compatibility": "localization",
    "evolutionary": "orthology",
    "direct_interaction": "fusion",
}


def score_category(
    category_name: str,
    normalized_score: float,
    available_weight: float = 1.0,
    *,
    configured_cap: float = 2.0,
) -> dict[str, Any]:
    """One PIH ScoreCategory (all schema-required fields)."""
    return {
        "category_name": category_name,
        "raw_weighted_sum": round(normalized_score * available_weight, 8),
        "available_weight": available_weight,
        "normalized_score": normalized_score,
        "configured_cap": configured_cap,
    }


def integrated_score(query_id: str, candidate_id: str, categories: list[dict[str, Any]]) -> dict[str, Any]:
    """One PIH IntegratedScore (all schema-required fields) wrapping ``categories``."""
    active = [category for category in categories if category["available_weight"] > 0]
    raw_sum = sum(category["raw_weighted_sum"] for category in active)
    available = sum(category["available_weight"] for category in active)
    return {
        "query_protein_id": query_id,
        "candidate_protein_id": candidate_id,
        "raw_weighted_sum": raw_sum,
        "available_weight": available,
        "evidence_category_count": len(active),
        "evidence_component_count": len(active),
        "positive_component_count": len(active),
        "neutral_component_count": 0,
        "negative_component_count": 0,
        "sufficient_evidence": len(active) >= 2 and available >= 1.0,
        "calculation_rule_version": "mvp1k-integrated-scoring-v1",
        "category_scores": categories,
        "component_scores": [
            {
                "component_name": _PIH_CATEGORY_COMPONENTS.get(category["category_name"], "genome_context"),
                "category_name": category["category_name"],
                "evidence_status": "available",
                "configured_weight": category["available_weight"],
                "direction": "positive" if category["normalized_score"] > 0 else "neutral",
                "applied": True,
                "effective_weight": category["available_weight"],
                "normalized_value": category["normalized_score"],
                "weighted_contribution": category["raw_weighted_sum"],
            }
            for category in active
        ],
        "normalized_score": max(0.0, min(1.0, raw_sum / available)) if available else None,
        "output_score": round(100.0 * max(0.0, min(1.0, raw_sum / available)), 8) if available else None,
        "status": "available" if active else "not_run",
    }


def bundle_record(
    query_id: str,
    candidate_id: str,
    integrated_scoring: list[dict[str, Any]] | None,
    *,
    run_id: str = "synthetic-run",
) -> dict[str, Any]:
    """One CandidateEvidenceBundle line. ``integrated_scoring`` is PIH's list (empty = not scored)."""
    return {
        "run_id": run_id,
        "schema_version": "1.0",
        "query_id": query_id,
        "candidate_id": candidate_id,
        "candidate_disposition": "included",
        "predicted_relationship_type": "insufficient_evidence",
        "engine_statuses": {"orthology": "available", "localization": "available"},
        "integrated_scoring": [] if integrated_scoring is None else integrated_scoring,
        # PIH's separate per-candidate score summary: a different structure with
        # NO category breakdown. It must never be mistaken for integrated_scoring.
        "score": {"gene_context_score": 0.8, "total_ranking_score": 0.5},
    }


def write_bundle(path: Path, records: list[dict[str, Any]]) -> Path:
    """Write ``records`` as PIH writes them: compact, key-sorted JSON, one per line."""
    path.write_text(
        "".join(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n" for record in records),
        encoding="utf-8",
    )
    return path


def write_pih_bundle(
    path: Path, *, query_id: str, candidate_id: str, category_scores: list[dict[str, Any]]
) -> Path:
    """A one-record bundle whose pair was scored by PIH, from ``{category_name, normalized_score, available_weight}`` dicts."""
    categories = [
        score_category(item["category_name"], item["normalized_score"], item["available_weight"])
        for item in category_scores
    ]
    return write_bundle(
        path, [bundle_record(query_id, candidate_id, [integrated_score(query_id, candidate_id, categories)])]
    )
