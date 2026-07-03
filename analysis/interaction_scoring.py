"""Lightweight interaction candidate ranking helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

from core.fasta import read_fasta_as_components
from core.models import ProteinRecord


@dataclass(frozen=True)
class InteractionScoringResult:
    """Interaction query rows and per-source pair rows for Excel output."""

    query_rows: list[dict[str, Any]]
    source_rows: dict[str, list[dict[str, Any]]]
    warnings: list[str]


CANDIDATE_SOURCE_MAP: dict[str, tuple[str, str]] = {
    "candidates": ("Candidates", "positive_only_records"),
    "candidates_relaxed": ("Candidates_relaxed", "candidates_relaxed_records"),
    "positive_all_sources": ("Positive_all_sources", "positive_all_sources_records"),
    "negative_unmatched": ("Negative_unmatched", "negative_unmatched_records"),
    "no_hit": ("No_hit", "no_hit_records"),
    "negative_hit": ("Negative_hit", "negative_hit_records"),
    "negative_strong_hit": ("Negative_strong_hit", "negative_strong_hit_records"),
    "negative_medium_hit": ("Negative_medium_hit", "negative_medium_hit_records"),
    "negative_weak_hit": ("Negative_weak_hit", "negative_weak_hit_records"),
}

INTERACTION_SHEET_DESCRIPTIONS: dict[str, tuple[str, str, str]] = {
    "Interaction_query": (
        "user-provided query proteins resolved by protein_id, old_locus_tag, and/or sequence",
        "confirms which target proteins are used as interaction queries",
        "check query resolution before interpreting interaction rankings",
    ),
    "Interaction_Candidates": (
        "interaction ranking against strict Candidates",
        "ranks possible partners among high-confidence candidate proteins",
        "first-pass interaction partner review",
    ),
    "Interaction_Candidates_relaxed": (
        "interaction ranking against Candidates_relaxed",
        "ranks possible partners while avoiding over-filtering by medium/weak negative hits",
        "review additional candidates that may have been missed by strict filtering",
    ),
    "Interaction_Positive_all_sources": (
        "interaction ranking against Positive_all_sources",
        "ranks possible partners among broadly positive-conserved candidates",
        "inspect interaction candidates from broadly conserved positive-source hits",
    ),
    "Interaction_Negative_unmatched": (
        "interaction ranking against Negative_unmatched",
        "ranks possible partners among proteins without negative hits",
        "review broad negative-unmatched interaction candidates",
    ),
    "Interaction_No_hit": (
        "interaction ranking against No_hit proteins",
        "ranks possible partners among poorly conserved or target-specific proteins",
        "explore novel or lineage-specific interaction candidates",
    ),
    "Interaction_Negative_hit": (
        "interaction ranking against Negative_hit proteins",
        "ranks possible partners among proteins with any negative hit",
        "inspect lower-priority or cautionary interaction candidates",
    ),
    "Interaction_Negative_strong_hit": (
        "interaction ranking against Negative_strong_hit proteins",
        "ranks possible partners among proteins with strong negative hits",
        "generally lower-priority review",
    ),
    "Interaction_Negative_medium_hit": (
        "interaction ranking against Negative_medium_hit proteins",
        "ranks possible partners among ambiguous homolog candidates",
        "review only when medium negative hits should not be automatically excluded",
    ),
    "Interaction_Negative_weak_hit": (
        "interaction ranking against Negative_weak_hit proteins",
        "ranks possible partners among weak negative-hit candidates",
        "review as cautionary candidates that may represent distant homologs or shared domains",
    ),
}

INTERACTION_QUERY_COLUMNS: tuple[str, ...] = (
    "query_id",
    "input_protein_id",
    "input_old_locus_tag",
    "resolved_protein_id",
    "resolved_old_locus_tag",
    "sequence_length",
    "resolution_status",
    "description",
    "notes",
)

INTERACTION_PAIR_COLUMNS: tuple[str, ...] = (
    "query_id",
    "query_protein_id",
    "query_old_locus_tag",
    "candidate_rank",
    "candidate_protein_id",
    "candidate_old_locus_tag",
    "candidate_source",
    "candidate_description",
    "interaction_priority_score",
    "interaction_score_reasons",
    "candidate_priority_score",
    "same_gene_neighborhood_score",
    "distance_bp",
    "co_occurrence_score",
    "domain_complementarity_score",
    "alphafold_readiness_score",
    "pair_total_length",
    "alphafold_recommended",
)

SEQUENCE_COLUMNS: tuple[str, ...] = (
    "query_sequence",
    "candidate_sequence",
)

CANDIDATE_PRIORITY_BASE: dict[str, float] = {
    "Candidates": 30.0,
    "Candidates_relaxed": 25.0,
    "Positive_all_sources": 30.0,
    "Negative_unmatched": 15.0,
    "No_hit": 20.0,
    "Negative_weak_hit": 10.0,
    "Negative_medium_hit": 5.0,
    "Negative_strong_hit": 0.0,
    "Negative_hit": 0.0,
}


def run_interaction_scoring(
    config: Any,
    blast_classification: Any,
) -> InteractionScoringResult | None:
    """Rank candidate proteins for configured query proteins."""
    scoring_config = config.interaction_scoring
    if not scoring_config.enabled:
        return None

    warnings: list[str] = []
    queries = _load_query_specs(scoring_config, warnings)
    if not queries:
        warnings.append("interaction_scoring enabled but no query protein was provided")
        return InteractionScoringResult(query_rows=[], source_rows={}, warnings=warnings)

    resolved_queries = [
        _resolve_query(query, index, blast_classification.all_records)
        for index, query in enumerate(queries, start=1)
    ]
    query_rows = [_query_row(query) for query in resolved_queries]
    source_rows: dict[str, list[dict[str, Any]]] = {}

    for source_key, enabled in scoring_config.candidate_sources.items():
        if not enabled:
            continue
        source_info = CANDIDATE_SOURCE_MAP.get(source_key)
        if source_info is None:
            warnings.append(f"unsupported interaction candidate source skipped: {source_key}")
            continue

        source_label, attr_name = source_info
        candidate_records = getattr(blast_classification, attr_name, None)
        if not candidate_records:
            warnings.append(f"interaction candidate source has no records: {source_key}")
            continue

        rows = _rank_source_candidates(
            resolved_queries=resolved_queries,
            candidate_records=candidate_records,
            candidate_source=source_label,
            scoring_config=scoring_config,
        )
        if rows:
            source_rows[f"Interaction_{source_label}"] = rows
        else:
            warnings.append(f"interaction candidate source produced no pairs: {source_key}")

    return InteractionScoringResult(
        query_rows=query_rows,
        source_rows=source_rows,
        warnings=warnings,
    )


def interaction_pair_columns(include_sequences: bool) -> tuple[str, ...]:
    """Return interaction pair columns for Excel output."""
    if include_sequences:
        return (*INTERACTION_PAIR_COLUMNS, *SEQUENCE_COLUMNS)

    return INTERACTION_PAIR_COLUMNS


def interaction_index_rows(
    sheet_names: list[str],
) -> list[tuple[str, str, str, str]]:
    """Return Index rows for actually created interaction sheets."""
    rows: list[tuple[str, str, str, str]] = []
    for sheet_name in sheet_names:
        descriptions = INTERACTION_SHEET_DESCRIPTIONS.get(
            sheet_name,
            (
                f"interaction ranking for {sheet_name}",
                "lightweight interaction prioritization output",
                "review candidate pairs before manual follow-up",
            ),
        )
        rows.append((sheet_name, *descriptions))

    return rows


def _load_query_specs(scoring_config: Any, warnings: list[str]) -> list[dict[str, str]]:
    query_specs: list[dict[str, str]] = []
    for query in scoring_config.query_proteins:
        query_specs.append(
            {
                "protein_id": query.protein_id,
                "old_locus_tag": query.old_locus_tag,
                "sequence": query.sequence,
                "description": "",
                "source": "config",
            }
        )

    if scoring_config.query_fasta is not None:
        try:
            ids, descriptions, sequences = read_fasta_as_components(
                Path(scoring_config.query_fasta)
            )
        except Exception as exc:
            warnings.append(f"query_fasta could not be read: {exc}")
        else:
            for protein_id in ids:
                query_specs.append(
                    {
                        "protein_id": protein_id,
                        "old_locus_tag": "",
                        "sequence": sequences.get(protein_id, ""),
                        "description": descriptions.get(protein_id, ""),
                        "source": "query_fasta",
                    }
                )

    return [
        query
        for query in query_specs
        if query["protein_id"] or query["old_locus_tag"] or query["sequence"]
    ]


def _resolve_query(
    query: dict[str, str],
    index: int,
    records: dict[str, ProteinRecord],
) -> dict[str, Any]:
    matched_record = _find_query_record(query, records)
    sequence = query["sequence"]
    if not sequence and matched_record is not None:
        sequence = matched_record.sequence

    query_id = query["protein_id"] or query["old_locus_tag"] or f"query_{index}"
    if matched_record is not None:
        resolved_protein_id = matched_record.protein_id
        resolved_old_locus_tag = matched_record.old_locus_tag or ""
        description = matched_record.description
        status = "resolved"
        notes = "resolved from target records"
    elif sequence:
        resolved_protein_id = query["protein_id"]
        resolved_old_locus_tag = query["old_locus_tag"]
        description = query.get("description", "")
        status = "resolved"
        notes = "using explicit sequence"
    else:
        resolved_protein_id = ""
        resolved_old_locus_tag = ""
        description = query.get("description", "")
        status = "unresolved"
        notes = "no matching target record or sequence"

    return {
        "query_id": query_id,
        "input_protein_id": query["protein_id"],
        "input_old_locus_tag": query["old_locus_tag"],
        "resolved_protein_id": resolved_protein_id,
        "resolved_old_locus_tag": resolved_old_locus_tag,
        "sequence": sequence,
        "sequence_length": len(sequence) if sequence else None,
        "resolution_status": status,
        "description": description,
        "notes": notes,
        "record": matched_record,
    }


def _find_query_record(
    query: dict[str, str],
    records: dict[str, ProteinRecord],
) -> ProteinRecord | None:
    protein_id = query["protein_id"]
    if protein_id and protein_id in records:
        return records[protein_id]

    old_locus_tag = query["old_locus_tag"]
    if old_locus_tag:
        for record in records.values():
            if record.old_locus_tag == old_locus_tag:
                return record

    return None


def _query_row(query: dict[str, Any]) -> dict[str, Any]:
    return {
        "query_id": query["query_id"],
        "input_protein_id": query["input_protein_id"],
        "input_old_locus_tag": query["input_old_locus_tag"],
        "resolved_protein_id": query["resolved_protein_id"],
        "resolved_old_locus_tag": query["resolved_old_locus_tag"],
        "sequence_length": query["sequence_length"],
        "resolution_status": query["resolution_status"],
        "description": query["description"],
        "notes": query["notes"],
    }


def _rank_source_candidates(
    resolved_queries: list[dict[str, Any]],
    candidate_records: dict[str, ProteinRecord],
    candidate_source: str,
    scoring_config: Any,
) -> list[dict[str, Any]]:
    all_rows: list[dict[str, Any]] = []
    for query in resolved_queries:
        if query["resolution_status"] != "resolved":
            continue
        query_rows: list[dict[str, Any]] = []
        for candidate in candidate_records.values():
            if _is_self_pair(query, candidate):
                continue
            query_rows.append(
                _score_pair(
                    query=query,
                    candidate=candidate,
                    candidate_source=candidate_source,
                    scoring_config=scoring_config,
                )
            )

        query_rows.sort(
            key=lambda row: (
                -float(row["interaction_priority_score"]),
                not bool(row["alphafold_recommended"]),
                str(row["candidate_protein_id"]),
            )
        )
        for rank, row in enumerate(
            query_rows[: scoring_config.max_candidates_per_query],
            start=1,
        ):
            row["candidate_rank"] = rank
            all_rows.append(row)

    all_rows.sort(
        key=lambda row: (
            str(row["query_id"]),
            int(row["candidate_rank"]),
            str(row["candidate_protein_id"]),
        )
    )
    return all_rows


def _score_pair(
    query: dict[str, Any],
    candidate: ProteinRecord,
    candidate_source: str,
    scoring_config: Any,
) -> dict[str, Any]:
    weights = scoring_config.scoring_weights
    reasons: list[str] = [f"candidate source: {candidate_source}"]
    candidate_priority_score = _candidate_priority_score(
        candidate_source,
        weights.candidate_priority,
    )
    same_gene_neighborhood_score = 0.0
    distance_bp = None
    reasons.append("genomic distance unavailable")

    co_occurrence_score = _co_occurrence_score(query["record"], candidate, weights)
    if co_occurrence_score:
        reasons.append("similar source pattern")
    else:
        reasons.append("no source pattern evidence")

    domain_complementarity_score = _domain_complementarity_score(
        query,
        candidate,
        weights.domain_complementarity,
        reasons,
    )
    alphafold_readiness_score, pair_total_length, alphafold_recommended = (
        _alphafold_readiness(query, candidate, scoring_config)
    )
    if alphafold_recommended:
        reasons.append("compatible for manual AlphaFold")
    else:
        reasons.append("missing sequence or length too large for manual AlphaFold")

    total_score = (
        candidate_priority_score
        + same_gene_neighborhood_score
        + co_occurrence_score
        + domain_complementarity_score
        + alphafold_readiness_score
    )
    row = {
        "query_id": query["query_id"],
        "query_protein_id": query["resolved_protein_id"],
        "query_old_locus_tag": query["resolved_old_locus_tag"],
        "candidate_rank": 0,
        "candidate_protein_id": candidate.protein_id,
        "candidate_old_locus_tag": candidate.old_locus_tag or "",
        "candidate_source": candidate_source,
        "candidate_description": candidate.description,
        "interaction_priority_score": round(total_score, 3),
        "interaction_score_reasons": "; ".join(reasons),
        "candidate_priority_score": round(candidate_priority_score, 3),
        "same_gene_neighborhood_score": same_gene_neighborhood_score,
        "distance_bp": distance_bp,
        "co_occurrence_score": round(co_occurrence_score, 3),
        "domain_complementarity_score": round(domain_complementarity_score, 3),
        "alphafold_readiness_score": round(alphafold_readiness_score, 3),
        "pair_total_length": pair_total_length,
        "alphafold_recommended": alphafold_recommended,
    }
    if scoring_config.include_sequences_in_excel:
        row["query_sequence"] = query["sequence"]
        row["candidate_sequence"] = candidate.sequence

    return row


def _candidate_priority_score(candidate_source: str, weight: float) -> float:
    base = CANDIDATE_PRIORITY_BASE.get(candidate_source, 10.0)
    return base / 30.0 * weight


def _co_occurrence_score(
    query_record: ProteinRecord | None,
    candidate: ProteinRecord,
    weights: Any,
) -> float:
    if query_record is None:
        return 0.0

    query_sources = set(query_record.positive_sources_hit)
    candidate_sources = set(candidate.positive_sources_hit)
    if query_sources or candidate_sources:
        union = query_sources | candidate_sources
        if union:
            similarity = len(query_sources & candidate_sources) / len(union)
            return similarity * weights.co_occurrence

    if not query_record.negative_hits and not candidate.negative_hits:
        return weights.co_occurrence * 0.5

    return 0.0


def _domain_complementarity_score(
    query: dict[str, Any],
    candidate: ProteinRecord,
    weight: float,
    reasons: list[str],
) -> float:
    query_text = _record_text(query["record"], query["description"])
    candidate_text = _record_text(candidate, candidate.description)
    combined = f"{query_text} {candidate_text}".lower()

    keyword_pairs = (
        ("enzyme", "carrier"),
        ("atpase", "transfer"),
        ("radical", "sulfur"),
        ("kinase", "enzyme"),
        ("ligase", "enzyme"),
        ("transferase", "enzyme"),
    )
    for left, right in keyword_pairs:
        if left in combined and right in combined:
            reasons.append("generic domain/description complementarity")
            return weight

    shared_words = _informative_words(query_text) & _informative_words(candidate_text)
    if shared_words:
        reasons.append("shared pathway-like description words")
        return weight * 0.5

    reasons.append("no domain evidence available")
    return 0.0


def _alphafold_readiness(
    query: dict[str, Any],
    candidate: ProteinRecord,
    scoring_config: Any,
) -> tuple[float, int | None, bool]:
    if not query["sequence"] or not candidate.sequence:
        return 0.0, None, False

    pair_total_length = len(query["sequence"]) + len(candidate.sequence)
    recommended = pair_total_length <= scoring_config.alphafold.max_pair_total_length
    score = scoring_config.scoring_weights.alphafold_readiness if recommended else 0.0
    return score, pair_total_length, recommended


def _is_self_pair(query: dict[str, Any], candidate: ProteinRecord) -> bool:
    return bool(query["resolved_protein_id"]) and (
        query["resolved_protein_id"] == candidate.protein_id
    )


def _record_text(record: ProteinRecord | None, fallback_description: str) -> str:
    if record is None:
        return fallback_description

    domains = " ".join(
        f"{domain.name} {domain.description}" for domain in record.domains
    )
    return f"{record.description} {domains}"


def _informative_words(text: str) -> set[str]:
    stop_words = {
        "protein",
        "hypothetical",
        "putative",
        "probable",
        "domain",
        "family",
        "like",
    }
    return {
        word
        for word in re.findall(r"[a-zA-Z][a-zA-Z0-9_-]{3,}", text.lower())
        if word not in stop_words
    }


__all__: tuple[str, ...] = (
    "CANDIDATE_SOURCE_MAP",
    "INTERACTION_PAIR_COLUMNS",
    "INTERACTION_QUERY_COLUMNS",
    "INTERACTION_SHEET_DESCRIPTIONS",
    "InteractionScoringResult",
    "interaction_index_rows",
    "interaction_pair_columns",
    "run_interaction_scoring",
)
