"""Lightweight interaction candidate ranking helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Any

from annotation.gff import GffFeatureLocation, load_gff_feature_map
from core.fasta import read_fasta_as_components
from core.models import ProteinRecord


@dataclass(frozen=True)
class InteractionScoringResult:
    """Interaction query rows and per-source pair rows for Excel output."""

    query_rows: list[dict[str, Any]]
    source_rows: dict[str, list[dict[str, Any]]]
    neighborhood_rows: list[dict[str, Any]]
    warnings: list[str]


CANDIDATE_SOURCE_MAP: dict[str, tuple[str, str, str]] = {
    "candidates": ("Candidates", "positive_only_records", "Interaction_Candidates"),
    "candidates_relaxed": (
        "Candidates_relaxed",
        "candidates_relaxed_records",
        "Interaction_Candidates_relaxed",
    ),
    "positive_all_sources": (
        "Positive_all_sources",
        "positive_all_sources_records",
        "Interaction_Positive_all",
    ),
    "negative_unmatched": (
        "Negative_unmatched",
        "negative_unmatched_records",
        "Interaction_Neg_unmatched",
    ),
    "no_hit": ("No_hit", "no_hit_records", "Interaction_No_hit"),
    "negative_hit": ("Negative_hit", "negative_hit_records", "Interaction_Neg_hit"),
    "negative_strong_hit": (
        "Negative_strong_hit",
        "negative_strong_hit_records",
        "Interaction_Neg_strong",
    ),
    "negative_medium_hit": (
        "Negative_medium_hit",
        "negative_medium_hit_records",
        "Interaction_Neg_medium",
    ),
    "negative_weak_hit": (
        "Negative_weak_hit",
        "negative_weak_hit_records",
        "Interaction_Neg_weak",
    ),
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
    "Interaction_Positive_all": (
        "interaction ranking against Positive_all_sources",
        "ranks possible partners among broadly positive-conserved candidates",
        "inspect interaction candidates from broadly conserved positive-source hits",
    ),
    "Interaction_Neg_unmatched": (
        "interaction ranking against Negative_unmatched",
        "ranks possible partners among proteins without negative hits",
        "review broad negative-unmatched interaction candidates",
    ),
    "Interaction_No_hit": (
        "interaction ranking against No_hit proteins",
        "ranks possible partners among poorly conserved or target-specific proteins",
        "explore novel or lineage-specific interaction candidates",
    ),
    "Interaction_Neg_hit": (
        "interaction ranking against Negative_hit proteins",
        "ranks possible partners among proteins with any negative hit",
        "inspect lower-priority or cautionary interaction candidates",
    ),
    "Interaction_Neg_strong": (
        "interaction ranking against Negative_strong_hit proteins",
        "ranks possible partners among proteins with strong negative hits",
        "generally lower-priority review",
    ),
    "Interaction_Neg_medium": (
        "interaction ranking against Negative_medium_hit proteins",
        "ranks possible partners among ambiguous homolog candidates",
        "review only when medium negative hits should not be automatically excluded",
    ),
    "Interaction_Neg_weak": (
        "interaction ranking against Negative_weak_hit proteins",
        "ranks possible partners among weak negative-hit candidates",
        "review as cautionary candidates that may represent distant homologs or shared domains",
    ),
    "Interaction_Neighborhood": (
        "nearby interaction candidates around each query protein",
        "summarizes same-contig candidate neighborhoods using GFF coordinates",
        "review local genomic context around interaction queries",
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
    "same_contig",
    "query_start",
    "query_end",
    "query_strand",
    "candidate_start",
    "candidate_end",
    "candidate_strand",
    "distance_bp",
    "strand_relation",
    "same_gene_neighborhood_score",
    "interaction_priority_score",
    "distance_independent_score",
    "distance_independent_rank",
    "priority_group",
    "interaction_score_reasons",
    "candidate_priority_score",
    "co_occurrence_score",
    "domain_complementarity_score",
    "alphafold_readiness_score",
    "pair_total_length",
    "alphafold_recommended",
)

SEQUENCE_COLUMNS: tuple[str, ...] = ("query_sequence", "candidate_sequence")

INTERACTION_NEIGHBORHOOD_SHEET = "Interaction_Neighborhood"

INTERACTION_NEIGHBORHOOD_COLUMNS: tuple[str, ...] = (
    "query_id",
    "query_protein_id",
    "query_old_locus_tag",
    "query_description",
    "query_contig",
    "query_start",
    "query_end",
    "query_strand",
    "candidate_rank_by_distance",
    "candidate_protein_id",
    "candidate_old_locus_tag",
    "candidate_description",
    "candidate_source",
    "candidate_contig",
    "candidate_start",
    "candidate_end",
    "candidate_strand",
    "distance_bp",
    "strand_relation",
    "neighborhood_band",
    "same_gene_neighborhood_score",
    "interaction_priority_score",
    "domain_complementarity_score",
    "candidate_priority_score",
    "co_occurrence_score",
    "alphafold_recommended",
    "interaction_score_reasons",
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

DESCRIPTION_STOPWORDS: set[str] = {
    "protein", "family", "domain", "containing", "domain-containing",
    "hypothetical", "putative", "predicted", "probable", "possible",
    "uncharacterized", "conserved", "archaeal", "bacterial", "homolog",
    "homologous", "like", "subunit", "chain", "component", "associated",
    "related", "membrane", "cytosolic", "methanosarcina", "acetivorans",
    "source", "candidate", "enzyme", "factor", "alpha", "beta",
}

MEANINGFUL_KEYWORDS: set[str] = {
    "atp-binding", "atpase", "ligase", "transferase", "methyltransferase",
    "aminotransferase", "nucleotidyltransferase", "hydrolase",
    "oxidoreductase", "dehydrogenase", "reductase", "radical", "sam",
    "radical sam", "sulfur", "sulphur", "thiol", "thioredoxin",
    "ferredoxin", "fe-s", "iron-sulfur", "molybdopterin", "moad",
    "this", "tusa", "dsre", "carrier", "kinase", "synthetase", "synthase",
    "deaminase", "amidase", "peptidase", "helicase", "nuclease", "rna",
    "trna", "rrna", "dna", "cofactor", "flavin", "fad", "nad", "plp",
    "biotin",
}

COMPLEMENTARY_TERM_PAIRS: tuple[tuple[str, str], ...] = (
    ("radical sam", "iron-sulfur"),
    ("radical sam", "ferredoxin"),
    ("radical", "iron-sulfur"),
    ("sulfur", "carrier"),
    ("sulphur", "carrier"),
    ("thiol", "carrier"),
    ("moad", "this"),
    ("moad", "sulfur"),
    ("this", "sulfur"),
    ("tusa", "sulfur"),
    ("dsre", "sulfur"),
    ("atpase", "ligase"),
    ("atp-binding", "ligase"),
    ("nucleotidyltransferase", "rna"),
    ("methyltransferase", "rna"),
    ("aminotransferase", "plp"),
    ("oxidoreductase", "ferredoxin"),
    ("dehydrogenase", "nad"),
    ("hydrolase", "amidase"),
    ("nuclease", "dna"),
    ("helicase", "dna"),
    ("transferase", "cofactor"),
)


def source_sheet_name(source_key: str) -> str:
    """Return the original classification sheet name for a source key."""
    return CANDIDATE_SOURCE_MAP[source_key][0]


def interaction_sheet_name(source_key: str) -> str:
    """Return the final Excel-safe interaction sheet name for a source key."""
    return CANDIDATE_SOURCE_MAP[source_key][2]


def run_interaction_scoring(config: Any, blast_classification: Any) -> InteractionScoringResult | None:
    """Rank candidate proteins for configured query proteins."""
    scoring_config = config.interaction_scoring
    if not scoring_config.enabled:
        return None

    warnings: list[str] = []
    queries = _load_query_specs(scoring_config, warnings)
    if not queries:
        warnings.append("interaction_scoring enabled but no query protein was provided")
        return InteractionScoringResult(query_rows=[], source_rows={}, neighborhood_rows=[], warnings=warnings)

    feature_map = _load_feature_map(config, warnings)
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

        source_label, attr_name, sheet_name = source_info
        candidate_records = getattr(blast_classification, attr_name, None)
        if not candidate_records:
            warnings.append(f"interaction candidate source has no records: {source_key}")
            continue

        rows = _rank_source_candidates(
            resolved_queries=resolved_queries,
            candidate_records=candidate_records,
            candidate_source=source_label,
            scoring_config=scoring_config,
            feature_map=feature_map,
        )
        if rows:
            source_rows[sheet_name] = rows
        else:
            warnings.append(f"interaction candidate source produced no pairs: {source_key}")

    neighborhood_rows = _build_neighborhood_rows(
        resolved_queries=resolved_queries,
        source_rows=source_rows,
        scoring_config=scoring_config,
    )

    return InteractionScoringResult(
        query_rows=query_rows,
        source_rows=source_rows,
        neighborhood_rows=neighborhood_rows,
        warnings=warnings,
    )


def interaction_neighborhood_columns() -> tuple[str, ...]:
    """Return columns for the optional neighborhood summary sheet."""
    return INTERACTION_NEIGHBORHOOD_COLUMNS


def interaction_pair_columns(include_sequences: bool) -> tuple[str, ...]:
    """Return interaction pair columns for Excel output."""
    if include_sequences:
        return (*INTERACTION_PAIR_COLUMNS, *SEQUENCE_COLUMNS)
    return INTERACTION_PAIR_COLUMNS


def interaction_index_rows(sheet_names: list[str]) -> list[tuple[str, str, str, str]]:
    """Return Index rows for actually created interaction sheets."""
    rows: list[tuple[str, str, str, str]] = []
    for sheet_name in sheet_names:
        descriptions = INTERACTION_SHEET_DESCRIPTIONS.get(
            sheet_name,
            (f"interaction ranking for {sheet_name}", "lightweight interaction prioritization output", "review candidate pairs before manual follow-up"),
        )
        rows.append((sheet_name, *descriptions))
    return rows


def _load_feature_map(config: Any, warnings: list[str]) -> dict[str, GffFeatureLocation]:
    gff_path = getattr(getattr(config, "paths", None), "gff_file", None)
    if gff_path is None:
        return {}
    if not Path(gff_path).exists():
        warnings.append(f"interaction GFF coordinate file was not found: {gff_path}")
        return {}
    try:
        return load_gff_feature_map(gff_path)
    except Exception as exc:
        warnings.append(f"interaction GFF coordinate map could not be loaded: {exc}")
        return {}


def _load_query_specs(scoring_config: Any, warnings: list[str]) -> list[dict[str, str]]:
    query_specs: list[dict[str, str]] = []
    for query in scoring_config.query_proteins:
        query_specs.append({"protein_id": query.protein_id, "old_locus_tag": query.old_locus_tag, "sequence": query.sequence, "description": "", "source": "config"})

    if scoring_config.query_fasta is not None:
        try:
            ids, descriptions, sequences = read_fasta_as_components(Path(scoring_config.query_fasta))
        except Exception as exc:
            warnings.append(f"query_fasta could not be read: {exc}")
        else:
            for protein_id in ids:
                query_specs.append({"protein_id": protein_id, "old_locus_tag": "", "sequence": sequences.get(protein_id, ""), "description": descriptions.get(protein_id, ""), "source": "query_fasta"})

    return [query for query in query_specs if query["protein_id"] or query["old_locus_tag"] or query["sequence"]]


def _resolve_query(query: dict[str, str], index: int, records: dict[str, ProteinRecord]) -> dict[str, Any]:
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

    return {"query_id": query_id, "input_protein_id": query["protein_id"], "input_old_locus_tag": query["old_locus_tag"], "resolved_protein_id": resolved_protein_id, "resolved_old_locus_tag": resolved_old_locus_tag, "sequence": sequence, "sequence_length": len(sequence) if sequence else None, "resolution_status": status, "description": description, "notes": notes, "record": matched_record}


def _find_query_record(query: dict[str, str], records: dict[str, ProteinRecord]) -> ProteinRecord | None:
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
    return {key: query[key] for key in ("query_id", "input_protein_id", "input_old_locus_tag", "resolved_protein_id", "resolved_old_locus_tag", "sequence_length", "resolution_status", "description", "notes")}


def _rank_source_candidates(resolved_queries: list[dict[str, Any]], candidate_records: dict[str, ProteinRecord], candidate_source: str, scoring_config: Any, feature_map: dict[str, GffFeatureLocation]) -> list[dict[str, Any]]:
    all_rows: list[dict[str, Any]] = []
    for query in resolved_queries:
        if query["resolution_status"] != "resolved":
            continue
        query_rows: list[dict[str, Any]] = []
        for candidate in candidate_records.values():
            if _is_self_pair(query, candidate):
                continue
            query_rows.append(_score_pair(query, candidate, candidate_source, scoring_config, feature_map))
        _assign_distance_independent_ranks(query_rows)
        query_rows.sort(key=lambda row: (-float(row["interaction_priority_score"]), not bool(row["alphafold_recommended"]), str(row["candidate_protein_id"])))
        for rank, row in enumerate(query_rows[: scoring_config.max_candidates_per_query], start=1):
            row["candidate_rank"] = rank
            all_rows.append(row)
    all_rows.sort(key=lambda row: (str(row["query_id"]), int(row["candidate_rank"]), str(row["candidate_protein_id"])))
    return all_rows


def _score_pair(query: dict[str, Any], candidate: ProteinRecord, candidate_source: str, scoring_config: Any, feature_map: dict[str, GffFeatureLocation]) -> dict[str, Any]:
    weights = scoring_config.scoring_weights
    reasons: list[str] = [f"candidate source: {candidate_source}"]
    candidate_priority_score = _candidate_priority_score(candidate_source, weights.candidate_priority)
    neighborhood = _gene_neighborhood(query, candidate, feature_map, weights.gene_neighborhood)
    reasons.append(neighborhood["reason"])

    co_occurrence_score = _co_occurrence_score(query["record"], candidate, weights)
    reasons.append("similar source pattern" if co_occurrence_score else "no source pattern evidence")

    domain_complementarity_score = _domain_complementarity_score(query, candidate, weights.domain_complementarity, reasons)
    alphafold_readiness_score, pair_total_length, alphafold_recommended = _alphafold_readiness(query, candidate, scoring_config)
    reasons.append("compatible for manual AlphaFold" if alphafold_recommended else "missing sequence or length too large for manual AlphaFold")

    distance_independent_score = (
        candidate_priority_score
        + co_occurrence_score
        + domain_complementarity_score
    )
    priority_group = _priority_group(
        candidate_source=candidate_source,
        same_gene_neighborhood_score=neighborhood["same_gene_neighborhood_score"],
        co_occurrence_score=co_occurrence_score,
        domain_complementarity_score=domain_complementarity_score,
    )
    if priority_group in {"distant_cooccurrence_candidate", "distant_domain_candidate"}:
        reasons.append("distant candidate retained by co-occurrence/domain evidence")

    total_score = candidate_priority_score + neighborhood["same_gene_neighborhood_score"] + co_occurrence_score + domain_complementarity_score + alphafold_readiness_score
    row = {
        "query_id": query["query_id"],
        "query_protein_id": query["resolved_protein_id"],
        "query_old_locus_tag": query["resolved_old_locus_tag"],
        "candidate_rank": 0,
        "candidate_protein_id": candidate.protein_id,
        "candidate_old_locus_tag": candidate.old_locus_tag or "",
        "candidate_source": candidate_source,
        "candidate_description": candidate.description,
        **{key: value for key, value in neighborhood.items() if key != "reason"},
        "interaction_priority_score": round(total_score, 3),
        "distance_independent_score": round(distance_independent_score, 3),
        "distance_independent_rank": 0,
        "priority_group": priority_group,
        "interaction_score_reasons": "; ".join(reasons),
        "candidate_priority_score": round(candidate_priority_score, 3),
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


def _assign_distance_independent_ranks(rows: list[dict[str, Any]]) -> None:
    """Assign distance-independent ranks within one query/source sheet."""
    ranked_rows = sorted(
        rows,
        key=lambda row: (
            -float(row["distance_independent_score"]),
            -float(row["domain_complementarity_score"]),
            -float(row["co_occurrence_score"]),
            str(row["candidate_protein_id"]),
        ),
    )
    for rank, row in enumerate(ranked_rows, start=1):
        row["distance_independent_rank"] = rank


def _priority_group(
    candidate_source: str,
    same_gene_neighborhood_score: float,
    co_occurrence_score: float,
    domain_complementarity_score: float,
) -> str:
    if same_gene_neighborhood_score > 0:
        return "nearby_candidate"
    if co_occurrence_score > 0:
        return "distant_cooccurrence_candidate"
    if domain_complementarity_score > 0:
        return "distant_domain_candidate"
    if candidate_source == "No_hit":
        return "no_hit_candidate"
    return "general_candidate"


def _gene_neighborhood(query: dict[str, Any], candidate: ProteinRecord, feature_map: dict[str, GffFeatureLocation], max_score: float) -> dict[str, Any]:
    query_location = _record_location(query["resolved_protein_id"], query["resolved_old_locus_tag"], feature_map)
    candidate_location = _record_location(candidate.protein_id, candidate.old_locus_tag or "", feature_map)
    base = {"same_contig": None, "query_contig": None, "query_start": None, "query_end": None, "query_strand": None, "candidate_contig": None, "candidate_start": None, "candidate_end": None, "candidate_strand": None, "distance_bp": None, "strand_relation": "unknown", "same_gene_neighborhood_score": 0.0}
    if query_location is None or candidate_location is None:
        return {**base, "reason": "genomic distance unavailable"}

    base.update({"query_contig": query_location.contig, "query_start": query_location.start, "query_end": query_location.end, "query_strand": query_location.strand, "candidate_contig": candidate_location.contig, "candidate_start": candidate_location.start, "candidate_end": candidate_location.end, "candidate_strand": candidate_location.strand})
    base["strand_relation"] = _strand_relation(query_location.strand, candidate_location.strand)
    if query_location.contig != candidate_location.contig:
        return {**base, "same_contig": False, "reason": "different contig"}

    distance = _interval_distance(query_location, candidate_location)
    if distance <= 5000:
        score = max_score
        reason = f"close genomic neighborhood: {distance} bp"
    elif distance <= 20000:
        score = max_score * 15.0 / 25.0
        reason = f"moderate genomic neighborhood: {distance} bp"
    elif distance <= 100000:
        score = max_score * 5.0 / 25.0
        reason = f"weak genomic neighborhood: {distance} bp"
    else:
        score = 0.0
        reason = "distant genomic neighborhood"
    return {**base, "same_contig": True, "distance_bp": distance, "same_gene_neighborhood_score": round(score, 3), "reason": reason}


def _build_neighborhood_rows(
    resolved_queries: list[dict[str, Any]],
    source_rows: dict[str, list[dict[str, Any]]],
    scoring_config: Any,
) -> list[dict[str, Any]]:
    neighborhood_config = getattr(scoring_config, "neighborhood", None)
    if neighborhood_config is None or not getattr(neighborhood_config, "enabled", True):
        return []
    if not any(query["resolution_status"] == "resolved" for query in resolved_queries):
        return []

    max_distance = getattr(neighborhood_config, "max_distance_bp", 100000)
    max_rows_per_query = getattr(neighborhood_config, "max_rows_per_query", 200)
    query_by_id = {query["query_id"]: query for query in resolved_queries}
    seen: set[tuple[str, str, str]] = set()
    rows: list[dict[str, Any]] = []

    for interaction_rows in source_rows.values():
        for row in interaction_rows:
            distance = row.get("distance_bp")
            if distance is None or row.get("same_contig") is not True:
                continue
            if int(distance) > int(max_distance):
                continue
            key = (
                str(row.get("query_id", "")),
                str(row.get("candidate_protein_id", "")),
                str(row.get("candidate_source", "")),
            )
            if key in seen:
                continue
            seen.add(key)
            query = query_by_id.get(str(row.get("query_id", "")), {})
            rows.append({
                "query_id": row.get("query_id", ""),
                "query_protein_id": row.get("query_protein_id", ""),
                "query_old_locus_tag": row.get("query_old_locus_tag", ""),
                "query_description": query.get("description", ""),
                "query_contig": row.get("query_contig"),
                "query_start": row.get("query_start"),
                "query_end": row.get("query_end"),
                "query_strand": row.get("query_strand"),
                "candidate_rank_by_distance": 0,
                "candidate_protein_id": row.get("candidate_protein_id", ""),
                "candidate_old_locus_tag": row.get("candidate_old_locus_tag", ""),
                "candidate_description": row.get("candidate_description", ""),
                "candidate_source": row.get("candidate_source", ""),
                "candidate_contig": row.get("candidate_contig"),
                "candidate_start": row.get("candidate_start"),
                "candidate_end": row.get("candidate_end"),
                "candidate_strand": row.get("candidate_strand"),
                "distance_bp": distance,
                "strand_relation": row.get("strand_relation", "unknown"),
                "neighborhood_band": _neighborhood_band(distance, row.get("same_contig")),
                "same_gene_neighborhood_score": row.get("same_gene_neighborhood_score", 0.0),
                "interaction_priority_score": row.get("interaction_priority_score", 0.0),
                "domain_complementarity_score": row.get("domain_complementarity_score", 0.0),
                "candidate_priority_score": row.get("candidate_priority_score", 0.0),
                "co_occurrence_score": row.get("co_occurrence_score", 0.0),
                "alphafold_recommended": row.get("alphafold_recommended", False),
                "interaction_score_reasons": row.get("interaction_score_reasons", ""),
            })

    rows.sort(
        key=lambda row: (
            str(row["query_id"]),
            int(row["distance_bp"]),
            -float(row["interaction_priority_score"]),
            str(row["candidate_protein_id"]),
            str(row["candidate_source"]),
        )
    )

    limited_rows: list[dict[str, Any]] = []
    counts_by_query: dict[str, int] = {}
    for row in rows:
        query_id = str(row["query_id"])
        count = counts_by_query.get(query_id, 0)
        if count >= int(max_rows_per_query):
            continue
        row["candidate_rank_by_distance"] = count + 1
        counts_by_query[query_id] = count + 1
        limited_rows.append(row)
    return limited_rows


def _neighborhood_band(distance_bp: int | None, same_contig: bool | None) -> str:
    if same_contig is False:
        return "different_contig"
    if distance_bp is None:
        return "unknown"
    if distance_bp == 0:
        return "overlap"
    if distance_bp <= 5000:
        return "<=5kb"
    if distance_bp <= 20000:
        return "<=20kb"
    if distance_bp <= 100000:
        return "<=100kb"
    return ">100kb"


def _record_location(protein_id: str, old_locus_tag: str, feature_map: dict[str, GffFeatureLocation]) -> GffFeatureLocation | None:
    for key in (protein_id, _without_version(protein_id), old_locus_tag):
        if key and key in feature_map:
            return feature_map[key]
    return None


def _interval_distance(left: GffFeatureLocation, right: GffFeatureLocation) -> int:
    if left.end < right.start:
        return right.start - left.end
    if right.end < left.start:
        return left.start - right.end
    return 0


def _strand_relation(query_strand: str | None, candidate_strand: str | None) -> str:
    if not query_strand or not candidate_strand:
        return "unknown"
    return "same_strand" if query_strand == candidate_strand else "opposite_strand"


def _candidate_priority_score(candidate_source: str, weight: float) -> float:
    return CANDIDATE_PRIORITY_BASE.get(candidate_source, 10.0) / 30.0 * weight


def _co_occurrence_score(query_record: ProteinRecord | None, candidate: ProteinRecord, weights: Any) -> float:
    if query_record is None:
        return 0.0
    query_sources = set(query_record.positive_sources_hit)
    candidate_sources = set(candidate.positive_sources_hit)
    if query_sources or candidate_sources:
        union = query_sources | candidate_sources
        if union:
            return len(query_sources & candidate_sources) / len(union) * weights.co_occurrence
    if not query_record.negative_hits and not candidate.negative_hits:
        return weights.co_occurrence * 0.5
    return 0.0


def _domain_complementarity_score(query: dict[str, Any], candidate: ProteinRecord, weight: float, reasons: list[str]) -> float:
    query_text = _record_text(query["record"], query["description"])
    candidate_text = _record_text(candidate, candidate.description)
    if not query_text.strip() or not candidate_text.strip():
        reasons.append("no domain/description evidence")
        return 0.0

    query_terms = _meaningful_terms(query_text)
    candidate_terms = _meaningful_terms(candidate_text)
    if _has_domain_functional_terms(query["record"]) or _has_domain_functional_terms(candidate):
        reasons.append("Pfam/CDD functional terms used")
    shared = sorted((query_terms & candidate_terms) & MEANINGFUL_KEYWORDS)
    complementary = _complementary_terms(query_terms, candidate_terms)
    if complementary:
        reasons.append(f"complementary terms: {complementary[0]} + {complementary[1]}")
        return weight
    if len(shared) >= 2:
        reasons.append("meaningful shared terms: " + ", ".join(shared[:4]))
        return weight * 8.0 / 15.0
    if len(shared) == 1:
        reasons.append("meaningful shared terms: " + shared[0])
        return weight * 3.0 / 15.0

    generic_overlap = _all_terms(query_text) & _all_terms(candidate_text)
    if generic_overlap:
        reasons.append("generic-only description overlap ignored")
    else:
        reasons.append("no domain/description evidence")
    return 0.0


def _meaningful_terms(text: str) -> set[str]:
    normalized = _normalize_text(text)
    terms = {term for term in _all_terms(normalized) if term not in DESCRIPTION_STOPWORDS}
    for keyword in MEANINGFUL_KEYWORDS:
        if " " in keyword and keyword in normalized:
            terms.add(keyword)
    return terms


def _has_domain_functional_terms(record: ProteinRecord | None) -> bool:
    if record is None:
        return False
    for domain in record.domains:
        domain_text = _record_domain_text(domain)
        if _meaningful_terms(domain_text) & MEANINGFUL_KEYWORDS:
            return True
    return False


def _all_terms(text: str) -> set[str]:
    normalized = _normalize_text(text)
    return set(re.findall(r"[a-z][a-z0-9-]{2,}", normalized))


def _normalize_text(text: str) -> str:
    return re.sub(r"[^a-z0-9+-]+", " ", text.lower()).replace("+", "-")


def _complementary_terms(query_terms: set[str], candidate_terms: set[str]) -> tuple[str, str] | None:
    for left, right in COMPLEMENTARY_TERM_PAIRS:
        if (left in query_terms and right in candidate_terms) or (right in query_terms and left in candidate_terms):
            return left, right
    return None


def _alphafold_readiness(query: dict[str, Any], candidate: ProteinRecord, scoring_config: Any) -> tuple[float, int | None, bool]:
    if not query["sequence"] or not candidate.sequence:
        return 0.0, None, False
    pair_total_length = len(query["sequence"]) + len(candidate.sequence)
    recommended = pair_total_length <= scoring_config.alphafold.max_pair_total_length
    score = scoring_config.scoring_weights.alphafold_readiness if recommended else 0.0
    return score, pair_total_length, recommended


def _is_self_pair(query: dict[str, Any], candidate: ProteinRecord) -> bool:
    return bool(query["resolved_protein_id"]) and query["resolved_protein_id"] == candidate.protein_id


def _record_text(record: ProteinRecord | None, fallback_description: str) -> str:
    if record is None:
        return fallback_description
    domains = " ".join(_record_domain_text(domain) for domain in record.domains)
    return f"{record.description} {domains}"


def _record_domain_text(domain: Any) -> str:
    return " ".join(
        str(value)
        for value in (
            getattr(domain, "source", ""),
            getattr(domain, "accession", ""),
            getattr(domain, "name", ""),
            getattr(domain, "description", ""),
        )
        if value
    )


def _without_version(protein_id: str) -> str:
    return re.sub(r"\.\d+$", "", protein_id)


__all__: tuple[str, ...] = (
    "CANDIDATE_SOURCE_MAP",
    "INTERACTION_NEIGHBORHOOD_COLUMNS",
    "INTERACTION_NEIGHBORHOOD_SHEET",
    "INTERACTION_PAIR_COLUMNS",
    "INTERACTION_QUERY_COLUMNS",
    "INTERACTION_SHEET_DESCRIPTIONS",
    "InteractionScoringResult",
    "interaction_index_rows",
    "interaction_neighborhood_columns",
    "interaction_pair_columns",
    "interaction_sheet_name",
    "run_interaction_scoring",
    "source_sheet_name",
)
