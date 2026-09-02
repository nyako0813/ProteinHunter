"""Lightweight interaction candidate ranking helpers."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
import re
from typing import Any

from annotation.gff import GffFeatureLocation, load_gff_feature_map
from core.evidence import EvidenceComponent, EvidenceStatus, linear_normalize
from core.fasta import read_fasta_as_components
from core.models import CandidateScore, ProteinRecord
from analysis.candidates import get_best_hit
from analysis.scoring import build_candidate_score
from analysis.functional_complementarity_rules import (
    FunctionalComplementarityRuleset,
    load_functional_complementarity_ruleset,
)
from analysis.pih_evidence_bridge import (
    BRIDGED_PIH_CATEGORIES,
    PIH_CATEGORY_WEIGHTS,
    PihEvidenceBundle,
    load_pih_evidence_bundle,
)
from analysis.pih_evidence_bridge import without_version as _pih_without_version
from analysis.scoring_engine import ScoreBreakdown, rank_candidates, score_candidate
from analysis.scoring_engine_config import (
    ScoringEngineConfig,
    SequenceEvidenceConfig,
    load_scoring_engine_config,
)

V2_SCORING_MODEL = "v2_evidence_based"

# Per-component weight budget within each shared category for scoring model
# v2. genomic_context holds exactly one component today, so its weight value
# only needs to be > 0. source_classification now holds two components
# (source_classification itself and sequence_evidence, weighted evenly) that
# share its category cap, the same pattern the two functional_annotation
# components (co_occurrence, domain_complementarity) already use; see
# docs/implementation_plan_sequence_evidence.md for the reasoning.
V2_COMPONENT_WEIGHTS: dict[str, float] = {
    "source_classification": 1.0,
    "sequence_evidence": 1.0,
    "genomic_context": 1.0,
    "co_occurrence": 10.0,
    "domain_complementarity": 10.0,
    # Negative-evidence weight, expressed directly in output_scale points
    # (see analysis/scoring_engine.py::_negative_penalty): a "strong"
    # negative BLAST hit alone can fully consume the default
    # negative_penalty_cap (30). This reuses ortholog_filter.py's existing
    # strong/medium/weak classification -- it does not add a new
    # biological rule, only routes an existing v5 signal through the
    # auditable evidence model.
    "negative_hit_strength": 30.0,
}

_NEGATIVE_HIT_STRENGTH_VALUES: dict[str, float] = {
    "strong": 1.0,
    "medium": 0.5,
    "weak": 0.2,
}


@dataclass(frozen=True)
class InteractionScoringResult:
    """Interaction query rows and per-source pair rows for Excel output."""

    query_rows: list[dict[str, Any]]
    source_rows: dict[str, list[dict[str, Any]]]
    neighborhood_rows: list[dict[str, Any]]
    warnings: list[str]
    # Component-level (v2) or per-pair (legacy) breakdown rows for the
    # optional Interaction_Evidence_Detail sheet. Empty unless
    # interaction_scoring.evidence_detail_sheet actually produced rows for
    # at least one enabled candidate source.
    evidence_detail_rows: list[dict[str, Any]] = field(default_factory=list)
    evidence_detail_scoring_model: str = "legacy_additive"


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
    "Interaction_Evidence_Detail": (
        "evidence breakdown behind interaction_priority_score for the same pairs "
        "already shown in the Interaction_* sheets (scoring_model: "
        "v2_evidence_based -> one row per evidence component; legacy_additive "
        "-> one row per pair, same columns as the main sheets)",
        "shows exactly which evidence category/component drove a high or low score",
        "audit why a specific candidate ranked where it did before trusting the score",
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
    # scoring model v2 only (blank/NaN for scoring_model: legacy_additive)
    "scoring_model",
    "evidence_tier",
    "formal_score_available",
    "evidence_category_count",
    "evidence_component_count",
    "available_weight_total",
    # protein_hunter_score reference columns (M1/M2, design spec section 22):
    # the candidate's own query-independent "is this generally a good
    # candidate" score -- see resolve_protein_hunter_scores. Present for
    # both scoring models; never affects interaction_priority_score,
    # candidate_rank, or sort order.
    "protein_hunter_score",
    "protein_hunter_score_components",
    "protein_hunter_score_reasons",
    # interaction_score reference columns (M3/M4, design spec section 22):
    # query-specific evidence only (genomic_context + domain_complementarity
    # for v2, an equivalent re-normalized sum for legacy_additive;
    # co_occurrence deliberately excluded, see INTERACTION_SCORE_COMPONENT_NAMES).
    # scoring_model: v2_evidence_based only for interaction_evidence_tier --
    # legacy_additive has no per-category tiering concept, so it stays blank.
    # Like protein_hunter_score, purely additive: never affects
    # interaction_priority_score, candidate_rank, or sort order.
    "interaction_score",
    "interaction_evidence_tier",
)

SEQUENCE_COLUMNS: tuple[str, ...] = ("query_sequence", "candidate_sequence")

INTERACTION_EVIDENCE_DETAIL_SHEET = "Interaction_Evidence_Detail"

#: Long format: one row per (query, candidate, category, component). Only
#: produced for scoring_model: v2_evidence_based.
INTERACTION_EVIDENCE_DETAIL_V2_COLUMNS: tuple[str, ...] = (
    "query_id",
    "query_protein_id",
    "query_old_locus_tag",
    "candidate_protein_id",
    "candidate_old_locus_tag",
    "candidate_source",
    "candidate_rank",
    "category",
    "component_name",
    "status",
    "raw_value",
    "normalized_value",
    "weight",
    "category_cap",
    "is_negative",
    "explanation",
)

#: Wide format: one row per (query, candidate), projecting the same
#: breakdown columns already present on scoring_model: legacy_additive rows
#: in the main Interaction_* sheets (see INTERACTION_PAIR_COLUMNS above).
INTERACTION_EVIDENCE_DETAIL_LEGACY_COLUMNS: tuple[str, ...] = (
    "query_id",
    "query_protein_id",
    "query_old_locus_tag",
    "candidate_protein_id",
    "candidate_old_locus_tag",
    "candidate_source",
    "candidate_rank",
    "candidate_priority_score",
    "same_gene_neighborhood_score",
    "co_occurrence_score",
    "domain_complementarity_score",
    "alphafold_readiness_score",
    "interaction_score_reasons",
)


def interaction_evidence_detail_columns(scoring_model: str) -> tuple[str, ...]:
    """Return Interaction_Evidence_Detail columns for the given scoring model."""
    if scoring_model == V2_SCORING_MODEL:
        return INTERACTION_EVIDENCE_DETAIL_V2_COLUMNS
    return INTERACTION_EVIDENCE_DETAIL_LEGACY_COLUMNS


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
    evidence_detail_rows: list[dict[str, Any]] = []

    scoring_model = getattr(scoring_config, "scoring_model", "legacy_additive")
    evidence_detail_config = getattr(scoring_config, "evidence_detail_sheet", None)
    include_no_hit_detail = bool(getattr(evidence_detail_config, "include_no_hit", False))
    engine_config: ScoringEngineConfig | None = None
    ruleset: FunctionalComplementarityRuleset | None = None
    pih_bundle: PihEvidenceBundle | None = None
    if scoring_model == V2_SCORING_MODEL:
        engine_config = load_scoring_engine_config(
            getattr(scoring_config, "scoring_engine_config", None)
        )
        ruleset = load_functional_complementarity_ruleset(
            getattr(scoring_config, "functional_complementarity_ruleset", None)
        )
        pih_bundle_path = getattr(scoring_config, "pih_evidence_bundle", None)
        if pih_bundle_path is not None:
            pih_bundle = load_pih_evidence_bundle(pih_bundle_path)
            warnings.extend(pih_bundle.warnings)

    protein_hunter_scores = resolve_protein_hunter_scores(config, blast_classification)

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

        collect_evidence_detail = source_key != "no_hit" or include_no_hit_detail

        if scoring_model == V2_SCORING_MODEL:
            rows, detail_rows = _rank_source_candidates_v2(
                resolved_queries=resolved_queries,
                candidate_records=candidate_records,
                candidate_source=source_label,
                scoring_config=scoring_config,
                feature_map=feature_map,
                engine_config=engine_config,
                ruleset=ruleset,
                pih_bundle=pih_bundle,
                collect_evidence_detail=collect_evidence_detail,
            )
        else:
            rows, detail_rows = _rank_source_candidates(
                resolved_queries=resolved_queries,
                candidate_records=candidate_records,
                candidate_source=source_label,
                scoring_config=scoring_config,
                feature_map=feature_map,
                collect_evidence_detail=collect_evidence_detail,
            )
        if rows:
            _attach_protein_hunter_score_columns(rows, protein_hunter_scores)
            source_rows[sheet_name] = rows
            evidence_detail_rows.extend(detail_rows)
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
        evidence_detail_rows=evidence_detail_rows,
        evidence_detail_scoring_model=scoring_model,
    )


def _interaction_scoring_target_records(
    scoring_config: Any, blast_classification: Any
) -> dict[str, ProteinRecord]:
    """Return the de-duplicated union of every record interaction_scoring touches.

    Shared target-set resolution for both CDD annotation scope
    (``resolve_cdd_annotation_targets``) and protein_hunter_score scope
    (``resolve_protein_hunter_scores``): both need "every record reachable
    through an enabled ``candidate_sources`` bucket, plus every resolved
    query record." Returns an empty dict when ``scoring_config.enabled`` is
    false.
    """
    if not scoring_config.enabled:
        return {}

    targets: dict[str, ProteinRecord] = {}

    for source_key, enabled in scoring_config.candidate_sources.items():
        if not enabled:
            continue
        source_info = CANDIDATE_SOURCE_MAP.get(source_key)
        if source_info is None:
            continue
        _source_label, attr_name, _sheet_name = source_info
        candidate_records = getattr(blast_classification, attr_name, None)
        if not candidate_records:
            continue
        for protein_id, record in candidate_records.items():
            targets.setdefault(protein_id, record)

    discarded_warnings: list[str] = []
    queries = _load_query_specs(scoring_config, discarded_warnings)
    for index, query in enumerate(queries, start=1):
        resolved = _resolve_query(query, index, blast_classification.all_records)
        record = resolved.get("record")
        if record is not None:
            targets.setdefault(record.protein_id, record)

    return targets


def resolve_cdd_annotation_targets(
    config: Any, blast_classification: Any
) -> dict[str, ProteinRecord]:
    """Return extra records that should be CDD-annotated for interaction_scoring.

    CDD annotation (analysis/domain_annotator.py::annotate_records_cdd)
    normally only ever sees ``positive_only_records`` ("Candidates"), so
    interaction_scoring's own domain_complementarity evaluation could never
    see a real CDD hit for a query or candidate outside that one bucket
    (see docs/implementation_plan_sequence_evidence.md's CDD investigation
    notes for how this was discovered). This function returns the
    de-duplicated union of:

    - every record reachable through an enabled
      ``interaction_scoring.candidate_sources`` bucket (the same
      ``CANDIDATE_SOURCE_MAP`` lookup ``run_interaction_scoring`` itself
      uses), and
    - every interaction query that actually resolved to a target
      ``ProteinRecord`` (a query given only an explicit ``sequence``, with
      no matching target record, has no record object to annotate and is
      excluded).

    Returns an empty dict when ``interaction_scoring.enabled`` is false, so
    CDD's default scope (``positive_only_records`` only) is unaffected for
    runs that do not use interaction_scoring. The caller is expected to
    merge this into ``positive_only_records`` itself, which is always
    CDD-annotated regardless of interaction_scoring.
    """
    return _interaction_scoring_target_records(config.interaction_scoring, blast_classification)


def resolve_protein_hunter_scores(
    config: Any, blast_classification: Any
) -> dict[str, CandidateScore]:
    """Return protein_hunter_score for every record interaction_scoring touches.

    ``analysis/scoring.py::score_records`` (the "Candidate scoring" pipeline
    step) only ever scores ``positive_only_records`` ("Candidates"), so any
    interaction_scoring candidate source beyond that bucket -- Candidates_relaxed,
    No_hit, etc. -- previously had no protein_hunter_score at all (Excel showed
    a misleading ``total_score`` of 0, indistinguishable from "scored, and
    scored zero"). This reuses the same target-set logic as
    ``resolve_cdd_annotation_targets`` (see the design-spec section 22
    protein_hunter_score/interaction_score split) and scores every one of
    them with the exact same query-independent formula
    (``analysis/scoring.py::build_candidate_score``, unmodified).

    Deliberately does not mutate ``ProteinRecord.score``: those ProteinRecord
    objects are shared with the plain classification sheets (Candidates_relaxed,
    No_hit, ...), which must keep showing exactly what they showed before --
    only interaction_scoring's own reference columns should reflect this
    wider scope. Returns an empty dict when ``interaction_scoring.enabled``
    is false.
    """
    scoring_config = config.interaction_scoring
    targets = _interaction_scoring_target_records(scoring_config, blast_classification)
    return {
        protein_id: build_candidate_score(record) for protein_id, record in targets.items()
    }


def _attach_protein_hunter_score_columns(
    rows: list[dict[str, Any]], protein_hunter_scores: dict[str, CandidateScore]
) -> None:
    """Attach protein_hunter_score reference columns to each pair row in place.

    Purely additive/read-only relative to everything else on the row:
    never touches candidate_rank, interaction_priority_score, or any other
    existing field.
    """
    for row in rows:
        score = protein_hunter_scores.get(row["candidate_protein_id"])
        if score is None:
            row["protein_hunter_score"] = None
            row["protein_hunter_score_components"] = ""
            row["protein_hunter_score_reasons"] = ""
            continue
        row["protein_hunter_score"] = score.total_score
        row["protein_hunter_score_components"] = "; ".join(
            f"{name}={value}" for name, value in score.components.items()
        )
        row["protein_hunter_score_reasons"] = "; ".join(score.reasons)


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


def _rank_source_candidates(
    resolved_queries: list[dict[str, Any]],
    candidate_records: dict[str, ProteinRecord],
    candidate_source: str,
    scoring_config: Any,
    feature_map: dict[str, GffFeatureLocation],
    collect_evidence_detail: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    all_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    ranking_metric = getattr(scoring_config, "ranking_metric", "interaction_priority_score")
    # M5: legacy's interaction_score is always a plain float (never None --
    # see _legacy_interaction_score), so ranking by it needs no special
    # eligibility handling the way v2's rank_candidates does.
    sort_field = "interaction_score" if ranking_metric == "interaction_score" else "interaction_priority_score"

    for query in resolved_queries:
        if query["resolution_status"] != "resolved":
            continue
        query_rows: list[dict[str, Any]] = []
        for candidate in candidate_records.values():
            if _is_self_pair(query, candidate):
                continue
            query_rows.append(_score_pair(query, candidate, candidate_source, scoring_config, feature_map))
        _assign_distance_independent_ranks(query_rows)
        query_rows.sort(key=lambda row: (-float(row[sort_field]), not bool(row["alphafold_recommended"]), str(row["candidate_protein_id"])))
        for rank, row in enumerate(query_rows[: scoring_config.max_candidates_per_query], start=1):
            row["candidate_rank"] = rank
            all_rows.append(row)
            if collect_evidence_detail:
                detail_rows.append(_evidence_detail_row_legacy(row))
    all_rows.sort(key=lambda row: (str(row["query_id"]), int(row["candidate_rank"]), str(row["candidate_protein_id"])))
    return all_rows, detail_rows


def _evidence_detail_row_legacy(row: dict[str, Any]) -> dict[str, Any]:
    """Project one legacy_additive pair row onto Interaction_Evidence_Detail columns."""
    return {
        "query_id": row["query_id"],
        "query_protein_id": row["query_protein_id"],
        "query_old_locus_tag": row["query_old_locus_tag"],
        "candidate_protein_id": row["candidate_protein_id"],
        "candidate_old_locus_tag": row["candidate_old_locus_tag"],
        "candidate_source": row["candidate_source"],
        "candidate_rank": row["candidate_rank"],
        "candidate_priority_score": row["candidate_priority_score"],
        "same_gene_neighborhood_score": row["same_gene_neighborhood_score"],
        "co_occurrence_score": row["co_occurrence_score"],
        "domain_complementarity_score": row["domain_complementarity_score"],
        "alphafold_readiness_score": row["alphafold_readiness_score"],
        "interaction_score_reasons": row["interaction_score_reasons"],
    }


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
    interaction_score = _legacy_interaction_score(
        neighborhood["same_gene_neighborhood_score"], domain_complementarity_score, weights
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
        "interaction_score": interaction_score,
        # legacy_additive has no per-category tiering concept (no evidence
        # categories/status to count) -- left blank, same as v2's own
        # evidence_tier is left blank on legacy rows.
        "interaction_evidence_tier": None,
    }
    if scoring_config.include_sequences_in_excel:
        row["query_sequence"] = query["sequence"]
        row["candidate_sequence"] = candidate.sequence
    return row


def _legacy_interaction_score(
    same_gene_neighborhood_score: float, domain_complementarity_score: float, weights: Any
) -> float | None:
    """legacy_additive counterpart of interaction_score: query-specific evidence only.

    Same INTERACTION_SCORE_COMPONENT_NAMES scope as v2 (genomic_context +
    domain_complementarity, co_occurrence excluded) re-normalized to 0-100
    against the configured weight budget for those two components. Unlike
    v2, legacy_additive has no MISSING/AVAILABLE evidence-status concept --
    every legacy sub-score is always a plain number (0 when there is no
    evidence, never None) -- so this can only ever be a numeric 0-100
    value, not the "no evidence at all" None that v2's interaction_score
    can report. Returns None only when gene_neighborhood/domain_complementarity
    both carry zero weight in scoring_weights (nothing to normalize against).
    """
    max_points = weights.gene_neighborhood + weights.domain_complementarity
    if max_points <= 0:
        return None
    raw = same_gene_neighborhood_score + domain_complementarity_score
    return round(raw / max_points * 100, 3)


def _rank_source_candidates_v2(
    resolved_queries: list[dict[str, Any]],
    candidate_records: dict[str, ProteinRecord],
    candidate_source: str,
    scoring_config: Any,
    feature_map: dict[str, GffFeatureLocation],
    engine_config: ScoringEngineConfig,
    ruleset: FunctionalComplementarityRuleset,
    pih_bundle: PihEvidenceBundle | None = None,
    collect_evidence_detail: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Evidence-based (scoring model v2) counterpart of _rank_source_candidates."""
    all_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    ranking_metric = getattr(scoring_config, "ranking_metric", "interaction_priority_score")

    for query in resolved_queries:
        if query["resolution_status"] != "resolved":
            continue

        pairs: list[tuple[str, dict[str, Any], ScoreBreakdown, ScoreBreakdown]] = []
        for candidate in candidate_records.values():
            if _is_self_pair(query, candidate):
                continue
            row, breakdown, interaction_breakdown = _score_pair_v2(
                query, candidate, candidate_source, scoring_config, feature_map, engine_config, ruleset,
                pih_bundle=pih_bundle,
            )
            pairs.append((candidate.protein_id, row, breakdown, interaction_breakdown))

        # M5: candidate_rank/row order can be driven by either the full
        # composite breakdown (default, unchanged since Phase 5) or the
        # query-specific-only breakdown -- interaction_priority_score,
        # evidence_tier, and interaction_score themselves never change
        # based on this setting, only which one ranks.
        ranking_breakdown_by_id = {
            candidate_id: (interaction_breakdown if ranking_metric == "interaction_score" else breakdown)
            for candidate_id, _row, breakdown, interaction_breakdown in pairs
        }
        ranked = rank_candidates(
            list(ranking_breakdown_by_id.items()),
            tie_precision=engine_config.tie_precision,
        )
        row_by_id = {candidate_id: row for candidate_id, row, _breakdown, _interaction_breakdown in pairs}
        breakdown_by_id = {
            candidate_id: breakdown for candidate_id, _row, breakdown, _interaction_breakdown in pairs
        }
        for ranked_item in ranked[: scoring_config.max_candidates_per_query]:
            row = row_by_id[ranked_item.candidate_id]
            row["candidate_rank"] = ranked_item.rank if ranked_item.rank is not None else 0
            row["distance_independent_rank"] = row["candidate_rank"]
            all_rows.append(row)
            if collect_evidence_detail:
                # Evidence_Detail always audits the full breakdown, regardless
                # of which score is currently driving candidate_rank.
                full_breakdown = breakdown_by_id[ranked_item.candidate_id]
                detail_rows.extend(_evidence_detail_rows_v2(row, full_breakdown, engine_config))

    all_rows.sort(
        key=lambda row: (str(row["query_id"]), int(row["candidate_rank"] or 0), str(row["candidate_protein_id"]))
    )
    return all_rows, detail_rows


def _evidence_detail_rows_v2(
    row: dict[str, Any], breakdown: ScoreBreakdown, engine_config: ScoringEngineConfig
) -> list[dict[str, Any]]:
    """Expand one v2 pair's ScoreBreakdown into one Interaction_Evidence_Detail row per component."""
    detail_rows: list[dict[str, Any]] = []
    for component in breakdown.components:
        detail_rows.append(
            {
                "query_id": row["query_id"],
                "query_protein_id": row["query_protein_id"],
                "query_old_locus_tag": row["query_old_locus_tag"],
                "candidate_protein_id": row["candidate_protein_id"],
                "candidate_old_locus_tag": row["candidate_old_locus_tag"],
                "candidate_source": row["candidate_source"],
                "candidate_rank": row["candidate_rank"],
                "category": component.category,
                "component_name": component.name,
                "status": component.status.value,
                "raw_value": component.raw_value,
                "normalized_value": component.normalized_value,
                "weight": component.weight,
                "category_cap": engine_config.category_caps.get(component.category),
                "is_negative": component.is_negative,
                "explanation": component.explanation,
            }
        )
    return detail_rows


#: Components that constitute query-specific "does this candidate actually
#: interact with THIS query" evidence (design spec section 22's
#: interaction_score). Deliberately excludes source_classification and
#: sequence_evidence (candidate-only conservation quality, see
#: protein_hunter_score) and negative_hit_strength (candidate-only
#: penalty). co_occurrence is also excluded: despite technically taking the
#: query as an input, it measures each protein's own BLAST hit pattern
#: against the positive reference genomes independently, not any
#: relationship between query and candidate -- and with only a handful of
#: configured positive sources its Jaccard value is coarse-grained enough
#: (this project's default config has exactly two, so it can only be 0.0,
#: 0.5, or 1.0) that it behaves like a candidate-quality signal in
#: practice. co_occurrence is still shown, unchanged, in
#: Interaction_Evidence_Detail for audit -- it is only excluded from this
#: sum. PIH-bridged categories are intentionally left out of this first cut
#: (all pih_* categories are query-specific in principle, but the scope for
#: this phase was fixed to genomic_context + domain_complementarity only).
INTERACTION_SCORE_COMPONENT_NAMES: frozenset[str] = frozenset({"genomic_context", "domain_complementarity"})


def _interaction_only_breakdown(
    components: list[EvidenceComponent], engine_config: ScoringEngineConfig
) -> ScoreBreakdown:
    """Re-score a pair using only query-specific (interaction_score) components.

    Reuses analysis/scoring_engine.py::score_candidate unmodified: feeding it
    a filtered component list is enough to get a correctly cap-renormalized
    0-100 score, tier, and eligibility for the restricted view -- no new
    scoring engine code needed.
    """
    interaction_components = [c for c in components if c.name in INTERACTION_SCORE_COMPONENT_NAMES]
    return score_candidate(interaction_components, engine_config)


def _score_pair_v2(
    query: dict[str, Any],
    candidate: ProteinRecord,
    candidate_source: str,
    scoring_config: Any,
    feature_map: dict[str, GffFeatureLocation],
    engine_config: ScoringEngineConfig,
    ruleset: FunctionalComplementarityRuleset,
    pih_bundle: PihEvidenceBundle | None = None,
) -> tuple[dict[str, Any], ScoreBreakdown, ScoreBreakdown]:
    """Score one query/candidate pair with the evidence-based engine.

    Returns ``(row, breakdown, interaction_breakdown)``: the full composite
    breakdown (drives interaction_priority_score/evidence_tier, unchanged
    since Phase 5) and the query-specific-only breakdown (drives
    interaction_score/interaction_evidence_tier, and -- when
    ``ranking_metric: interaction_score`` -- candidate_rank itself; see
    ``_rank_source_candidates_v2``).
    """
    components, location_info = _build_evidence_components_v2(
        query, candidate, candidate_source, feature_map, ruleset, engine_config, pih_bundle=pih_bundle
    )
    breakdown = score_candidate(components, engine_config)
    interaction_breakdown = _interaction_only_breakdown(components, engine_config)

    alphafold_readiness_score, pair_total_length, alphafold_recommended = _alphafold_readiness(
        query, candidate, scoring_config
    )

    reasons: list[str] = [f"candidate source: {candidate_source}", "scoring model: v2_evidence_based"]
    for component in components:
        if component.explanation:
            reasons.append(f"{component.name}: {component.explanation}")
    reasons.append(
        "compatible for manual AlphaFold"
        if alphafold_recommended
        else "missing sequence or length too large for manual AlphaFold "
        "(reference only; not part of the v2 total score)"
    )
    if not breakdown.eligible:
        reasons.append("insufficient evidence for a formal score; see evidence_tier=Unclassified")

    final_score = breakdown.final_score
    genomic_context_category = breakdown.category_scores.get("genomic_context")
    source_category = breakdown.category_scores.get("source_classification")

    row: dict[str, Any] = {
        "query_id": query["query_id"],
        "query_protein_id": query["resolved_protein_id"],
        "query_old_locus_tag": query["resolved_old_locus_tag"],
        "candidate_rank": 0,
        "candidate_protein_id": candidate.protein_id,
        "candidate_old_locus_tag": candidate.old_locus_tag or "",
        "candidate_source": candidate_source,
        "candidate_description": candidate.description,
        **location_info,
        "same_gene_neighborhood_score": round(genomic_context_category.capped_score, 3)
        if genomic_context_category is not None
        else None,
        "interaction_priority_score": round(final_score, 3) if final_score is not None else None,
        "distance_independent_score": round(final_score, 3) if final_score is not None else None,
        "distance_independent_rank": 0,
        "priority_group": breakdown.tier,
        "interaction_score_reasons": "; ".join(reasons),
        "candidate_priority_score": round(source_category.capped_score, 3)
        if source_category is not None
        else None,
        "co_occurrence_score": _component_contribution(components, "co_occurrence"),
        "domain_complementarity_score": _component_contribution(components, "domain_complementarity"),
        "alphafold_readiness_score": round(alphafold_readiness_score, 3),
        "pair_total_length": pair_total_length,
        "alphafold_recommended": alphafold_recommended,
        "scoring_model": V2_SCORING_MODEL,
        "evidence_tier": breakdown.tier,
        "formal_score_available": breakdown.eligible,
        "evidence_category_count": breakdown.evidence_category_count,
        "evidence_component_count": breakdown.evidence_component_count,
        "available_weight_total": round(breakdown.available_weight_total, 3),
        "interaction_score": round(interaction_breakdown.final_score, 3)
        if interaction_breakdown.final_score is not None
        else None,
        "interaction_evidence_tier": interaction_breakdown.tier,
    }
    if scoring_config.include_sequences_in_excel:
        row["query_sequence"] = query["sequence"]
        row["candidate_sequence"] = candidate.sequence
    return row, breakdown, interaction_breakdown


def _component_contribution(components: list[EvidenceComponent], name: str) -> float | None:
    """Return a named component's weighted contribution, or None if unavailable."""
    component = next((c for c in components if c.name == name), None)
    if component is None or component.status is not EvidenceStatus.AVAILABLE:
        return None
    return round(component.contribution, 3)


def _build_evidence_components_v2(
    query: dict[str, Any],
    candidate: ProteinRecord,
    candidate_source: str,
    feature_map: dict[str, GffFeatureLocation],
    ruleset: FunctionalComplementarityRuleset,
    engine_config: ScoringEngineConfig,
    pih_bundle: PihEvidenceBundle | None = None,
) -> tuple[list[EvidenceComponent], dict[str, Any]]:
    """Build the evidence components for one pair, reusing v5's raw signals."""
    components: list[EvidenceComponent] = []

    source_value = min(1.0, CANDIDATE_PRIORITY_BASE.get(candidate_source, 10.0) / 30.0)
    components.append(
        EvidenceComponent.available(
            "source_classification",
            "source_classification",
            source_value,
            V2_COMPONENT_WEIGHTS["source_classification"],
            raw_value=candidate_source,
            source="blast_classification",
            explanation=f"candidate source: {candidate_source}",
        )
    )

    seq_status, seq_value, seq_reason = _sequence_evidence_status_and_value(
        candidate, engine_config.sequence_evidence
    )
    if seq_status is EvidenceStatus.AVAILABLE:
        components.append(
            EvidenceComponent.available(
                "sequence_evidence",
                "source_classification",
                seq_value,
                V2_COMPONENT_WEIGHTS["sequence_evidence"],
                raw_value=get_best_hit(candidate.positive_hits).percent_identity,
                source="blast_hit",
                explanation=seq_reason,
            )
        )
    else:
        components.append(
            EvidenceComponent.unavailable(
                "sequence_evidence",
                "source_classification",
                seq_status,
                source="blast_hit",
                explanation=seq_reason,
            )
        )

    location_info, geo_status, geo_value, geo_reason = _gene_neighborhood_v2(query, candidate, feature_map)
    if geo_status is EvidenceStatus.AVAILABLE:
        components.append(
            EvidenceComponent.available(
                "genomic_context",
                "genomic_context",
                geo_value,
                V2_COMPONENT_WEIGHTS["genomic_context"],
                raw_value=location_info.get("distance_bp"),
                source="gff",
                explanation=geo_reason,
            )
        )
    else:
        components.append(
            EvidenceComponent.unavailable(
                "genomic_context", "genomic_context", geo_status, source="gff", explanation=geo_reason
            )
        )

    co_status, co_value, co_reason = _co_occurrence_status_and_value(query["record"], candidate)
    if co_status is EvidenceStatus.AVAILABLE:
        components.append(
            EvidenceComponent.available(
                "co_occurrence",
                "functional_annotation",
                co_value,
                V2_COMPONENT_WEIGHTS["co_occurrence"],
                source="blast_sources",
                explanation=co_reason,
            )
        )
    else:
        components.append(
            EvidenceComponent.unavailable(
                "co_occurrence", "functional_annotation", co_status, source="blast_sources", explanation=co_reason
            )
        )

    dom_status, dom_value, dom_reason = _domain_complementarity_status_and_value(query, candidate, ruleset)
    if dom_status is EvidenceStatus.AVAILABLE:
        components.append(
            EvidenceComponent.available(
                "domain_complementarity",
                "functional_annotation",
                dom_value,
                V2_COMPONENT_WEIGHTS["domain_complementarity"],
                source="annotation_text",
                explanation=dom_reason,
            )
        )
    else:
        components.append(
            EvidenceComponent.unavailable(
                "domain_complementarity",
                "functional_annotation",
                dom_status,
                source="annotation_text",
                explanation=dom_reason,
            )
        )

    neg_status, neg_value, neg_reason = _negative_hit_status_and_value(candidate)
    if neg_status is EvidenceStatus.AVAILABLE:
        components.append(
            EvidenceComponent.available(
                "negative_hit_strength",
                "source_reliability",
                neg_value,
                V2_COMPONENT_WEIGHTS["negative_hit_strength"],
                raw_value=candidate.negative_hit_strength,
                is_negative=True,
                source="ortholog_filter",
                explanation=neg_reason,
            )
        )
    else:
        components.append(
            EvidenceComponent.unavailable(
                "negative_hit_strength",
                "source_reliability",
                neg_status,
                source="ortholog_filter",
                explanation=neg_reason,
            )
        )

    if pih_bundle is not None:
        query_keys = [
            query["resolved_protein_id"],
            query["resolved_old_locus_tag"],
            _pih_without_version(query["resolved_protein_id"]),
        ]
        candidate_keys = [
            candidate.protein_id,
            candidate.old_locus_tag or "",
            _pih_without_version(candidate.protein_id),
        ]
        pih_categories = pih_bundle.lookup(query_keys, candidate_keys)
        for category_name in BRIDGED_PIH_CATEGORIES:
            evidence = pih_categories.get(category_name)
            if evidence is None:
                # No "unavailable" placeholder is needed here: a category with
                # no available evidence simply never appears in this pair's
                # component list, and the scoring engine already excludes
                # inactive categories from the score denominator (see
                # analysis/scoring_engine.py::_score_categories).
                continue
            components.append(
                EvidenceComponent.available(
                    f"pih_{category_name}",
                    f"pih_{category_name}",
                    evidence.normalized_score,
                    PIH_CATEGORY_WEIGHTS[category_name],
                    raw_value=evidence.available_weight,
                    source="protein_interaction_hunter",
                    explanation=(
                        f"PIH {category_name} category "
                        f"(available_weight={evidence.available_weight:.1f})"
                    ),
                )
            )

    return components, location_info


def _sequence_evidence_status_and_value(
    candidate: ProteinRecord, cfg: SequenceEvidenceConfig
) -> tuple[EvidenceStatus, float | None, str]:
    """Normalize the candidate's best positive BLAST hit into a 0.0-1.0 strength value.

    Reuses analysis.candidates.get_best_hit's existing (bitscore, then
    lowest evalue) representative-hit rule -- the same rule already used by
    output/excel.py's best_positive_hit columns -- instead of inventing a
    new way to aggregate multiple positive_hits. A candidate with no
    positive BLAST hit at all (e.g. the No_hit source) is MISSING, not a
    scored zero: "evaluated, weak hit" and "never evaluated" must stay
    distinguishable (see docs/implementation_plan_sequence_evidence.md).
    """
    if not candidate.positive_hits:
        return EvidenceStatus.MISSING, None, "no positive BLAST hit"

    hit = get_best_hit(candidate.positive_hits)

    identity_score = linear_normalize(hit.percent_identity, cfg.identity_floor, cfg.identity_ceiling)
    evalue_score = _evalue_strength_score(hit.evalue, cfg)
    weighted_sum = cfg.identity_weight * identity_score + cfg.evalue_weight * evalue_score
    total_weight = cfg.identity_weight + cfg.evalue_weight

    coverage = hit.query_coverage
    if coverage is None:
        coverage_note = "coverage unavailable"
    else:
        coverage_score = linear_normalize(coverage, cfg.coverage_floor, cfg.coverage_ceiling)
        weighted_sum += cfg.coverage_weight * coverage_score
        total_weight += cfg.coverage_weight
        coverage_note = f"coverage={coverage:.1f}%"

    normalized_value = weighted_sum / total_weight
    reason = (
        f"best positive BLAST hit: identity={hit.percent_identity:.1f}%, "
        f"{coverage_note}, evalue={hit.evalue:.2e} -> strength={normalized_value:.2f}"
    )
    return EvidenceStatus.AVAILABLE, normalized_value, reason


def _evalue_strength_score(evalue: float, cfg: SequenceEvidenceConfig) -> float:
    """Map a BLAST e-value onto 0.0-1.0 strength via a log-scale normalization.

    e-values are exponentially distributed and already upper-bounded by the
    BLAST evalue cutoff used to produce positive_hits in the first place, so
    a plain linear scale would be meaningless; -log10(evalue) is normalized
    instead. e-value <= 0 (BLAST can report an exact 0.0 for extremely
    significant hits, a floating-point underflow) is treated directly as the
    strongest possible signal, since -log10(0) is mathematically undefined.
    """
    if evalue <= 0:
        return 1.0
    return linear_normalize(
        -math.log10(evalue),
        -math.log10(cfg.evalue_reference_ceiling),
        -math.log10(cfg.evalue_reference_floor),
    )


def _negative_hit_status_and_value(candidate: ProteinRecord) -> tuple[EvidenceStatus, float | None, str]:
    """Reuse ortholog_filter.py's negative-hit strength as negative evidence.

    ortholog_filter.py itself is unchanged; this only reads its already
    computed classification (record.negative_hit_strength) and lets the v2
    scoring engine apply it as a capped penalty instead of an early hard
    filter, for candidate sources (e.g. Candidates_relaxed, No_hit) that
    intentionally retain records with a weak/medium negative hit.
    """
    strength = candidate.negative_hit_strength
    value = _NEGATIVE_HIT_STRENGTH_VALUES.get(strength)
    if value is None:
        return EvidenceStatus.NOT_APPLICABLE, None, "no negative BLAST hit"
    return EvidenceStatus.AVAILABLE, value, f"negative BLAST hit strength: {strength}"


def _gene_neighborhood_v2(
    query: dict[str, Any], candidate: ProteinRecord, feature_map: dict[str, GffFeatureLocation]
) -> tuple[dict[str, Any], EvidenceStatus, float | None, str]:
    """Evidence-status-aware counterpart of _gene_neighborhood."""
    query_location = _record_location(query["resolved_protein_id"], query["resolved_old_locus_tag"], feature_map)
    candidate_location = _record_location(candidate.protein_id, candidate.old_locus_tag or "", feature_map)
    base = {
        "same_contig": None,
        "query_contig": None,
        "query_start": None,
        "query_end": None,
        "query_strand": None,
        "candidate_contig": None,
        "candidate_start": None,
        "candidate_end": None,
        "candidate_strand": None,
        "distance_bp": None,
        "strand_relation": "unknown",
    }
    if query_location is None or candidate_location is None:
        return base, EvidenceStatus.MISSING, None, "genomic coordinates unavailable"

    base.update(
        {
            "query_contig": query_location.contig,
            "query_start": query_location.start,
            "query_end": query_location.end,
            "query_strand": query_location.strand,
            "candidate_contig": candidate_location.contig,
            "candidate_start": candidate_location.start,
            "candidate_end": candidate_location.end,
            "candidate_strand": candidate_location.strand,
        }
    )
    base["strand_relation"] = _strand_relation(query_location.strand, candidate_location.strand)
    if query_location.contig != candidate_location.contig:
        return {**base, "same_contig": False}, EvidenceStatus.NOT_APPLICABLE, None, "different contig"

    distance = _interval_distance(query_location, candidate_location)
    if distance <= 5000:
        normalized_value = 1.0
        reason = f"close genomic neighborhood: {distance} bp"
    elif distance <= 20000:
        normalized_value = 0.6
        reason = f"moderate genomic neighborhood: {distance} bp"
    elif distance <= 100000:
        normalized_value = 0.2
        reason = f"weak genomic neighborhood: {distance} bp"
    else:
        normalized_value = 0.0
        reason = "distant genomic neighborhood"
    return {**base, "same_contig": True, "distance_bp": distance}, EvidenceStatus.AVAILABLE, normalized_value, reason


def _co_occurrence_status_and_value(
    query_record: ProteinRecord | None, candidate: ProteinRecord
) -> tuple[EvidenceStatus, float | None, str]:
    """Evidence-status-aware counterpart of _co_occurrence_score."""
    if query_record is None:
        return EvidenceStatus.MISSING, None, "query has no BLAST classification record"

    query_sources = set(query_record.positive_sources_hit)
    candidate_sources = set(candidate.positive_sources_hit)
    union = query_sources | candidate_sources
    if union:
        jaccard = len(query_sources & candidate_sources) / len(union)
        return EvidenceStatus.AVAILABLE, jaccard, f"positive-source overlap (Jaccard): {jaccard:.2f}"
    if not query_record.negative_hits and not candidate.negative_hits:
        return (
            EvidenceStatus.AVAILABLE,
            0.5,
            "no positive-source overlap evidence, but neither side has a negative BLAST hit",
        )
    return (
        EvidenceStatus.AVAILABLE,
        0.0,
        "evaluated: no positive-source overlap and at least one side has a negative BLAST hit",
    )


def _domain_complementarity_status_and_value(
    query: dict[str, Any], candidate: ProteinRecord, ruleset: FunctionalComplementarityRuleset
) -> tuple[EvidenceStatus, float | None, str]:
    """Evidence-status-aware counterpart of _domain_complementarity_score."""
    query_text = _record_text(query["record"], query["description"])
    candidate_text = _record_text(candidate, candidate.description)
    if not query_text.strip() or not candidate_text.strip():
        return EvidenceStatus.MISSING, None, "no domain/description evidence"

    query_terms = _meaningful_terms_v2(query_text, ruleset)
    candidate_terms = _meaningful_terms_v2(candidate_text, ruleset)
    rule_note = ""
    if _has_domain_functional_terms_v2(query["record"], ruleset) or _has_domain_functional_terms_v2(
        candidate, ruleset
    ):
        rule_note = "Pfam/CDD functional terms used; "

    match = ruleset.find_match(query_terms, candidate_terms)
    if match is not None:
        return EvidenceStatus.AVAILABLE, 1.0, f"{rule_note}complementary rule matched: {match.rule_id}"

    shared = sorted((query_terms & candidate_terms) & ruleset.meaningful_keywords)
    if len(shared) >= 2:
        return EvidenceStatus.AVAILABLE, 8.0 / 15.0, f"{rule_note}meaningful shared terms: {', '.join(shared[:4])}"
    if len(shared) == 1:
        return EvidenceStatus.AVAILABLE, 3.0 / 15.0, f"{rule_note}meaningful shared term: {shared[0]}"

    generic_overlap = _all_terms(query_text) & _all_terms(candidate_text)
    if generic_overlap:
        return EvidenceStatus.AVAILABLE, 0.0, f"{rule_note}generic-only description overlap ignored"
    return EvidenceStatus.AVAILABLE, 0.0, f"{rule_note}no domain/description match"


def _meaningful_terms_v2(text: str, ruleset: FunctionalComplementarityRuleset) -> set[str]:
    normalized = _normalize_text(text)
    terms = {term for term in _all_terms(normalized) if term not in ruleset.stopwords}
    for keyword in ruleset.meaningful_keywords:
        if " " in keyword and keyword in normalized:
            terms.add(keyword)
    return terms


def _has_domain_functional_terms_v2(record: ProteinRecord | None, ruleset: FunctionalComplementarityRuleset) -> bool:
    if record is None:
        return False
    for domain in record.domains:
        domain_text = _record_domain_text(domain)
        if _meaningful_terms_v2(domain_text, ruleset) & ruleset.meaningful_keywords:
            return True
    return False


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
    "resolve_cdd_annotation_targets",
    "run_interaction_scoring",
    "source_sheet_name",
)
