"""Deterministic, LLM-free narrative text for the Word report (Phase 6-8 Stage 2).

Every function here is a pure function of a row dict (already shaped like a
``04_Score_Breakdown`` row, see ``output/report_v2.py``/``output/excel.py``)
plus a small amount of run-level context (category cap values, whether a
PIH evidence bundle was configured) -- no LLM call, no network call, no
randomness, so the same input always produces the same string (design spec
section 45, reproducibility). See
``claude/phase678_stage2_word_report_investigation.md`` section 3 for the
design rationale and the original template drafts this module implements.

This module deliberately has no knowledge of ``python-docx``, ``config.py``,
or ``analysis/scoring_engine_config.py`` -- callers (``output/word_report.py``)
resolve category caps and the PIH-bundle flag once per run and pass in plain
values, keeping this module importable and testable without either
``python-docx`` or a real pipeline config.
"""

from __future__ import annotations

from typing import Any, NamedTuple, Sequence


class CategoryRef(NamedTuple):
    """One evidence category as it appears on a pair row.

    ``kind`` groups row fields that represent the same underlying evidence
    category under different names in the two scoring models (e.g. v2's
    single ``functional_domain_score`` vs. legacy_additive's separate
    ``co_occurrence_score``/``domain_complementarity_score``) -- used to
    avoid printing the same "biological color" sentence twice for one
    category. ``label``/``cap`` are display-only.
    """

    row_key: str
    kind: str
    label: str
    cap: float


_TIER_OPENING: dict[str, str] = {
    "Tier1_VeryStrong": (
        "This candidate reached the highest confidence tier for this query "
        "(Tier 1 — Very Strong), with a Final Score of {final_score:.1f}/100 "
        "supported by evidence from {evidence_category_count} independent "
        "evidence categories."
    ),
    "Tier2_Strong": (
        "This candidate reached the second-highest confidence tier for this "
        "query (Tier 2 — Strong), with a Final Score of {final_score:.1f}/100 "
        "supported by evidence from {evidence_category_count} independent "
        "evidence categories."
    ),
    "Tier3_Moderate": (
        "This candidate reached a moderate confidence tier for this query "
        "(Tier 3 — Moderate), with a Final Score of {final_score:.1f}/100. "
        "The supporting evidence is present but comes from fewer independent "
        "categories ({evidence_category_count}) than higher-tier candidates."
    ),
    "Tier4_Weak": (
        "This candidate's Final Score of {final_score:.1f}/100 places it in "
        "the weakest confidence tier evaluated (Tier 4 — Weak). Its position "
        "in this list reflects a relative rank among the candidates scored "
        "for this query, not strong independent support."
    ),
}

_UNCLASSIFIED_OPENING = (
    "This candidate did not receive a formal Final Score for this query "
    "(insufficient evidence to meet the pipeline's minimum evidence "
    "requirement). It appears here only because it was among the "
    "higher-ranked candidates by whatever partial evidence was available; "
    "treat its inclusion as informational, not as a ranked, scored result."
)

_CANDIDATE_SOURCE_SENTENCES: dict[str, str] = {
    "Candidates": (
        "This protein was classified as a strict positive candidate — a "
        "BLAST hit to a positive reference sequence with no negative-reference "
        "hit at all, the most stringent classification this pipeline produces."
    ),
    "Positive_all_sources": (
        "This protein hit every configured positive reference source, with "
        "no negative-reference hit."
    ),
    "Candidates_relaxed": (
        "This protein was classified as a positive candidate under relaxed "
        "criteria — it has a positive-reference hit, and tolerates a "
        "non-strong (medium or weak) negative-reference hit."
    ),
    "No_hit": (
        "This protein had no BLAST hit to either the positive or negative "
        "reference sets. It is a lineage-specific or novel candidate not "
        "detected by sequence homology alone — its ranking here depends "
        "entirely on non-sequence evidence (genomic context, domain "
        "annotation, interaction evidence), not on a BLAST match."
    ),
    "Negative_unmatched": (
        "This protein has no negative-reference hit, though it did not meet "
        "the stricter positive-match criteria above."
    ),
    "Negative_hit": (
        "This protein has at least one negative-reference BLAST hit "
        "(strength: {negative_hit_strength}). It is included in this "
        "ranking despite that hit because of the other evidence enumerated "
        "above; interpret its position with additional caution."
    ),
}

_NEGATIVE_HIT_CAVEAT = (
    "Note: this candidate also has a {negative_hit_strength} BLAST hit to a "
    "negative-reference sequence, which does not change its candidate_source "
    "classification but is a relevant caveat."
)

#: "Biological color" sentences for the Biological Interpretation section,
#: keyed by CategoryRef.kind. Deliberately hedged ("consistent with",
#: never "confirms"/"proves"/"shows that X is Y") per design spec section
#: 35. "sequence" has no entry -- BLAST-based classification is procedural,
#: not itself biological-interpretation material (it is already covered in
#: "why this candidate ranks highly").
_CATEGORY_COLOR_SENTENCES: dict[str, str] = {
    "genomic_context": (
        "Its genomic proximity to the query gene is consistent with — but "
        "does not by itself establish — a shared operon or functional "
        "module."
    ),
    "functional_domain": (
        "Shared or complementary functional domain annotations were found "
        "between this candidate and the query, consistent with a related or "
        "complementary biochemical role; domain similarity alone does not "
        "establish direct interaction."
    ),
    "interaction": (
        "External interaction evidence (protein-protein interaction "
        "database and/or coexpression data) supports a functional "
        "association between this candidate and the query, independent of "
        "sequence similarity or genomic position."
    ),
    "evolutionary": (
        "Evolutionary/phylogenetic profile evidence from the optional PIH "
        "bridge supports consistency between this candidate and the query "
        "across the genomes compared."
    ),
    "cellular_compatibility": (
        "Cellular-compatibility evidence from the optional PIH bridge "
        "supports this candidate and the query being compatible with the "
        "same subcellular context."
    ),
}

_NEGATIVE_EVIDENCE_CLOSER = (
    "Negative (biological contradiction) evidence is a reserved category "
    "with no implemented signal in this pipeline version (see Evidence "
    "Architecture, section 5.7) — its absence here does not mean this "
    "candidate was checked for contradicting evidence and cleared."
)


def _is_contributing(value: Any) -> bool:
    """Whether a category score counts as "contributed to this ranking" for display.

    None means "not evaluated" (the established convention for every
    category reference column, see output/excel.py's
    INTERACTION_SCORE_EXPLANATIONS). A present-but-zero value is treated as
    "evaluated, no signal" for legacy_additive rows (which cannot otherwise
    distinguish missing from zero) and simply as "no meaningful
    contribution" for v2 rows either way -- in both cases, not worth
    enumerating as a reason this candidate ranks highly.
    """
    return value is not None and value > 0


def build_why_ranks_highly(row: dict[str, Any], category_refs: Sequence[CategoryRef]) -> str:
    """Build the deterministic "why this candidate ranks highly" paragraph.

    ``category_refs`` should already be resolved for this row's
    ``scoring_model`` (see output/word_report.py::category_refs_for_scoring_model)
    -- this function only reads ``row`` and the refs, no model branching.
    """
    tier = row.get("final_score_tier")
    if tier in _TIER_OPENING:
        opening = _TIER_OPENING[tier].format(
            final_score=row.get("final_score") or 0.0,
            evidence_category_count=row.get("evidence_category_count") or 0,
        )
    else:
        opening = _UNCLASSIFIED_OPENING

    contributing = [ref for ref in category_refs if _is_contributing(row.get(ref.row_key))]
    if contributing:
        listing = ", ".join(
            f"{ref.label} ({row[ref.row_key]:.1f}/{ref.cap:.0f})" for ref in contributing
        )
        evidence_sentence = f"Contributing evidence categories: {listing}."
    else:
        evidence_sentence = (
            "No individual evidence category scored above zero for this "
            "candidate under the current run's configuration."
        )

    candidate_source = str(row.get("candidate_source") or "")
    negative_hit_strength = str(row.get("negative_hit_strength") or "none")
    source_sentence = _CANDIDATE_SOURCE_SENTENCES.get(candidate_source, "").format(
        negative_hit_strength=negative_hit_strength
    )

    parts = [opening, evidence_sentence]
    if source_sentence:
        parts.append(source_sentence)
    if negative_hit_strength != "none" and candidate_source != "Negative_hit":
        parts.append(_NEGATIVE_HIT_CAVEAT.format(negative_hit_strength=negative_hit_strength))
    return " ".join(parts)


def build_evolutionary_closer(scoring_model: str, pih_bundle_configured: bool) -> str:
    """Return the run-level sentence explaining Evolutionary/Cellular Compatibility coverage.

    Constant for every candidate in one run (depends only on
    ``scoring_model``/``pih_bundle_configured``, never on the row) -- callers
    may compute it once per run and reuse it rather than recomputing per
    candidate, though calling it repeatedly is also correct (it is a pure
    function).
    """
    if scoring_model != "v2_evidence_based":
        return (
            "This run used the legacy_additive scoring model, which does "
            "not evaluate Evolutionary or Cellular Compatibility evidence "
            "at all — those categories exist only under the "
            "v2_evidence_based scoring model."
        )
    if pih_bundle_configured:
        return (
            "Evolutionary and Cellular Compatibility evidence were "
            "evaluated for this run using a supplied ProteinInteractionHunter "
            "(PIH) data bundle; a blank value for either means no matching "
            "evidence was found for this specific candidate, not that the "
            "category was skipped."
        )
    return (
        "Evolutionary and Cellular Compatibility evidence require an "
        "optional external data bundle (ProteinInteractionHunter/PIH) that "
        "was not supplied for this run, so neither category was evaluated "
        "for any candidate in this report."
    )


def build_biological_interpretation(
    row: dict[str, Any],
    rank: int,
    n_candidates: int,
    category_refs: Sequence[CategoryRef],
    evolutionary_closer: str,
) -> str:
    """Build the deterministic "Biological Interpretation" paragraph.

    The opening hedge sentence (design spec section 35: never assert this
    candidate IS the target enzyme/interaction partner) is mandatory and
    never varies in structure. ``evolutionary_closer`` is
    build_evolutionary_closer's output for this run, passed in rather than
    recomputed here so a caller writing many candidates can compute it once.
    """
    candidate_id = str(row.get("candidate_protein_id") or "")
    query_id = str(row.get("query_id") or "")
    opening = (
        f"Based on the evidence currently available to this pipeline, "
        f"{candidate_id} ranked {rank} of {n_candidates} evaluated "
        f"candidates for query {query_id}. This reflects the evidence "
        "categories described below as currently available to the "
        f"pipeline — it is not a confirmed identification of "
        f"{candidate_id} as the target enzyme or interaction partner, and "
        "should be read as 'best-supported by currently available "
        "evidence,' not as a settled conclusion."
    )

    color_sentences: list[str] = []
    seen_kinds: set[str] = set()
    for ref in category_refs:
        if ref.kind in seen_kinds or ref.kind not in _CATEGORY_COLOR_SENTENCES:
            continue
        if not _is_contributing(row.get(ref.row_key)):
            continue
        color_sentences.append(_CATEGORY_COLOR_SENTENCES[ref.kind])
        seen_kinds.add(ref.kind)

    closer = (
        "Categories not listed above were not evaluated for this candidate "
        "in this run — a blank category means no evidence was available "
        "to assess, not that the evidence was checked and found absent. "
        f"{evolutionary_closer} {_NEGATIVE_EVIDENCE_CLOSER}"
    )

    return " ".join([opening, *color_sentences, closer])


__all__: tuple[str, ...] = (
    "CategoryRef",
    "build_biological_interpretation",
    "build_evolutionary_closer",
    "build_why_ranks_highly",
)
