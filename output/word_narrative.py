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

from dataclasses import dataclass
from typing import Any, NamedTuple, Sequence

from output.report_i18n import DEFAULT_LANGUAGE, STRINGS, bilingual, normalize_language, t


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


# ---------------------------------------------------------------------------
# Intermediate representation shared by every report renderer
# ---------------------------------------------------------------------------
# output/report_sections.py assembles the whole report as a list of these;
# output/word_report.py renders them to .docx, and the planned Notion export
# (claude/notion_export_design.md) renders the same list to Notion blocks.
# Nothing here knows about python-docx or Notion.
#
# Kinds and fields (decided while implementing the Japanese localization):
#   heading      text, level (0 = document title, 1-4), anchor (optional stable id;
#                docx renders it as a bookmark the Excel workbook links to)
#   paragraph    text, label (optional bold lead-in, rendered before the text)
#   bullet_list  items
#   table        header, rows (first row of cells is the header; all text)
#   toc          text (placeholder shown until the reader updates the field);
#                the docx renderer draws a TOC field, the Notion renderer its
#                native table-of-contents block
#
# Structural hints for renderers that split the report (Notion turns each
# candidate into its own database page); the docx renderer ignores them:
#   role         "details" on the "8. Candidate Details" heading, "query" on a
#                per-query heading inside it, "candidate" on a candidate title
#   meta         key/value pairs of the data behind a "query"/"candidate"
#                heading (query_id, candidate_id, rank, final_score, tier,
#                source, description) -- data values, never translated


@dataclass(frozen=True)
class NarrativeSection:
    kind: str
    text: str = ""
    level: int = 0
    label: str = ""
    items: tuple[str, ...] = ()
    header: tuple[str, ...] = ()
    rows: tuple[tuple[str, ...], ...] = ()
    anchor: str = ""
    role: str = ""
    meta: tuple[tuple[str, str], ...] = ()

    def meta_dict(self) -> dict[str, str]:
        return dict(self.meta)


def heading(
    text: str,
    level: int,
    anchor: str = "",
    role: str = "",
    meta: Sequence[tuple[str, str]] = (),
) -> NarrativeSection:
    return NarrativeSection("heading", text=text, level=level, anchor=anchor, role=role, meta=tuple(meta))


def paragraph(text: str, label: str = "") -> NarrativeSection:
    return NarrativeSection("paragraph", text=text, label=label)


def bullet_list(items: Sequence[str]) -> NarrativeSection:
    return NarrativeSection("bullet_list", items=tuple(items))


def table(header: Sequence[str], rows: Sequence[Sequence[str]]) -> NarrativeSection:
    return NarrativeSection("table", header=tuple(header), rows=tuple(tuple(row) for row in rows))


def toc(text: str) -> NarrativeSection:
    return NarrativeSection("toc", text=text)


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


def build_why_ranks_highly(
    row: dict[str, Any],
    category_refs: Sequence[CategoryRef],
    language: str = DEFAULT_LANGUAGE,
) -> str:
    """Build the deterministic "why this candidate ranks highly" paragraph.

    ``category_refs`` should already be resolved for this row's
    ``scoring_model`` (see output/word_report.py::category_refs_for_scoring_model)
    -- this function only reads ``row`` and the refs, no model branching.
    ``language`` picks the wording (output/report_i18n.py); in Japanese,
    classification labels appear as ``English(日本語)`` and data values such as
    scores are unchanged.
    """
    language = normalize_language(language)
    tier = row.get("final_score_tier")
    tier_key = f"tier_opening.{tier}"
    if tier is not None and tier_key in STRINGS[language]:
        opening = t(
            tier_key,
            language,
            final_score=row.get("final_score") or 0.0,
            evidence_category_count=row.get("evidence_category_count") or 0,
        )
    else:
        opening = t("tier_opening.unclassified", language)

    contributing = [ref for ref in category_refs if _is_contributing(row.get(ref.row_key))]
    if contributing:
        listing = t("why.listing_separator", language).join(
            f"{bilingual(ref.label, language)} ({row[ref.row_key]:.1f}/{ref.cap:.0f})" for ref in contributing
        )
        evidence_sentence = t("why.contributing", language, listing=listing)
    else:
        evidence_sentence = t("why.none", language)

    candidate_source = str(row.get("candidate_source") or "")
    negative_hit_strength = str(row.get("negative_hit_strength") or "none")
    strength_text = bilingual(negative_hit_strength, language)
    source_key = f"source.{candidate_source}"
    source_sentence = (
        t(source_key, language, negative_hit_strength=strength_text) if source_key in STRINGS[language] else ""
    )

    parts = [opening, evidence_sentence]
    if source_sentence:
        parts.append(source_sentence)
    if negative_hit_strength != "none" and candidate_source != "Negative_hit":
        parts.append(t("negative_hit_caveat", language, negative_hit_strength=strength_text))
    return t("sentence_joiner", language).join(parts)


def build_evolutionary_closer(
    scoring_model: str,
    pih_bundle_configured: bool,
    language: str = DEFAULT_LANGUAGE,
) -> str:
    """Return the run-level sentence explaining Evolutionary/Cellular Compatibility coverage.

    Constant for every candidate in one run (depends only on
    ``scoring_model``/``pih_bundle_configured``, never on the row) -- callers
    may compute it once per run and reuse it rather than recomputing per
    candidate, though calling it repeatedly is also correct (it is a pure
    function).
    """
    if scoring_model != "v2_evidence_based":
        return t("evolutionary_closer.legacy", language)
    if pih_bundle_configured:
        return t("evolutionary_closer.pih", language)
    return t("evolutionary_closer.no_pih", language)


def build_biological_interpretation(
    row: dict[str, Any],
    rank: int,
    n_candidates: int,
    category_refs: Sequence[CategoryRef],
    evolutionary_closer: str,
    language: str = DEFAULT_LANGUAGE,
) -> str:
    """Build the deterministic "Biological Interpretation" paragraph.

    The opening hedge sentence (design spec section 35: never assert this
    candidate IS the target enzyme/interaction partner) is mandatory and
    never varies in structure, in either language. ``evolutionary_closer`` is
    build_evolutionary_closer's output for this run, passed in rather than
    recomputed here so a caller writing many candidates can compute it once.
    """
    language = normalize_language(language)
    candidate_id = str(row.get("candidate_protein_id") or "")
    query_id = str(row.get("query_id") or "")
    opening = t(
        "interpretation.opening",
        language,
        candidate_id=candidate_id,
        rank=rank,
        n_candidates=n_candidates,
        query_id=query_id,
    )

    color_sentences: list[str] = []
    seen_kinds: set[str] = set()
    for ref in category_refs:
        color_key = f"color.{ref.kind}"
        if ref.kind in seen_kinds or color_key not in STRINGS[language]:
            continue
        if not _is_contributing(row.get(ref.row_key)):
            continue
        color_sentences.append(t(color_key, language))
        seen_kinds.add(ref.kind)

    closer = t(
        "interpretation.closer",
        language,
        evolutionary_closer=evolutionary_closer,
        negative_evidence_closer=t("negative_evidence_closer", language),
    )

    return t("sentence_joiner", language).join([opening, *color_sentences, closer])


__all__: tuple[str, ...] = (
    "CategoryRef",
    "NarrativeSection",
    "bullet_list",
    "build_biological_interpretation",
    "build_evolutionary_closer",
    "build_why_ranks_highly",
    "heading",
    "paragraph",
    "table",
    "toc",
)
