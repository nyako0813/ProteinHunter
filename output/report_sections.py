"""The whole candidate report as a renderer-neutral list of ``NarrativeSection``.

This is the "what does the report say" half of the Word report: title page,
"5. Evidence Architecture", "7. Candidate Ranking" and "8. Candidate Details",
with every wording taken from ``output/report_i18n.py`` in the requested
language and the per-candidate sentences from ``output/word_narrative.py``.
It has no knowledge of ``python-docx``: ``output/word_report.py`` renders the
list to a ``.docx`` and the planned Notion export renders the same list to
Notion blocks (claude/notion_export_design.md), so wording and structure are
decided exactly once.

Data values are never translated: locus tags, protein ids, scores, tier names
and ``candidate_source`` values in tables come straight from the rows, so the
report keeps matching the Excel workbook whatever the language.
"""

from __future__ import annotations

import re
from datetime import datetime
from typing import Any, Mapping, Sequence

from core.provenance import RunProvenance, code_version_text, provenance_sidecar_filename
from output.report_i18n import DEFAULT_LANGUAGE, bilingual, normalize_language, t
from output.report_v2 import bookmark_name
from output.word_narrative import (
    CategoryRef,
    NarrativeSection,
    build_biological_interpretation,
    build_evolutionary_closer,
    build_why_ranks_highly,
    heading,
    paragraph,
    table,
    toc,
)

_CONSERVED_NOTE = re.compile(
    r"^Query (?P<query_id>\S+) itself has a (?P<strength>\w+) negative-reference BLAST hit "
    r"\(negative_hit_strength=(?P<strength_raw>\w+)\)\."
)


def localize_conserved_query_note(note: str, language: str = DEFAULT_LANGUAGE) -> str:
    """The conserved-query visibility note (analysis/interaction_scoring.py) in ``language``.

    The analysis module words that warning in English for the run log, and
    the Word report reuses it. For Japanese the note is rebuilt from its query
    id and strength; if the message ever stops matching the expected shape the
    English original is kept rather than guessed at.
    """
    if normalize_language(language) == "en":
        return note
    match = _CONSERVED_NOTE.match(note)
    if match is None:
        return note
    return t(
        "conserved_query.note",
        language,
        query_id=match["query_id"],
        strength=bilingual(match["strength"], language),
        strength_raw=match["strength_raw"],
    )


def _title_sections(
    language: str,
    generated_at: datetime,
    scoring_model: str,
    n_queries: int,
    max_per_query: int,
    excel_filename: str,
    provenance: RunProvenance | None,
) -> list[NarrativeSection]:
    sections = [
        heading(t("report.title", language), 0),
        paragraph(t("report.generated", language, timestamp=f"{generated_at:%Y-%m-%d %H:%M}")),
        paragraph(t("report.scoring_model", language, scoring_model=scoring_model)),
    ]
    if provenance is not None:
        location = (
            provenance_sidecar_filename(excel_filename)
            if excel_filename
            else t("report.sidecar_unknown_location", language)
        )
        sections += [
            paragraph(t("report.code_version", language, version=code_version_text(provenance))),
            paragraph(t("report.config_fingerprint", language, fingerprint=provenance.config_hash)),
            paragraph(t("report.sidecar", language, location=location)),
        ]
        if provenance.reference_genome_check_passed is False:
            sections.append(
                paragraph(
                    t(
                        "report.reference_genome_check",
                        language,
                        text=t("report.reference_genome_warning", language),
                    )
                )
            )
    sections.append(paragraph(t("report.queries_evaluated", language, n_queries=n_queries, max_per_query=max_per_query)))
    sections.append(toc(t("report.toc_placeholder", language)))
    return sections


def _architecture_sections(
    language: str,
    scoring_model: str,
    caps: Mapping[str, float],
    pih_bundle_configured: bool,
    conserved_query_notes: Sequence[str],
) -> list[NarrativeSection]:
    closer = build_evolutionary_closer(scoring_model, pih_bundle_configured, language)
    interaction_cap = (
        caps.get("external_ppi_evidence", 0.0)
        + caps.get("coexpression_evidence", 0.0)
        + caps.get("pih_direct_interaction", 0.0)
    )
    sections = [
        heading(t("arch.title", language), 1),
        # repr() keeps the quotes the report has always shown around the model name.
        paragraph(t("arch.intro", language, scoring_model=repr(scoring_model))),
        heading(t("arch.seq.heading", language, cap=caps.get("source_classification", 0.0)), 2),
        paragraph(t("arch.seq.body", language)),
        heading(t("arch.func.heading", language, cap=caps.get("functional_annotation", 0.0)), 2),
        paragraph(t("arch.func.body", language)),
        heading(t("arch.genomic.heading", language, cap=caps.get("genomic_context", 0.0)), 2),
        paragraph(t("arch.genomic.body", language)),
        heading(t("arch.interaction.heading", language, cap=interaction_cap), 2),
        paragraph(t("arch.interaction.body", language)),
        heading(t("arch.evolutionary.heading", language, cap=caps.get("pih_evolutionary", 0.0)), 2),
        paragraph(t("arch.evolutionary.body", language, closer=closer)),
        heading(t("arch.cellular.heading", language, cap=caps.get("pih_cellular_compatibility", 0.0)), 2),
        paragraph(t("arch.cellular.body", language, closer=closer)),
        heading(t("arch.negative.heading", language), 2),
        paragraph(t("arch.negative.body", language)),
    ]
    # A known side effect of that absence, not a Negative Evidence signal:
    # see patches/conserved_query_visibility_design.md.
    sections += [paragraph(localize_conserved_query_note(note, language)) for note in conserved_query_notes]
    return sections


def _ranking_sections(language: str, grouped: Sequence[tuple[str, Sequence[dict[str, Any]]]]) -> list[NarrativeSection]:
    sections = [heading(t("ranking.title", language), 1)]
    if not grouped:
        sections.append(paragraph(t("ranking.empty", language)))
        return sections

    header = (
        t("ranking.col.rank", language),
        t("ranking.col.candidate", language),
        t("ranking.col.final_score", language),
        t("ranking.col.tier", language),
        t("ranking.col.source", language),
    )
    for index, (query_id, rows) in enumerate(grouped, start=1):
        sections.append(heading(t("ranking.query_heading", language, index=index, query_id=query_id), 2))
        body = []
        for row in rows:
            final_score = row.get("final_score")
            body.append(
                (
                    str(row.get("candidate_rank") or ""),
                    str(row.get("candidate_protein_id") or ""),
                    f"{final_score:.1f}" if final_score is not None else "—",
                    str(row.get("final_score_tier") or "—"),
                    str(row.get("candidate_source") or ""),
                )
            )
        sections.append(table(header, body))
    return sections


def _details_sections(
    language: str,
    grouped: Sequence[tuple[str, Sequence[dict[str, Any]]]],
    category_refs: Sequence[CategoryRef],
    evolutionary_closer: str,
    excel_filename: str,
) -> list[NarrativeSection]:
    sections = [heading(t("details.title", language), 1)]
    if not grouped:
        sections.append(paragraph(t("details.empty", language)))
        return sections

    for index, (query_id, rows) in enumerate(grouped, start=1):
        sections.append(heading(t("details.query_heading", language, index=index, query_id=query_id), 2))
        n_candidates = len(rows)
        for row in rows:
            candidate_id = str(row.get("candidate_protein_id") or "")
            description = str(row.get("candidate_description") or "").strip()
            title = candidate_id if not description else f"{candidate_id} — {description[:80]}"
            sections.append(heading(title, 4, anchor=bookmark_name(query_id, candidate_id)))
            sections.append(
                paragraph(
                    build_why_ranks_highly(row, category_refs, language),
                    label=t("details.why_label", language),
                )
            )
            sections.append(
                paragraph(
                    build_biological_interpretation(
                        row,
                        rank=int(row.get("candidate_rank") or 0),
                        n_candidates=n_candidates,
                        category_refs=category_refs,
                        evolutionary_closer=evolutionary_closer,
                        language=language,
                    ),
                    label=t("details.interpretation_label", language),
                )
            )
            sections.append(
                paragraph(
                    t(
                        "details.excel_ref",
                        language,
                        excel_filename=excel_filename,
                        query_id=query_id,
                        candidate_id=candidate_id,
                    )
                )
            )
    return sections


def build_report_sections(
    *,
    language: str = DEFAULT_LANGUAGE,
    generated_at: datetime,
    scoring_model: str,
    category_caps: Mapping[str, float],
    pih_bundle_configured: bool,
    grouped: Sequence[tuple[str, Sequence[dict[str, Any]]]],
    category_refs: Sequence[CategoryRef],
    max_per_query: int,
    excel_filename: str = "",
    provenance: RunProvenance | None = None,
    conserved_query_notes: Sequence[str] = (),
) -> list[NarrativeSection]:
    """Every section of the report, in reading order.

    ``grouped`` is the already-selected ``(query_id, rows)`` list (see
    output/word_report._group_by_query); ``category_refs`` and
    ``category_caps`` are resolved once per run by the caller from the run's
    scoring model and engine config.
    """
    language = normalize_language(language)
    evolutionary_closer = build_evolutionary_closer(scoring_model, pih_bundle_configured, language)
    return [
        *_title_sections(language, generated_at, scoring_model, len(grouped), max_per_query, excel_filename, provenance),
        *_architecture_sections(language, scoring_model, category_caps, pih_bundle_configured, conserved_query_notes),
        *_ranking_sections(language, grouped),
        *_details_sections(language, grouped, category_refs, evolutionary_closer, excel_filename),
    ]


__all__: tuple[str, ...] = (
    "build_report_sections",
    "localize_conserved_query_note",
)
