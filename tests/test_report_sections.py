"""Tests for the renderer-neutral report structure (output/report_sections.py)."""

from __future__ import annotations

import re
from datetime import datetime

import pytest

from core.provenance import RunProvenance
from output.report_sections import build_report_sections, localize_conserved_query_note
from output.report_v2 import bookmark_name
from output.word_narrative import CategoryRef

_JAPANESE = re.compile(r"[぀-ヿ㐀-鿿]")
NOW = datetime(2026, 9, 20, 12, 30)
REFS = (
    CategoryRef("candidate_priority_score", "sequence", "Sequence/Source Classification", 30.0),
    CategoryRef("same_gene_neighborhood_score", "genomic_context", "Genomic Context", 25.0),
    CategoryRef("interaction_evidence_score", "interaction", "Interaction", 47.0),
)
CAPS = {"source_classification": 30.0, "genomic_context": 25.0, "functional_annotation": 20.0, "external_ppi_evidence": 15.0}


def row(query: str, candidate: str, rank: int, tier: str = "Tier2_Strong", source: str = "Candidates", **extra) -> dict:
    base = {
        "query_id": query,
        "candidate_protein_id": candidate,
        "candidate_description": "a hypothetical enzyme",
        "candidate_rank": rank,
        "candidate_source": source,
        "negative_hit_strength": "none",
        "final_score": 61.25,
        "final_score_tier": tier,
        "evidence_category_count": 3,
        "candidate_priority_score": 30.0,
        "same_gene_neighborhood_score": 12.5,
        "interaction_evidence_score": 9.0,
    }
    base.update(extra)
    return base


GROUPED = [
    ("MA_4115", [row("MA_4115", "MA_0363", 1), row("MA_4115", "MA_4110", 2, tier="Tier3_Moderate", source="Negative_hit", negative_hit_strength="medium")]),
    ("MA_0688", [row("MA_0688", "MA_0687", 1, tier="Tier1_VeryStrong", source="No_hit")]),
]


def build(language: str, **overrides):
    kwargs = dict(
        language=language,
        generated_at=NOW,
        scoring_model="v2_evidence_based",
        category_caps=CAPS,
        pih_bundle_configured=False,
        grouped=GROUPED,
        category_refs=REFS,
        max_per_query=15,
        excel_filename="MA_4115_2.xlsx",
        provenance=None,
        conserved_query_notes=(),
    )
    kwargs.update(overrides)
    return build_report_sections(**kwargs)


def test_structure_is_identical_in_both_languages() -> None:
    """Only wording differs: same sections, same order, same levels, anchors and table data."""
    en, ja = build("en"), build("ja")

    assert [(s.kind, s.level, s.anchor) for s in en] == [(s.kind, s.level, s.anchor) for s in ja]
    assert [s.rows for s in en if s.kind == "table"] == [s.rows for s in ja if s.kind == "table"]
    assert [s.header for s in en if s.kind == "table"] != [s.header for s in ja if s.kind == "table"]


def test_every_candidate_gets_one_bookmark_anchor_named_like_the_excel_link_expects() -> None:
    sections = build("ja")

    anchors = [s.anchor for s in sections if s.kind == "heading" and s.anchor]
    assert anchors == [bookmark_name(q, r["candidate_protein_id"]) for q, rows in GROUPED for r in rows]


def test_default_language_is_english_and_matches_explicit_en() -> None:
    assert build_report_sections(
        generated_at=NOW, scoring_model="v2_evidence_based", category_caps=CAPS, pih_bundle_configured=False,
        grouped=GROUPED, category_refs=REFS, max_per_query=15, excel_filename="x.xlsx",
    ) == build("en", excel_filename="x.xlsx")


def test_english_title_page_wording() -> None:
    texts = [s.text for s in build("en")[:5]]

    assert texts[0] == "ProteinHunter Candidate Report"
    assert texts[1] == "Report generated: 2026-09-20 12:30"
    assert texts[2] == "Scoring model: v2_evidence_based"
    assert texts[3] == (
        "Queries evaluated: 2. Candidates shown per query: up to 15, plus any additional "
        "Tier1_VeryStrong/Tier2_Strong candidate regardless of rank."
    )


def test_japanese_report_has_no_untranslated_prose() -> None:
    """Every prose paragraph and heading (except candidate titles, which are ids + gene descriptions) contains Japanese."""
    for section in build("ja"):
        if section.kind == "heading" and section.level == 4:
            continue
        for text in (section.text, section.label):
            if text:
                assert _JAPANESE.search(text), (section.kind, text)


def test_data_values_are_not_translated() -> None:
    ja = build("ja")

    header, body = next(s for s in ja if s.kind == "table").header, next(s for s in ja if s.kind == "table").rows
    assert body[0] == ("1", "MA_0363", "61.2", "Tier2_Strong", "Candidates")
    assert body[1][3:] == ("Tier3_Moderate", "Negative_hit")
    text = " ".join(s.text for s in ja)
    assert "MA_0363" in text and "MA_4115" in text and "61.2/100" in text
    assert "Tier2_Strong" not in header  # headers are labels and are translated


def test_narrative_glosses_labels_as_english_then_japanese() -> None:
    text = " ".join(s.text for s in build("ja"))

    assert "Genomic Context(ゲノム文脈) (12.5/25)" in text  # category label, then the unchanged score
    assert "medium(中)" in text  # negative_hit_strength keeps its own spelling
    assert "Tier 1 — Very Strong(非常に強い)" in text
    assert "Negative_hit(陰性ヒット)" in text


def test_scores_and_ids_in_narrative_are_identical_across_languages() -> None:
    def numbers(sections) -> list[str]:
        return sorted(re.findall(r"\d+\.\d+/\d+|MA_\d+", " ".join(s.text for s in sections)))

    # Same numbers and locus tags appear, whichever language the prose is in.
    assert numbers(build("en")) == numbers(build("ja"))


def test_empty_report_is_localized() -> None:
    en = build("en", grouped=[])
    ja = build("ja", grouped=[])

    assert [s.text for s in en if s.text == "No query-specific candidates were produced by this run."] != []
    assert [s.text for s in ja if s.text == "この実行では、クエリ別の候補は得られなかった。"] != []
    assert [s.kind for s in en] == [s.kind for s in ja]


def _provenance(check_passed: bool | None) -> RunProvenance:
    return RunProvenance("5.0", "1a2b3c4", True, "deadbeef", check_passed, NOW)


def test_provenance_lines_and_genome_warning_follow_the_language() -> None:
    en = " | ".join(s.text for s in build("en", provenance=_provenance(False)))
    ja = " | ".join(s.text for s in build("ja", provenance=_provenance(False)))

    assert "Code version: 5.0 (git 1a2b3c4+dirty)" in en and "Reference genome check: WARNING" in en
    assert "コードのバージョン: 5.0 (git 1a2b3c4+dirty)" in ja and "参照ゲノムのチェック: 警告" in ja
    assert "MA_4115_2.run_provenance.yaml" in en and "MA_4115_2.run_provenance.yaml" in ja
    assert "Reference genome check" not in " | ".join(s.text for s in build("en", provenance=_provenance(True)))


CONSERVED_EN = (
    "Query WP_011024006.1 itself has a strong negative-reference BLAST hit "
    "(negative_hit_strength=strong). Interaction partners of broadly conserved proteins may be conserved too."
)


def test_conserved_query_note_is_translated_from_its_parts_and_english_passes_through() -> None:
    assert localize_conserved_query_note(CONSERVED_EN, "en") == CONSERVED_EN

    ja = localize_conserved_query_note(CONSERVED_EN, "ja")
    assert "WP_011024006.1" in ja and "strong(強)" in ja and "negative_hit_strength=strong" in ja
    assert "candidate_sources.negative_hit" in ja and _JAPANESE.search(ja)


def test_conserved_query_note_in_an_unexpected_shape_is_kept_verbatim() -> None:
    odd = "Something the analysis module says differently in a future version."

    assert localize_conserved_query_note(odd, "ja") == odd


def test_notes_appear_after_the_negative_evidence_section_only() -> None:
    sections = build("ja", conserved_query_notes=[CONSERVED_EN])

    texts = [s.text for s in sections]
    negative_heading = next(i for i, s in enumerate(sections) if s.kind == "heading" and s.text.startswith("5.7"))
    note_index = next(i for i, text in enumerate(texts) if "WP_011024006.1" in text)
    ranking_index = next(i for i, s in enumerate(sections) if s.kind == "heading" and s.text.startswith("7."))
    assert negative_heading < note_index < ranking_index
