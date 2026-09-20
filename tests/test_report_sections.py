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


def test_structural_hints_mark_the_details_query_and_candidate_headings() -> None:
    sections = build("en")

    details = [s for s in sections if s.role == "details"]
    queries = [s for s in sections if s.role == "query"]
    candidates = [s for s in sections if s.role == "candidate"]
    assert len(details) == 1 and details[0].level == 1
    assert [q.meta_dict() for q in queries] == [{"query_id": "MA_4115"}, {"query_id": "MA_0688"}]
    assert [c.meta_dict()["candidate_id"] for c in candidates] == ["MA_0363", "MA_4110", "MA_0687"]
    assert candidates[0].meta_dict() == {
        "query_id": "MA_4115",
        "candidate_id": "MA_0363",
        "rank": "1",
        "final_score": "61.2",
        "tier": "Tier2_Strong",
        "source": "Candidates",
        "description": "a hypothetical enzyme",
    }


def test_structural_hints_are_data_values_identical_in_both_languages() -> None:
    def hints(language: str):
        return [(s.role, s.meta) for s in build(language) if s.role]

    assert hints("en") == hints("ja")



# ---------------------------------------------------------------------------
# Evidence-section roles and the per-candidate domain information
# ---------------------------------------------------------------------------

from types import SimpleNamespace

from output.report_sections import DomainEntry, DomainEvidence, collect_domain_data


def test_evidence_headings_are_marked_for_renderers_that_split_them() -> None:
    sections = build("en")

    roots = [s for s in sections if s.role == "evidence_root"]
    evidence = [s for s in sections if s.role == "evidence"]
    assert len(roots) == 1 and roots[0].level == 1
    assert [s.meta_dict()["section"] for s in evidence] == [
        "sequence", "functional_domain", "genomic_context", "interaction", "evolutionary", "cellular_compatibility", "negative",
    ]
    assert all(s.level == 2 for s in evidence)


def _candidate_block(sections, candidate_id: str):
    """The sections between a candidate's heading and the next heading."""
    start = next(i for i, s in enumerate(sections) if s.role == "candidate" and s.meta_dict()["candidate_id"] == candidate_id)
    end = next((i for i in range(start + 1, len(sections)) if sections[i].kind == "heading"), len(sections))
    return sections[start + 1 : end]


def test_domain_block_lists_hits_in_sequence_order_between_interpretation_and_excel_reference() -> None:
    domains = {"MA_0363": [
        DomainEntry("CDD", "cd02", "second", "", 200, 300, 2e-5),
        DomainEntry("Pfam", "PF01", "first", "the first domain", 10, 150, 1e-30),
        DomainEntry("CDD", "cdX", "no-position"),
    ]}
    evidence = {("MA_4115", "MA_0363"): DomainEvidence("domain family match (UniProt Pfam/InterPro): a x b", 1.0)}

    block = _candidate_block(build("en", domains_by_protein=domains, domain_evidence=evidence), "MA_0363")

    assert [s.kind for s in block] == ["paragraph", "paragraph", "paragraph", "bullet_list", "paragraph", "paragraph"]
    assert block[2].label == "Domain information: " and "3 domain hit(s)" in block[2].text
    assert block[3].items == (
        "Pfam PF01 — first: the first domain [aa 10–150, E=1.0e-30]",
        "CDD cd02 — second [aa 200–300, E=2.0e-05]",
        "CDD cdX — no-position",
    )
    assert block[4].label == "Domain evidence used for scoring: "
    assert block[4].text == "domain family match (UniProt Pfam/InterPro): a x b (domain_complementarity = 1.00)"
    assert block[5].text.startswith("Full data for this candidate")


def test_long_domain_lists_are_capped_with_a_count_of_the_rest() -> None:
    many = [DomainEntry("CDD", f"cd{i:03d}", f"d{i}", "", i, i + 5, 1e-3) for i in range(30)]

    block = _candidate_block(build("en", domains_by_protein={"MA_0363": many}), "MA_0363")

    items = next(s for s in block if s.kind == "bullet_list").items
    assert len(items) == 26 and items[-1] == "...and 5 more domain hit(s) not listed here."


def test_domain_block_distinguishes_no_hits_from_no_record() -> None:
    sections = build("en", domains_by_protein={"MA_0363": []})

    no_hits = _candidate_block(sections, "MA_0363")
    unknown = _candidate_block(sections, "MA_4110")
    assert "No domain hits are recorded" in no_hits[2].text and no_hits[2].label == "Domain information: "
    assert [s.kind for s in unknown] == ["paragraph"] * 3  # no annotation record: no domain block at all


def test_without_domain_data_the_report_is_exactly_as_before() -> None:
    plain = build("en")

    assert plain == build("en", domains_by_protein=None, domain_evidence=None)
    assert all("Domain information" not in s.label for s in plain)


def test_domain_block_is_localized_but_carries_the_same_data() -> None:
    domains = {"MA_0363": [DomainEntry("CDD", "cd01", "HUP", "ATP pyrophosphatase", 5, 120, 1e-20)]}
    evidence = {("MA_4115", "MA_0363"): DomainEvidence("domain family match: a x b", 0.5)}

    en = _candidate_block(build("en", domains_by_protein=domains, domain_evidence=evidence), "MA_0363")
    ja = _candidate_block(build("ja", domains_by_protein=domains, domain_evidence=evidence), "MA_0363")

    assert ja[2].label == "ドメイン情報: " and "1 件のドメインヒット" in ja[2].text
    assert ja[3].items == en[3].items  # data lines are not translated
    assert ja[4].label == "スコアリングで用いたドメイン証拠: " and ja[4].text == en[4].text


def test_collect_domain_data_reads_records_and_the_domain_complementarity_evidence() -> None:
    hit = SimpleNamespace(source="CDD", accession="cd01", name="HUP", description="d", start=1, end=9, evalue=1e-5)
    classification = SimpleNamespace(all_records={"MA_0363": SimpleNamespace(domains=[hit]), "MA_4110": SimpleNamespace(domains=[])})
    detail_rows = [
        {"query_id": "MA_4115", "candidate_protein_id": "MA_0363", "component_name": "co_occurrence", "status": "AVAILABLE", "explanation": "overlap", "normalized_value": 1.0},
        {"query_id": "MA_4115", "candidate_protein_id": "MA_0363", "component_name": "domain_complementarity", "status": "MISSING", "explanation": "missing"},
        {"query_id": "MA_4115", "candidate_protein_id": "MA_0363", "component_name": "domain_complementarity", "status": "AVAILABLE", "explanation": "family match", "normalized_value": 0.75},
        {"query_id": "MA_4115", "candidate_protein_id": "OTHER", "component_name": "domain_complementarity", "status": "AVAILABLE", "explanation": "not on the report", "normalized_value": 1.0},
    ]

    domains, evidence = collect_domain_data(classification, SimpleNamespace(evidence_detail_rows=detail_rows), GROUPED)

    assert domains == {"MA_0363": [DomainEntry("CDD", "cd01", "HUP", "d", 1, 9, 1e-5)], "MA_4110": []}  # MA_0687 has no record
    assert evidence == {("MA_4115", "MA_0363"): DomainEvidence("family match", 0.75)}


def test_collect_domain_data_tolerates_missing_attributes() -> None:
    assert collect_domain_data(SimpleNamespace(), None, GROUPED) == ({}, {})


def test_placeholder_descriptions_are_not_printed() -> None:
    block = _candidate_block(
        build("en", domains_by_protein={"MA_0363": [DomainEntry("CDD", "cl00959", "Nitrate_red_gam superfamily", "-", 96, 254, 1.9e-4)]}),
        "MA_0363",
    )

    assert next(s for s in block if s.kind == "bullet_list").items == ("CDD cl00959 — Nitrate_red_gam superfamily [aa 96–254, E=1.9e-04]",)

