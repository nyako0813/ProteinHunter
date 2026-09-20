"""Tests for the Phase 6-8 Stage 2 single-file Word report writer."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest
from docx import Document
from docx.opc.constants import RELATIONSHIP_TYPE as RT

from core.models import DomainHit, ProteinRecord
from config import (
    INTERACTION_ALPHAFOLD_DEFAULT,
    INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    INTERACTION_NEIGHBORHOOD_DEFAULT,
    INTERACTION_SCORING_WEIGHTS_DEFAULT,
    WORD_REPORT_DEFAULT,
    InteractionScoringConfig,
    WordReportConfig,
)
from analysis.interaction_scoring import CONSERVED_QUERY_WARNING_PREFIX
from output.report_v2 import bookmark_name
from output.word_report import category_refs_for_scoring_model, write_word_report


def blast_classification(**buckets) -> SimpleNamespace:
    defaults = {
        "all_records": {},
        "positive_only_records": {},
        "candidates_relaxed_records": {},
        "positive_all_sources_records": {},
        "negative_unmatched_records": {},
        "no_hit_records": {},
        "negative_hit_records": {},
    }
    defaults.update(buckets)
    return SimpleNamespace(**defaults)


def app_config(
    *,
    scoring_model: str = "v2_evidence_based",
    pih_evidence_bundle: str | None = None,
    word_report: WordReportConfig = WORD_REPORT_DEFAULT,
) -> SimpleNamespace:
    return SimpleNamespace(
        interaction_scoring=InteractionScoringConfig(
            enabled=True,
            query_proteins=(),
            query_fasta=None,
            candidate_sources={"candidates": True},
            max_candidates_per_query=200,
            include_sequences_in_excel=False,
            scoring_weights=INTERACTION_SCORING_WEIGHTS_DEFAULT,
            alphafold=INTERACTION_ALPHAFOLD_DEFAULT,
            neighborhood=INTERACTION_NEIGHBORHOOD_DEFAULT,
            scoring_model=scoring_model,
            scoring_engine_config=None,
            functional_complementarity_ruleset=None,
            pih_evidence_bundle=pih_evidence_bundle,
            evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
            word_report=word_report,
        )
    )


def _pair_row(query_id: str, candidate_id: str, candidate_source: str = "Candidates", **extra) -> dict:
    base = {
        "query_id": query_id,
        "candidate_protein_id": candidate_id,
        "candidate_source": candidate_source,
        "candidate_description": "a candidate protein",
        "negative_hit_strength": "none",
        "final_score": 10.0,
        "final_score_tier": "Tier4_Weak",
        "evidence_category_count": 1,
        "alphafold_recommended": False,
    }
    base.update(extra)
    return base


def _interaction_result(source_rows: dict[str, list[dict]], scoring_model: str = "v2_evidence_based") -> SimpleNamespace:
    return SimpleNamespace(
        source_rows=source_rows,
        evidence_detail_rows=[],
        evidence_detail_scoring_model=scoring_model,
        query_rows=[{"query_id": qid} for qid in {row["query_id"] for rows in source_rows.values() for row in rows}],
        neighborhood_rows=[],
    )


def test_write_word_report_creates_file(tmp_path: Path) -> None:
    output_path = tmp_path / "reports" / "report.docx"
    source_rows = {"Interaction_Candidates": [_pair_row("q1", "c1", final_score=80.0, final_score_tier="Tier1_VeryStrong")]}

    result = write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result(source_rows),
        excel_filename="results.xlsx",
    )

    assert result == output_path.resolve()
    assert result.exists()
    Document(str(result))  # must open without raising


def test_word_report_has_expected_top_level_sections(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    source_rows = {"Interaction_Candidates": [_pair_row("q1", "c1")]}

    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result(source_rows),
    )

    document = Document(str(output_path))
    headings = [p.text for p in document.paragraphs if p.style.name.startswith("Heading")]
    assert "5. Evidence Architecture" in headings
    assert "7. Candidate Ranking" in headings
    assert "8. Candidate Details" in headings
    assert "7.1 Query: q1" in headings
    assert "8.1 Query: q1" in headings
    assert any(h.startswith("5.7 Negative Evidence") for h in headings)


def test_word_report_multiple_queries_get_separate_subsections(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    source_rows = {
        "Interaction_Candidates": [
            _pair_row("q1", "c1"),
            _pair_row("q2", "c1"),
        ]
    }

    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result(source_rows),
    )

    document = Document(str(output_path))
    headings = [p.text for p in document.paragraphs if p.style.name.startswith("Heading")]
    assert "7.1 Query: q1" in headings
    assert "7.2 Query: q2" in headings
    assert "8.1 Query: q1" in headings
    assert "8.2 Query: q2" in headings


def test_word_report_settings_force_field_update_on_open(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]}),
    )

    document = Document(str(output_path))
    assert "updateFields" in document.settings.element.xml


def test_word_report_toc_field_present(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]}),
    )

    document = Document(str(output_path))
    assert "TOC" in document.element.xml


def test_word_report_candidate_bookmark_matches_report_v2_bookmark_name(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "cand_1")]}),
    )

    document = Document(str(output_path))
    expected_name = bookmark_name("q1", "cand_1")
    assert f'w:name="{expected_name}"' in document.element.xml


def test_word_report_respects_top_n_and_tier_safety_net(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    # build_workbook_sheets always reranks by final_score (rerank_final_score_rows),
    # so candidate_rank here is assigned from final_score descending: c1=1 .. c5=5.
    rows = [_pair_row("q1", f"c{i}", final_score=100.0 - i) for i in range(1, 6)]
    # rank 5 (c5) is beyond max_per_query=3, but Tier1 -> must still appear.
    rows[4]["final_score_tier"] = "Tier1_VeryStrong"

    word_report = WordReportConfig(enabled=True, max_candidates_per_query=3)
    write_word_report(
        config=app_config(word_report=word_report),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": rows}),
    )

    document = Document(str(output_path))
    headings = [p.text for p in document.paragraphs if p.style.name == "Heading 4"]
    # c1..c3 (top 3 by final_score) plus c5 (Tier1 safety net) = 4 candidates shown, c4 excluded.
    assert len(headings) == 4
    assert not any(h.startswith("c4") for h in headings)
    assert any(h.startswith("c5") for h in headings)


def test_word_report_legacy_additive_model_does_not_crash(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    row = _pair_row(
        "q1",
        "c1",
        candidate_priority_score=20.0,
        same_gene_neighborhood_score=5.0,
        co_occurrence_score=3.0,
        domain_complementarity_score=2.0,
        string_ppi_score=0.0,
    )

    write_word_report(
        config=app_config(scoring_model="legacy_additive"),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [row]}, scoring_model="legacy_additive"),
    )

    document = Document(str(output_path))
    text = "\n".join(p.text for p in document.paragraphs)
    assert "legacy_additive scoring model" in text


def test_word_report_evidence_architecture_states_pih_bundle_absence(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    write_word_report(
        config=app_config(pih_evidence_bundle=None),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]}),
    )

    document = Document(str(output_path))
    text = "\n".join(p.text for p in document.paragraphs)
    assert "was not supplied for this run" in text


def test_word_report_5_7_includes_conserved_query_warning_without_prefix(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    interaction_result = _interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]})
    interaction_result.warnings = [
        f"{CONSERVED_QUERY_WARNING_PREFIX}Query HdrD1 itself has a strong negative-reference BLAST hit.",
        "some unrelated warning",
    ]

    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=interaction_result,
    )

    text = "\n".join(p.text for p in Document(str(output_path)).paragraphs)
    assert "Query HdrD1 itself has a strong negative-reference BLAST hit." in text
    assert CONSERVED_QUERY_WARNING_PREFIX not in text
    assert "some unrelated warning" not in text


def test_word_report_5_7_unchanged_without_conserved_query_warning(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"

    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]}),
    )

    text = "\n".join(p.text for p in Document(str(output_path)).paragraphs)
    assert "itself has a" not in text


def test_word_report_prints_excel_cross_reference_as_plain_text_not_link(tmp_path: Path) -> None:
    output_path = tmp_path / "report.docx"
    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "cand_1")]}),
        excel_filename="ProteinHunter_results.xlsx",
    )

    document = Document(str(output_path))
    text = "\n".join(p.text for p in document.paragraphs)
    assert "ProteinHunter_results.xlsx" in text
    assert "query_id=q1" in text
    assert "candidate_protein_id=cand_1" in text
    hyperlink_targets = [r.target_ref for r in document.part.rels.values() if r.reltype == RT.HYPERLINK]
    assert "ProteinHunter_results.xlsx" not in hyperlink_targets


def test_category_refs_for_v2_include_all_six_categories() -> None:
    from analysis.scoring_engine_config import load_scoring_engine_config

    engine_config = load_scoring_engine_config(None)
    refs = category_refs_for_scoring_model("v2_evidence_based", engine_config, None)
    kinds = {ref.kind for ref in refs}
    assert kinds == {
        "sequence",
        "genomic_context",
        "functional_domain",
        "interaction",
        "evolutionary",
        "cellular_compatibility",
    }


def test_category_refs_for_legacy_use_scoring_weights() -> None:
    refs = category_refs_for_scoring_model("legacy_additive", None, INTERACTION_SCORING_WEIGHTS_DEFAULT)
    caps = {ref.row_key: ref.cap for ref in refs}
    assert caps["candidate_priority_score"] == INTERACTION_SCORING_WEIGHTS_DEFAULT.candidate_priority
    assert caps["co_occurrence_score"] == INTERACTION_SCORING_WEIGHTS_DEFAULT.co_occurrence


# ---------------------------------------------------------------------------
# Run provenance on the title page (core/provenance.py)
# ---------------------------------------------------------------------------


def _make_provenance(*, check_passed: bool | None = None, commit: str | None = "1a2b3c4", dirty: bool | None = False):
    from datetime import datetime

    from core.provenance import RunProvenance

    return RunProvenance("5.0", commit, dirty, "deadbeef", check_passed, datetime(2026, 9, 20, 12, 30))


def _report_text(tmp_path: Path, **kwargs) -> str:
    output_path = tmp_path / "report.docx"
    write_word_report(
        config=app_config(),
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]}),
        **kwargs,
    )
    return "\n".join(p.text for p in Document(str(output_path)).paragraphs)


def test_title_page_shows_code_version_fingerprint_and_sidecar_name(tmp_path: Path) -> None:
    text = _report_text(tmp_path, excel_filename="MA_4115_2.xlsx", provenance=_make_provenance(dirty=True))

    assert "Code version: 5.0 (git 1a2b3c4+dirty)" in text
    assert "Config fingerprint: deadbeef" in text
    assert "saved alongside this report as MA_4115_2.run_provenance.yaml" in text
    assert "Reference genome check" not in text


def test_title_page_without_provenance_is_unchanged(tmp_path: Path) -> None:
    text = _report_text(tmp_path, excel_filename="MA_4115_2.xlsx")

    assert "Code version" not in text
    assert "Config fingerprint" not in text
    assert "run_provenance" not in text
    assert "Scoring model: v2_evidence_based" in text


def test_title_page_warns_only_when_reference_genome_check_failed(tmp_path: Path) -> None:
    failed = _report_text(tmp_path, excel_filename="r.xlsx", provenance=_make_provenance(check_passed=False))
    passed = _report_text(tmp_path, excel_filename="r.xlsx", provenance=_make_provenance(check_passed=True))
    unknown = _report_text(tmp_path, excel_filename="r.xlsx", provenance=_make_provenance(check_passed=None))

    assert "Reference genome check: WARNING" in failed
    assert "Reference genome check" not in passed
    assert "Reference genome check" not in unknown


def test_title_page_handles_missing_excel_filename_and_unknown_git(tmp_path: Path) -> None:
    text = _report_text(tmp_path, provenance=_make_provenance(commit=None, dirty=None))

    assert "Code version: 5.0 (git unknown)" in text
    assert "<Excel workbook stem>.run_provenance.yaml" in text


# ---------------------------------------------------------------------------
# Japanese localization (report_language: ja; output/report_i18n.py)
# ---------------------------------------------------------------------------

import re as _re

from docx.oxml.ns import qn as _qn

from output.word_narrative import bullet_list as _bullet_list, heading as _heading, NarrativeSection as _Section
from output.word_report import JAPANESE_FONT, _render_sections

_JA_CHARS = _re.compile(r"[぀-ヿ㐀-鿿]")


def _localized_document(tmp_path: Path, language: str | None, *, name: str = "report.docx", **kwargs) -> Document:
    config = app_config()
    if language is not None:
        config.report_language = language
    rows = [
        _pair_row("q1", "c1", final_score=80.0, final_score_tier="Tier1_VeryStrong", candidate_rank=1),
        _pair_row("q1", "c2", candidate_source="Negative_hit", negative_hit_strength="medium", candidate_rank=2),
    ]
    output_path = tmp_path / name
    write_word_report(
        config=config,
        blast_classification=blast_classification(),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": rows}),
        excel_filename="results.xlsx",
        **kwargs,
    )
    return Document(str(output_path))


def _bookmark_names(document: Document) -> list[str]:
    return [element.get(_qn("w:name")) for element in document.element.body.iter(_qn("w:bookmarkStart"))]


def test_japanese_report_has_japanese_headings(tmp_path: Path) -> None:
    document = _localized_document(tmp_path, "ja")

    headings = [p.text for p in document.paragraphs if p.style.name.startswith("Heading")]
    for expected in ("5. 証拠の構成", "7. 候補の順位", "8. 候補の詳細", "7.1 クエリ: q1", "8.1 クエリ: q1"):
        assert expected in headings
    assert any(h.startswith("5.7 負の証拠") for h in headings)
    assert "ProteinHunter 候補レポート" in [p.text for p in document.paragraphs]


def test_japanese_report_leaves_no_english_prose(tmp_path: Path) -> None:
    document = _localized_document(tmp_path, "ja")

    for paragraph in document.paragraphs:
        if paragraph.style.name == "Heading 4" or not paragraph.text.strip():
            continue  # candidate titles are ids and gene descriptions
        assert _JA_CHARS.search(paragraph.text), paragraph.text
    for table in document.tables:
        assert all(_JA_CHARS.search(cell.text) for cell in table.rows[0].cells)  # header row is translated


def test_japanese_report_keeps_ids_scores_and_classification_values(tmp_path: Path) -> None:
    en = _localized_document(tmp_path, "en", name="en.docx")
    ja = _localized_document(tmp_path, "ja", name="ja.docx")

    def cells(document: Document) -> list[list[str]]:
        return [[cell.text for cell in row.cells] for table in document.tables for row in table.rows[1:]]

    assert cells(en) == cells(ja)  # candidate ids, "80.0", tier and candidate_source values are identical
    assert cells(ja)[0][3] == "Tier1_VeryStrong"
    ja_text = "\n".join(p.text for p in ja.paragraphs)
    assert "80.0/100" in ja_text and "c1" in ja_text and "Negative_hit(陰性ヒット)" in ja_text


def test_bookmarks_for_the_excel_links_do_not_depend_on_the_language(tmp_path: Path) -> None:
    en = _localized_document(tmp_path, "en", name="en.docx")
    ja = _localized_document(tmp_path, "ja", name="ja.docx")

    assert _bookmark_names(en) == _bookmark_names(ja) != []


def test_english_is_the_default_and_identical_to_an_explicit_en(tmp_path: Path) -> None:
    default = _localized_document(tmp_path, None, name="default.docx")
    explicit = _localized_document(tmp_path, "en", name="explicit.docx")

    def structure(document: Document) -> list[str]:
        return [p.text for p in document.paragraphs if not p.text.startswith("Report generated")] + [
            cell.text for table in document.tables for row in table.rows for cell in row.cells
        ]

    assert structure(default) == structure(explicit)
    assert "Candidate Details" in " ".join(structure(default))


def test_japanese_output_sets_east_asian_fonts_and_language_only_in_japanese(tmp_path: Path) -> None:
    ja = _localized_document(tmp_path, "ja", name="ja.docx")
    en = _localized_document(tmp_path, "en", name="en.docx")

    for style_name in ("Normal", "Heading 1", "Title"):
        fonts = ja.styles[style_name].element.rPr.find(_qn("w:rFonts"))
        assert fonts.get(_qn("w:eastAsia")) == JAPANESE_FONT
        assert _qn("w:eastAsiaTheme") not in fonts.attrib
    lang = ja.styles["Normal"].element.rPr.find(_qn("w:lang"))
    assert lang.get(_qn("w:eastAsia")) == "ja-JP"

    normal_fonts = en.styles["Normal"].element.rPr.find(_qn("w:rFonts")) if en.styles["Normal"].element.rPr is not None else None
    assert normal_fonts is None or normal_fonts.get(_qn("w:eastAsia")) != JAPANESE_FONT


def test_japanese_title_page_provenance_and_genome_warning(tmp_path: Path) -> None:
    document = _localized_document(tmp_path, "ja", provenance=_make_provenance(check_passed=False))

    text = "\n".join(p.text for p in document.paragraphs)
    assert "コードのバージョン: 5.0 (git 1a2b3c4)" in text
    assert "設定フィンガープリント: deadbeef" in text
    assert "results.run_provenance.yaml" in text
    assert "参照ゲノムのチェック: 警告" in text


def test_conserved_query_note_is_localized_in_the_japanese_report(tmp_path: Path) -> None:
    note = (
        "conserved query visibility: Query WP_1 itself has a strong negative-reference BLAST hit "
        "(negative_hit_strength=strong). Interaction partners of broadly conserved proteins may be conserved too."
    )
    config = app_config()
    config.report_language = "ja"
    interaction_result = _interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]})
    interaction_result.warnings = [note]
    write_word_report(config=config, blast_classification=blast_classification(), output_path=tmp_path / "r.docx", interaction_result=interaction_result)

    text = "\n".join(p.text for p in Document(str(tmp_path / "r.docx")).paragraphs)
    assert "クエリ WP_1 自身が" in text and "strong(強)" in text


def test_unsupported_report_language_is_rejected(tmp_path: Path) -> None:
    config = app_config()
    config.report_language = "fr"

    with pytest.raises(ValueError, match="unsupported report language"):
        write_word_report(config=config, blast_classification=blast_classification(), output_path=tmp_path / "r.docx", interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]}))


def test_renderer_draws_bullet_lists_and_rejects_unknown_kinds() -> None:
    document = Document()

    _render_sections(document, [_heading("Title", 1), _bullet_list(["first", "second"])])

    assert [p.text for p in document.paragraphs] == ["Title", "first", "second"]
    assert [p.style.name for p in document.paragraphs][1:] == ["List Bullet", "List Bullet"]
    with pytest.raises(ValueError, match="unknown report section kind"):
        _render_sections(document, [_Section("carousel")])


# ---------------------------------------------------------------------------
# Domain information in the candidate details (shared with the Notion export)
# ---------------------------------------------------------------------------


def test_candidate_details_list_the_domains_and_the_pipelines_domain_evidence(tmp_path: Path) -> None:
    hit = DomainHit(source="CDD", accession="cd01", name="HUP domain", description="ATP pyrophosphatase", start=5, end=120, evalue=1e-20)
    classification = blast_classification(all_records={"c1": ProteinRecord(protein_id="c1", sequence="MKV", description="d", domains=[hit])})
    interaction_result = _interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1"), _pair_row("q1", "c2")]})
    interaction_result.evidence_detail_rows = [
        {"query_id": "q1", "candidate_protein_id": "c1", "component_name": "domain_complementarity", "status": "AVAILABLE",
         "explanation": "domain family match: a x b", "normalized_value": 1.0}
    ]
    output_path = tmp_path / "report.docx"

    write_word_report(config=app_config(), blast_classification=classification, output_path=output_path, interaction_result=interaction_result, excel_filename="r.xlsx")

    document = Document(str(output_path))
    texts = [p.text for p in document.paragraphs]
    assert any(t.startswith("Domain information: 1 domain hit(s)") for t in texts)
    bullets = [p.text for p in document.paragraphs if p.style.name == "List Bullet"]
    assert bullets == ["CDD cd01 — HUP domain: ATP pyrophosphatase [aa 5–120, E=1.0e-20]"]
    assert any(t == "Domain evidence used for scoring: domain family match: a x b (domain_complementarity = 1.00)" for t in texts)
    assert sum(t.startswith("Domain information") for t in texts) == 1  # c2 has no annotation record: no block


def test_japanese_candidate_details_have_the_domain_block_too(tmp_path: Path) -> None:
    hit = DomainHit(source="Pfam", accession="PF01", name="X")
    config = app_config()
    config.report_language = "ja"
    output_path = tmp_path / "report.docx"

    write_word_report(
        config=config,
        blast_classification=blast_classification(all_records={"c1": ProteinRecord(protein_id="c1", sequence="MKV", description="d", domains=[hit])}),
        output_path=output_path,
        interaction_result=_interaction_result({"Interaction_Candidates": [_pair_row("q1", "c1")]}),
    )

    texts = [p.text for p in Document(str(output_path)).paragraphs]
    assert any(t.startswith("ドメイン情報: この候補には 1 件のドメインヒット") for t in texts)
    assert "Pfam PF01 — X" in [p.text for p in Document(str(output_path)).paragraphs if p.style.name == "List Bullet"]
