"""Phase 6-8 Stage 2: single-file Word report generation.

Builds one ``.docx`` per pipeline run (design spec section 29: one Word
file per run, not per query -- multiple queries become subsections within
the report's fixed "7. Candidate Ranking" / "8. Candidate Details"
sections instead, see
``claude/phase678_stage2_word_report_investigation.md`` item 1).

This module owns everything ``python-docx``-specific: low-level OOXML
helpers for bookmarks, external hyperlinks, and the Table of Contents
field (``python-docx`` has no high-level API for any of the three), plus
rendering. What the report *says* -- structure and every wording, in English
or Japanese (``report_language``) -- is built by ``output/report_sections.py``
as a renderer-neutral list of ``NarrativeSection``; this module only draws
that list into a document. It reuses ``output/report_v2.py`` for row
shaping/selection (the same consolidated rows the Excel workbook is built
from -- see ``build_workbook_sheets``). It has no knowledge of Excel/openpyxl -- the reverse direction
(Excel's ``word_report_link`` column) lives in ``output/excel.py`` and
depends on this module's ``bookmark_name``/``select_top_candidates_per_query``
re-exports from ``output/report_v2.py``, not on this module directly, so
the two writers stay independently callable.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from docx import Document
from docx.enum.style import WD_STYLE_TYPE
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import RGBColor
from docx.text.paragraph import Paragraph

from analysis.interaction_scoring import CONSERVED_QUERY_WARNING_PREFIX, is_conserved_query_visibility_warning
from analysis.scoring_engine_config import ScoringEngineConfig, load_scoring_engine_config
from core.exceptions import WordReportError
from core.provenance import RunProvenance
from output.report_i18n import DEFAULT_LANGUAGE, normalize_language
from output.report_sections import build_report_sections
from output.report_v2 import (
    TIER_SAFETY_NET,
    build_workbook_sheets,
    select_top_candidates_per_query,
)
from output.word_narrative import CategoryRef, NarrativeSection


# ---------------------------------------------------------------------------
# Low-level OOXML helpers -- python-docx has no high-level API for any of
# these three (bookmarks, external hyperlinks, TOC field), so this is the
# well-documented ~15-line-per-helper escape hatch into raw OOXML that
# claude/phase678_excel_word_redesign_investigation.md (item 3) anticipated.
# All confirmed working (create -> save -> reopen -> structurally verify)
# during this module's own M1 smoke test.
# ---------------------------------------------------------------------------


def _ensure_hyperlink_style(document: Document) -> None:
    """Add a "Hyperlink" character style if the document's template lacks one.

    python-docx's default template does not ship a "Hyperlink" style (only
    Heading 1-9 and a handful of others) -- without this, a hyperlink run
    styled "Hyperlink" would reference a style that does not exist, which
    Word tolerates but renders as plain unstyled text (still clickable,
    just visually indistinguishable from a normal run).
    """
    if "Hyperlink" in document.styles:
        return
    style = document.styles.add_style("Hyperlink", WD_STYLE_TYPE.CHARACTER)
    style.font.color.rgb = RGBColor(0x05, 0x63, 0xC1)
    style.font.underline = True


def _set_update_fields_on_open(document: Document) -> None:
    """Force Word to recompute fields (the TOC in particular) when the file is opened.

    Without this, python-docx's TOC field shows only its placeholder text
    ("Right-click and choose Update Field...") until a person manually
    updates it once -- the user explicitly asked to avoid requiring that
    manual step. ``<w:updateFields w:val="true"/>`` in settings.xml is the
    standard OOXML mechanism for "recalculate fields on open"; confirmed to
    round-trip correctly (present after save + reopen) during this
    module's M1 smoke test.
    """
    settings = document.settings.element
    update_fields = OxmlElement("w:updateFields")
    update_fields.set(qn("w:val"), "true")
    settings.insert(0, update_fields)


def _add_bookmark(paragraph: Paragraph, name: str, bookmark_id: int) -> None:
    """Wrap ``paragraph``'s existing content in a named bookmark.

    Word bookmark ids only need to be unique within one document; callers
    pass a running counter. ``name`` must already be Word-legal (see
    ``output/report_v2.py::bookmark_name``) -- not re-validated here.
    """
    start = OxmlElement("w:bookmarkStart")
    start.set(qn("w:id"), str(bookmark_id))
    start.set(qn("w:name"), name)
    end = OxmlElement("w:bookmarkEnd")
    end.set(qn("w:id"), str(bookmark_id))
    paragraph._p.insert(0, start)
    paragraph._p.append(end)


def _add_external_hyperlink(paragraph: Paragraph, url: str, text: str) -> None:
    """Append a clickable external hyperlink run to ``paragraph``.

    Not currently called by this module (the Word->Excel reference below
    is deliberately plain text, not a link -- see the module docstring's
    cross-reference to option D in
    claude/phase678_excel_word_redesign_investigation.md item 6), but kept
    here as the validated M1 primitive for a future direction that does
    want a real Word-side external link.
    """
    part = paragraph.part
    r_id = part.relate_to(url, RT.HYPERLINK, is_external=True)
    hyperlink = OxmlElement("w:hyperlink")
    hyperlink.set(qn("r:id"), r_id)
    run = OxmlElement("w:r")
    run_properties = OxmlElement("w:rPr")
    style = OxmlElement("w:rStyle")
    style.set(qn("w:val"), "Hyperlink")
    run_properties.append(style)
    run.append(run_properties)
    text_element = OxmlElement("w:t")
    text_element.text = text
    run.append(text_element)
    hyperlink.append(run)
    paragraph._p.append(hyperlink)


def _add_toc_field(document: Document, placeholder: str) -> None:
    """Insert a Word Table of Contents field covering heading levels 1-3.

    ``placeholder`` is the text shown until the field is updated.

    Candidate-detail headings are deliberately level 4 (see
    ``_write_candidate_details``) so the TOC lists report sections and
    per-query subsections only, not every individual candidate.
    """
    paragraph = document.add_paragraph()
    begin_run = paragraph.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = 'TOC \\o "1-3" \\h \\z \\u'
    begin_run._r.append(fld_begin)
    begin_run._r.append(instr)

    separate_run = paragraph.add_run()
    fld_separate = OxmlElement("w:fldChar")
    fld_separate.set(qn("w:fldCharType"), "separate")
    separate_run._r.append(fld_separate)

    placeholder_run = OxmlElement("w:r")
    placeholder_text = OxmlElement("w:t")
    placeholder_text.text = placeholder
    placeholder_run.append(placeholder_text)
    paragraph._p.append(placeholder_run)

    end_run = paragraph.add_run()
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    end_run._r.append(fld_end)


# ---------------------------------------------------------------------------
# Category references (which row columns feed "why ranks highly" /
# "Biological Interpretation", and what they are called/capped at) --
# resolved once per run from the run's actual config, per scoring_model.
# ---------------------------------------------------------------------------


def category_refs_for_scoring_model(
    scoring_model: str,
    engine_config: ScoringEngineConfig,
    legacy_weights: Any,
) -> tuple[CategoryRef, ...]:
    """Resolve the CategoryRef list output/word_narrative.py should enumerate for one run.

    v2_evidence_based uses the scoring engine's own category_caps (the
    same numbers 05_Sequence_Evidence..10_Negative_Evidence categorize
    Interaction_Evidence_Detail rows by). legacy_additive has no category
    concept at all -- its closest analogues are the fixed
    interaction_scoring.scoring_weights point budget, kept as separate
    line items (e.g. co_occurrence and domain_complementarity are not
    combined the way v2's functional_domain_score already is) since that
    is how legacy_additive actually computes and exposes them.
    """
    if scoring_model == "v2_evidence_based":
        caps = engine_config.category_caps
        interaction_cap = (
            caps.get("external_ppi_evidence", 0.0)
            + caps.get("coexpression_evidence", 0.0)
            + caps.get("pih_direct_interaction", 0.0)
        )
        return (
            CategoryRef("candidate_priority_score", "sequence", "Sequence/Source Classification", caps.get("source_classification", 0.0)),
            CategoryRef("same_gene_neighborhood_score", "genomic_context", "Genomic Context", caps.get("genomic_context", 0.0)),
            CategoryRef("functional_domain_score", "functional_domain", "Functional/Domain", caps.get("functional_annotation", 0.0)),
            CategoryRef("interaction_evidence_score", "interaction", "Interaction", interaction_cap),
            CategoryRef("evolutionary_score", "evolutionary", "Evolutionary", caps.get("pih_evolutionary", 0.0)),
            CategoryRef("cellular_compatibility_score", "cellular_compatibility", "Cellular Compatibility", caps.get("pih_cellular_compatibility", 0.0)),
        )

    weights = legacy_weights
    return (
        CategoryRef("candidate_priority_score", "sequence", "Sequence/Source Classification", getattr(weights, "candidate_priority", 0.0)),
        CategoryRef("same_gene_neighborhood_score", "genomic_context", "Genomic Context", getattr(weights, "gene_neighborhood", 0.0)),
        CategoryRef("co_occurrence_score", "functional_domain", "Co-occurrence", getattr(weights, "co_occurrence", 0.0)),
        CategoryRef("domain_complementarity_score", "functional_domain", "Domain Complementarity", getattr(weights, "domain_complementarity", 0.0)),
        CategoryRef("string_ppi_score", "interaction", "Interaction (STRING PPI)", getattr(weights, "external_ppi", 0.0)),
    )


# ---------------------------------------------------------------------------
# Sections 7 / 8: per-query candidate ranking and details
# ---------------------------------------------------------------------------


def _group_by_query(rows: list[dict[str, Any]]) -> list[tuple[str, list[dict[str, Any]]]]:
    """Group already query_id-then-rank-sorted rows into (query_id, rows) pairs.

    Preserves the incoming order (see report_v2.rerank_final_score_rows)
    rather than re-sorting -- query section order is a deliberate,
    reproducible property of the upstream sort, not decided here.
    """
    groups: list[tuple[str, list[dict[str, Any]]]] = []
    current_query: str | None = None
    current_rows: list[dict[str, Any]] = []
    for row in rows:
        query_id = str(row.get("query_id") or "")
        if query_id != current_query:
            if current_query is not None:
                groups.append((current_query, current_rows))
            current_query = query_id
            current_rows = []
        current_rows.append(row)
    if current_query is not None:
        groups.append((current_query, current_rows))
    return groups


# ---------------------------------------------------------------------------
# Rendering: NarrativeSection list -> python-docx
# ---------------------------------------------------------------------------

#: East Asian font applied to the document's styles for Japanese output. Without
#: an explicit ``eastAsia`` font Word falls back to whatever the system default
#: is (often a serif face); Yu Gothic ships with current Windows and macOS Office.
JAPANESE_FONT = "Yu Gothic"
_JAPANESE_STYLE_NAMES = ("Normal", "Title", "Heading 1", "Heading 2", "Heading 3", "Heading 4", "Table Grid")


def _apply_japanese_fonts(document: Document) -> None:
    """Point the document's main styles at a Japanese font and mark the language as ja-JP."""
    for name in _JAPANESE_STYLE_NAMES:
        if name not in document.styles:
            continue
        run_properties = document.styles[name].element.get_or_add_rPr()
        fonts = run_properties.find(qn("w:rFonts"))
        if fonts is None:
            fonts = OxmlElement("w:rFonts")
            run_properties.insert(0, fonts)
        fonts.set(qn("w:eastAsia"), JAPANESE_FONT)
        theme_attribute = qn("w:eastAsiaTheme")
        if theme_attribute in fonts.attrib:  # a theme font would override the explicit one
            del fonts.attrib[theme_attribute]
        if name == "Normal":
            language = OxmlElement("w:lang")
            language.set(qn("w:eastAsia"), "ja-JP")
            run_properties.append(language)


def _render_sections(document: Document, sections: list[NarrativeSection]) -> None:
    """Draw the report sections into ``document`` in order.

    Headings that carry an ``anchor`` become Word bookmarks (the targets of the
    Excel workbook's ``word_report_link`` column); ids run 1, 2, ... in order.
    """
    bookmark_id = 1
    for section in sections:
        if section.kind == "heading":
            heading = document.add_heading(section.text, level=section.level)
            if section.anchor:
                _add_bookmark(heading, section.anchor, bookmark_id)
                bookmark_id += 1
        elif section.kind == "paragraph":
            paragraph = document.add_paragraph()
            if section.label:
                paragraph.add_run(section.label).bold = True
            paragraph.add_run(section.text)
        elif section.kind == "bullet_list":
            for item in section.items:
                document.add_paragraph(item, style="List Bullet")
        elif section.kind == "table":
            grid = document.add_table(rows=1, cols=len(section.header))
            grid.style = "Table Grid"
            for cell, text in zip(grid.rows[0].cells, section.header):
                cell.text = text
            for row in section.rows:
                for cell, text in zip(grid.add_row().cells, row):
                    cell.text = text
        elif section.kind == "toc":
            _add_toc_field(document, section.text)
        else:
            raise ValueError(f"unknown report section kind: {section.kind!r}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def write_word_report(
    config: Any,
    blast_classification: Any,
    output_path: str | Path,
    interaction_result: Any | None = None,
    excel_filename: str = "",
    provenance: RunProvenance | None = None,
) -> Path:
    """Write the Phase 6-8 Stage 2 single-file Word report and return its path.

    Mirrors ``output/excel.write_classification_workbook``'s signature
    (same ``config``/``blast_classification``/``interaction_result`` the
    Excel writer already receives at main.py's call site) so both writers
    can be called independently from the same pipeline state.
    ``excel_filename`` is printed as a plain-text cross-reference in each
    candidate's detail section (design spec's Excel cross-link
    requirement, option D -- see
    claude/phase678_excel_word_redesign_investigation.md item 6: reliable
    Word->Excel deep-linking to a specific row is not available, so this
    is a stable filename+id reference rather than a fragile link).

    ``provenance`` (run provenance, core/provenance.py) adds the code
    version, config fingerprint and sidecar filename to the title page (plus
    a warning when the reference-genome check failed). Left ``None`` (the
    default), the title page is unchanged.

    The report's language comes from ``config.report_language`` ("en", the
    default, or "ja"): one setting decides the language of the one Word file.
    Only wording is translated; locus tags, scores and classification values
    stay as in the Excel workbook (output/report_i18n.py).
    """
    language = normalize_language(getattr(config, "report_language", DEFAULT_LANGUAGE))
    resolved_output = Path(output_path).expanduser().resolve()
    resolved_output.parent.mkdir(parents=True, exist_ok=True)

    scoring_config = getattr(config, "interaction_scoring", None)
    scoring_model = str(getattr(scoring_config, "scoring_model", "legacy_additive"))
    engine_config = load_scoring_engine_config(getattr(scoring_config, "scoring_engine_config", None))
    legacy_weights = getattr(scoring_config, "scoring_weights", None)
    pih_bundle_configured = bool(getattr(scoring_config, "pih_evidence_bundle", None))
    word_report_config = getattr(scoring_config, "word_report", None)
    max_per_query = int(getattr(word_report_config, "max_candidates_per_query", 15))

    try:
        sheets_data = build_workbook_sheets(config, blast_classification, interaction_result)
        selected_rows = select_top_candidates_per_query(
            sheets_data["final_score_rows"], max_per_query, TIER_SAFETY_NET
        )
        grouped = _group_by_query(selected_rows)
        category_refs = category_refs_for_scoring_model(scoring_model, engine_config, legacy_weights)
        conserved_query_notes = [
            warning[len(CONSERVED_QUERY_WARNING_PREFIX):]
            for warning in getattr(interaction_result, "warnings", None) or []
            if is_conserved_query_visibility_warning(warning)
        ]
        sections = build_report_sections(
            language=language,
            generated_at=datetime.now(),
            scoring_model=scoring_model,
            category_caps=engine_config.category_caps,
            pih_bundle_configured=pih_bundle_configured,
            grouped=grouped,
            category_refs=category_refs,
            max_per_query=max_per_query,
            excel_filename=excel_filename,
            provenance=provenance,
            conserved_query_notes=conserved_query_notes,
        )

        document = Document()
        _ensure_hyperlink_style(document)
        _set_update_fields_on_open(document)
        if language == "ja":
            _apply_japanese_fonts(document)
        _render_sections(document, sections)

        document.save(str(resolved_output))
    except Exception as exc:
        message = (
            f"ProteinHunter could not write the Word report: {resolved_output}. "
            "Please check that the folder is writable and the file is not open."
        )
        raise WordReportError(message) from exc

    return resolved_output


__all__: tuple[str, ...] = (
    "WordReportError",
    "category_refs_for_scoring_model",
    "write_word_report",
)
