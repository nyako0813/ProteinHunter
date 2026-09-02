"""Minimal Excel output helpers for ProteinHunter candidate records."""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl.styles import Alignment, Font
from openpyxl.utils import get_column_letter
from openpyxl.worksheet.worksheet import Worksheet

from core.exceptions import ExcelOutputError
from core.models import BlastHit, ProteinRecord
from analysis.interaction_scoring import (
    INTERACTION_EVIDENCE_DETAIL_SHEET,
    INTERACTION_NEIGHBORHOOD_SHEET,
    INTERACTION_QUERY_COLUMNS,
    interaction_evidence_detail_columns,
    interaction_index_rows,
    interaction_neighborhood_columns,
    interaction_pair_columns,
)


EXCEL_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "description",
    "old_locus_tag",
    "total_score",
    "score_components",
    "score_reasons",
    "domain_sources",
    "domain_names",
    "domain_accessions",
    "domain_descriptions",
    "domain_count",
    "unique_domain_count",
    "unique_domain_accessions",
    "unique_domain_names",
    "sequence_length",
    "positive_hit_count",
    "negative_hit_count",
    "blast_status",
    "best_positive_hit",
    "best_positive_bitscore",
    "best_positive_evalue",
    "best_negative_hit",
    "best_negative_bitscore",
    "best_negative_evalue",
    "negative_best_identity",
    "negative_best_query_coverage",
    "negative_best_evalue",
    "negative_best_source",
    "negative_hit_strength",
    "negative_strong_hit_count",
    "negative_medium_hit_count",
    "negative_weak_hit_count",
    "negative_exclusion_reason",
    "motifs",
    "uniprot_accession",
    "alphafold_url",
    "notes",
    "positive_source_count",
    "positive_sources_hit",
    "positive_sources_missing",
)


POSITIVE_SOURCE_SUMMARY_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "description",
    "old_locus_tag",
    "negative_hit",
    "positive_source_count",
    "positive_sources_hit",
    "positive_sources_missing",
    "positive_hit_count",
    "negative_hit_count",
    "blast_status",
    "best_positive_hit",
    "best_positive_bitscore",
    "best_positive_evalue",
    "best_negative_hit",
    "best_negative_bitscore",
    "best_negative_evalue",
)


INDEX_ROWS: tuple[tuple[str, str, str, str], ...] = (
    (
        "Candidates",
        "positive hit present and no negative hit",
        "strict candidate set; negative-free positive-associated targets",
        "first-pass high-confidence candidate review",
    ),
    (
        "Candidates_relaxed",
        "positive hit present and no strong negative hit",
        "retains targets with only medium/weak negative hits",
        "avoid over-filtering by weak homolog/domain-level matches",
    ),
    (
        "Positive_all_sources",
        "hits all positive sources and no negative hit",
        "broadly conserved among cnm5U-positive organisms",
        "search for shared cnm5U-related factors",
    ),
    (
        "Positive_source_summary",
        "summarizes positive source hit distribution for each target",
        "shows how widely each target is conserved in positive sources",
        "compare positive source breadth",
    ),
    (
        "Negative_unmatched",
        "no negative hit",
        "targets not found in cnm5U-negative organisms",
        "review all negative-unmatched targets",
    ),
    (
        "No_hit",
        "no positive hit and no negative hit",
        "Methanosarcina acetivorans-specific or poorly conserved candidates",
        "search for novel thioamidation-related factors",
    ),
    (
        "Negative_hit",
        "any negative hit",
        "targets with at least one hit in cnm5U-negative organisms",
        "inspect what would be excluded by strict filtering",
    ),
    (
        "Negative_strong_hit",
        "at least one strong negative hit",
        "likely common/ortholog-like factor present in negative organisms",
        "generally lower priority or exclusion-oriented candidates",
    ),
    (
        "Negative_medium_hit",
        "medium negative hit but no strong negative hit",
        "ambiguous homolog candidates",
        "check conserved motifs, domains, and structure before excluding",
    ),
    (
        "Negative_weak_hit",
        "weak negative hit only",
        "possible distant homolog or shared domain",
        "do not exclude automatically; use as caution flag",
    ),
)


NEGATIVE_EVIDENCE_EXPLANATIONS: tuple[tuple[str, str], ...] = (
    ("negative_best_identity", "best identity among negative hits"),
    (
        "negative_best_query_coverage",
        "query coverage of the representative/best negative hit",
    ),
    ("negative_best_evalue", "E-value of the representative/best negative hit"),
    (
        "negative_best_source",
        "negative source where the representative/best negative hit was found",
    ),
    ("negative_hit_strength", "strong / medium / weak / none"),
    (
        "negative_strong_hit_count",
        "number of negative hits classified as strong",
    ),
    (
        "negative_medium_hit_count",
        "number of negative hits classified as medium",
    ),
    ("negative_weak_hit_count", "number of negative hits classified as weak"),
    (
        "negative_exclusion_reason",
        "reason why the record was retained or excluded under the current ortholog_filter mode",
    ),
)


INTERACTION_SCORE_EXPLANATIONS: tuple[tuple[str, str], ...] = (
    (
        "interaction_priority_score",
        "Overall interaction/functional priority score including candidate source, gene neighborhood, co-occurrence, domain complementarity, and modeling readiness.",
    ),
    (
        "candidate_priority_score",
        "Score based on which candidate source sheet the protein came from, such as Candidates, Candidates_relaxed, No_hit, etc.",
    ),
    (
        "same_gene_neighborhood_score",
        "Score based on genomic distance from the query protein using GFF coordinates. Close genes get positive score, but distant genes are not excluded.",
    ),
    (
        "distance_bp",
        "Genomic distance in base pairs between query and candidate genes when both coordinates are available.",
    ),
    (
        "co_occurrence_score",
        "Score based on similarity of positive/negative source hit patterns between query and candidate.",
    ),
    (
        "domain_complementarity_score",
        "Score based on meaningful functional terms from description/Pfam/CDD/domain annotations. Generic words alone are ignored.",
    ),
    (
        "alphafold_readiness_score",
        "Score indicating whether the pair is practical to model structurally based on sequence availability and length. This is not evidence of interaction.",
    ),
    (
        "distance_independent_score",
        "Score excluding gene neighborhood and AlphaFold readiness. Formula: candidate_priority_score + co_occurrence_score + domain_complementarity_score.",
    ),
    (
        "distance_independent_rank",
        "Rank based on distance_independent_score within each Interaction_* sheet.",
    ),
    (
        "priority_group",
        "Category such as nearby_candidate, distant_cooccurrence_candidate, distant_domain_candidate, no_hit_candidate, or general_candidate.",
    ),
    (
        "interaction_score_reasons",
        "Human-readable reasons explaining why the candidate received its score.",
    ),
    (
        "protein_hunter_score",
        "Query-independent 'is this generally a good candidate' score for the candidate "
        "protein alone (same formula as the Candidates sheet's total_score, but computed "
        "for every interaction_scoring candidate source, not just Candidates). Does not "
        "change if the query changes, and is never part of interaction_priority_score, "
        "candidate_rank, or sheet sort order -- reference only.",
    ),
    (
        "interaction_score",
        "Query-specific evidence only, re-normalized to 0-100: genomic_context "
        "(+ STRING's neighborhood channel) + domain_complementarity + external_ppi_evidence "
        "(STRING's cooccurrence channel). Populated for both scoring models -- v2_evidence_based "
        "computes it from EvidenceComponent categories (blank/None means no query-specific "
        "evidence was available at all, not a scored zero); legacy_additive computes an "
        "equivalent blended sum (always a number, 0 when there is no evidence, since legacy has "
        "no 'missing vs. evaluated-zero' concept). Deliberately excludes source_classification, "
        "sequence_evidence, and co_occurrence, which mainly reflect the candidate's own "
        "conservation profile rather than evidence specific to this query pair. Reference only -- "
        "does not affect interaction_priority_score, candidate_rank, or sheet sort order.",
    ),
    (
        "string_ppi_score",
        "legacy_additive only (blank for v2_evidence_based, which reports STRING evidence as "
        "separate string_cooccurrence/string_neighborhood rows in Interaction_Evidence_Detail "
        "instead). Averages STRING's cooccurrence and neighborhood channels for this pair into "
        "one 0-weights.external_ppi point score. 0 when interaction_scoring.string_ppi_ncbi_taxon_id "
        "is unset or STRING has no data for this pair -- reference only, folded into both "
        "interaction_priority_score and interaction_score.",
    ),
)

INTERACTION_SCORE_NOTES: tuple[str, ...] = (
    "interaction_priority_score is not a direct protein-protein interaction probability.",
    "Gene neighborhood is used as positive evidence only.",
    "Distant archaeal candidates should not be excluded solely because they are far from the query gene.",
    "AlphaFold readiness does not mean AlphaFold predicts interaction.",
    "string_cooccurrence/string_neighborhood/string_ppi_score are derived from STRING "
    "(https://string-db.org), used under the Creative Commons Attribution 4.0 "
    "International (CC BY 4.0) license. Please credit STRING if these values are "
    "published or redistributed.",
)


def records_to_dataframe(records: dict[str, ProteinRecord]) -> pd.DataFrame:
    """Convert ProteinRecord objects into a tabular pandas DataFrame."""
    rows = [_record_to_row(record) for record in records.values()]
    return pd.DataFrame(rows, columns=EXCEL_COLUMNS)


def write_records_to_excel(
    records: dict[str, ProteinRecord],
    output_path: str | Path,
    sheet_name: str = "Candidates",
) -> Path:
    """Write ProteinRecord objects to an Excel workbook and return its path."""
    resolved_output = Path(output_path).expanduser().resolve()
    resolved_output.parent.mkdir(parents=True, exist_ok=True)
    dataframe = records_to_dataframe(records)

    try:
        with pd.ExcelWriter(resolved_output, engine="openpyxl") as writer:
            dataframe.to_excel(writer, sheet_name=sheet_name, index=False)
            worksheet = writer.sheets[sheet_name]
            _format_worksheet(worksheet, dataframe)
    except Exception as exc:
        message = (
            f"ProteinHunter could not write the Excel file: {resolved_output}. "
            "Please check that the folder is writable and the file is not open."
        )
        raise ExcelOutputError(message) from exc

    return resolved_output


def write_classification_workbook(
    candidates: dict[str, ProteinRecord],
    output_path: str | Path,
    negative_unmatched: dict[str, ProteinRecord] | None = None,
    no_hit: dict[str, ProteinRecord] | None = None,
    negative_hit: dict[str, ProteinRecord] | None = None,
    positive_all_sources: dict[str, ProteinRecord] | None = None,
    positive_source_summary: dict[str, ProteinRecord] | None = None,
    candidates_relaxed: dict[str, ProteinRecord] | None = None,
    negative_strong_hit: dict[str, ProteinRecord] | None = None,
    negative_medium_hit: dict[str, ProteinRecord] | None = None,
    negative_weak_hit: dict[str, ProteinRecord] | None = None,
    interaction_result: Any | None = None,
    include_interaction_sequences: bool = False,
) -> Path:
    """Write candidate records and BLAST classification sheets to Excel."""
    resolved_output = Path(output_path).expanduser().resolve()
    resolved_output.parent.mkdir(parents=True, exist_ok=True)
    sheets = {
        "Candidates": records_to_dataframe(candidates),
        "Candidates_relaxed": records_to_dataframe(candidates_relaxed or {}),
        "Positive_all_sources": records_to_dataframe(positive_all_sources or {}),
        "Positive_source_summary": positive_source_summary_dataframe(
            positive_source_summary or {}
        ),
        "Negative_unmatched": records_to_dataframe(negative_unmatched or {}),
        "No_hit": records_to_dataframe(no_hit or {}),
        "Negative_hit": records_to_dataframe(negative_hit or {}),
        "Negative_strong_hit": records_to_dataframe(negative_strong_hit or {}),
        "Negative_medium_hit": records_to_dataframe(negative_medium_hit or {}),
        "Negative_weak_hit": records_to_dataframe(negative_weak_hit or {}),
    }
    interaction_sheets = _interaction_dataframes(
        interaction_result,
        include_sequences=include_interaction_sequences,
    )
    all_index_rows = (
        *INDEX_ROWS,
        *interaction_index_rows(list(interaction_sheets)),
    )

    try:
        with pd.ExcelWriter(resolved_output, engine="openpyxl") as writer:
            index_dataframe = _index_dataframe(all_index_rows)
            index_dataframe.to_excel(writer, sheet_name="Index", index=False)
            _format_index_worksheet(
                writer.sheets["Index"],
                index_dataframe,
                all_index_rows,
            )

            for sheet_name, dataframe in {**sheets, **interaction_sheets}.items():
                dataframe.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False,
                    startrow=1,
                )
                worksheet = writer.sheets[sheet_name]
                _add_back_to_index_link(worksheet)
                _format_worksheet(worksheet, dataframe, header_row=2)
    except Exception as exc:
        message = (
            f"ProteinHunter could not write the Excel file: {resolved_output}. "
            "Please check that the folder is writable and the file is not open."
        )
        raise ExcelOutputError(message) from exc

    return resolved_output


def positive_source_summary_dataframe(
    records: dict[str, ProteinRecord],
) -> pd.DataFrame:
    """Return the compact positive-source summary DataFrame."""
    rows = [_positive_source_summary_row(record) for record in records.values()]
    return pd.DataFrame(rows, columns=POSITIVE_SOURCE_SUMMARY_COLUMNS)


def _interaction_dataframes(
    interaction_result: Any | None,
    include_sequences: bool,
) -> dict[str, pd.DataFrame]:
    """Return Interaction_* DataFrames only when interaction scoring ran."""
    if interaction_result is None:
        return {}

    sheets: dict[str, pd.DataFrame] = {}
    sheets["Interaction_query"] = pd.DataFrame(
        interaction_result.query_rows,
        columns=INTERACTION_QUERY_COLUMNS,
    )
    pair_columns = interaction_pair_columns(include_sequences)
    for sheet_name, rows in interaction_result.source_rows.items():
        sheets[sheet_name] = pd.DataFrame(rows, columns=pair_columns)

    evidence_detail_rows = getattr(interaction_result, "evidence_detail_rows", [])
    if evidence_detail_rows:
        evidence_detail_scoring_model = getattr(
            interaction_result, "evidence_detail_scoring_model", "legacy_additive"
        )
        sheets[INTERACTION_EVIDENCE_DETAIL_SHEET] = pd.DataFrame(
            evidence_detail_rows,
            columns=interaction_evidence_detail_columns(evidence_detail_scoring_model),
        )

    neighborhood_rows = getattr(interaction_result, "neighborhood_rows", [])
    if neighborhood_rows:
        sheets[INTERACTION_NEIGHBORHOOD_SHEET] = pd.DataFrame(
            neighborhood_rows,
            columns=interaction_neighborhood_columns(),
        )

    return sheets


def _format_worksheet(
    worksheet: Worksheet,
    dataframe: pd.DataFrame,
    header_row: int = 1,
) -> None:
    """Apply simple readability formatting to an Excel worksheet."""
    worksheet.freeze_panes = f"A{header_row + 1}"
    if len(dataframe.columns) > 0:
        last_column = get_column_letter(len(dataframe.columns))
        last_row = max(header_row, header_row + len(dataframe.index))
        worksheet.auto_filter.ref = f"A{header_row}:{last_column}{last_row}"

    for header_cell in worksheet[header_row]:
        header_cell.font = Font(bold=True)

    wrap_columns = {
        "description",
        "domain_descriptions",
        "score_reasons",
        "notes",
        "positive_sources_hit",
        "positive_sources_missing",
        "negative_exclusion_reason",
        "explanation",
        "Selection rule",
        "Biological interpretation",
        "Recommended use",
    }

    for column_index, column_name in enumerate(dataframe.columns, start=1):
        column_letter = get_column_letter(column_index)
        width = _column_width(column_name, dataframe[column_name])

        if column_name in {"description", "domain_descriptions", "score_reasons"}:
            width = max(width, 35)
        elif column_name == "notes":
            width = max(width, 45)
        elif column_name in {
            "Selection rule",
            "Biological interpretation",
            "Recommended use",
        }:
            width = min(max(width, 28), 70)

        worksheet.column_dimensions[column_letter].width = width

        if column_name in wrap_columns:
            for cell in worksheet[column_letter]:
                cell.alignment = Alignment(wrap_text=True, vertical="top")


def _index_dataframe(
    index_rows: tuple[tuple[str, str, str, str], ...],
) -> pd.DataFrame:
    """Return the workbook navigation index DataFrame."""
    rows: list[dict[str, str]] = [
        {
            "Sheet": sheet,
            "Selection rule": selection_rule,
            "Biological interpretation": interpretation,
            "Recommended use": recommended_use,
        }
        for sheet, selection_rule, interpretation, recommended_use in index_rows
    ]
    rows.append(
        {
            "Sheet": "",
            "Selection rule": "",
            "Biological interpretation": "",
            "Recommended use": "",
        }
    )
    rows.append(
        {
            "Sheet": "Interaction scoring columns",
            "Selection rule": "",
            "Biological interpretation": "",
            "Recommended use": "",
        }
    )
    for column_name, explanation in INTERACTION_SCORE_EXPLANATIONS:
        rows.append(
            {
                "Sheet": column_name,
                "Selection rule": explanation,
                "Biological interpretation": "",
                "Recommended use": "",
            }
        )
    rows.append(
        {
            "Sheet": "Interaction scoring notes",
            "Selection rule": "; ".join(INTERACTION_SCORE_NOTES),
            "Biological interpretation": "",
            "Recommended use": "",
        }
    )
    rows.append(
        {
            "Sheet": "",
            "Selection rule": "",
            "Biological interpretation": "",
            "Recommended use": "",
        }
    )
    rows.append(
        {
            "Sheet": "Negative evidence columns",
            "Selection rule": "",
            "Biological interpretation": "",
            "Recommended use": "",
        }
    )
    for column_name, explanation in NEGATIVE_EVIDENCE_EXPLANATIONS:
        rows.append(
            {
                "Sheet": column_name,
                "Selection rule": explanation,
                "Biological interpretation": "",
                "Recommended use": "",
            }
        )

    return pd.DataFrame(
        rows,
        columns=(
            "Sheet",
            "Selection rule",
            "Biological interpretation",
            "Recommended use",
        ),
    )


def _format_index_worksheet(
    worksheet: Worksheet,
    dataframe: pd.DataFrame,
    index_rows: tuple[tuple[str, str, str, str], ...],
) -> None:
    """Apply navigation links and readable formatting to the Index sheet."""
    _format_worksheet(worksheet, dataframe)
    for row_index, sheet_name in enumerate(
        (row[0] for row in index_rows),
        start=2,
    ):
        cell = worksheet.cell(row=row_index, column=1)
        cell.hyperlink = f"#'{sheet_name}'!A1"
        cell.style = "Hyperlink"

    for row in worksheet.iter_rows(min_row=2):
        first_cell = row[0]
        if first_cell.value in {
            "Interaction scoring columns",
            "Interaction scoring notes",
            "Negative evidence columns",
        }:
            first_cell.font = Font(bold=True)


def _add_back_to_index_link(worksheet: Worksheet) -> None:
    """Add a consistent internal link back to the Index sheet."""
    worksheet["A1"] = "Back to Index"
    worksheet["A1"].hyperlink = "#'Index'!A1"
    worksheet["A1"].style = "Hyperlink"


def _column_width(column_name: str, values: pd.Series) -> int:
    """Return a practical Excel column width based on header and cell text."""
    max_length = len(column_name)

    for value in values:
        if pd.isna(value):
            continue

        max_length = max(max_length, len(str(value)))

    return min(max(max_length + 2, 10), 60)


def _record_to_row(record: ProteinRecord) -> dict[str, Any]:
    """Convert one ProteinRecord into one Excel row dictionary."""
    best_positive = _best_hit(record.positive_hits)
    best_negative = _best_hit(record.negative_hits)

    return {
        "protein_id": record.protein_id,
        "description": record.description,
        "old_locus_tag": record.old_locus_tag or "",
        "total_score": record.score.total_score if record.score else 0,
        "score_components": _score_components(record),
        "score_reasons": "; ".join(record.score.reasons) if record.score else "",
        "domain_sources": _unique_join(domain.source for domain in record.domains),
        "domain_names": _join(domain.name for domain in record.domains),
        "domain_accessions": _join(domain.accession for domain in record.domains),
        "domain_descriptions": _join(
            domain.description for domain in record.domains if domain.description
        ),
        "domain_count": len(record.domains),
        "unique_domain_count": _unique_domain_count(record),
        "unique_domain_accessions": _unique_join(
            domain.accession for domain in record.domains
        ),
        "unique_domain_names": _unique_domain_names(record),
        "sequence_length": record.length,
        "positive_hit_count": len(record.positive_hits),
        "negative_hit_count": len(record.negative_hits),
        "blast_status": _blast_status(record),
        "best_positive_hit": best_positive.subject_id if best_positive else None,
        "best_positive_bitscore": best_positive.bitscore if best_positive else None,
        "best_positive_evalue": best_positive.evalue if best_positive else None,
        "best_negative_hit": best_negative.subject_id if best_negative else None,
        "best_negative_bitscore": best_negative.bitscore if best_negative else None,
        "best_negative_evalue": best_negative.evalue if best_negative else None,
        "negative_best_identity": record.negative_best_identity,
        "negative_best_query_coverage": record.negative_best_query_coverage,
        "negative_best_evalue": record.negative_best_evalue,
        "negative_best_source": record.negative_best_source,
        "negative_hit_strength": record.negative_hit_strength or "",
        "negative_strong_hit_count": record.negative_strong_hit_count,
        "negative_medium_hit_count": record.negative_medium_hit_count,
        "negative_weak_hit_count": record.negative_weak_hit_count,
        "negative_exclusion_reason": record.negative_exclusion_reason,
        "motifs": "; ".join(record.motifs),
        "uniprot_accession": record.uniprot_accession,
        "alphafold_url": record.alphafold_url,
        "notes": "; ".join(record.notes),
        "positive_source_count": record.positive_source_count,
        "positive_sources_hit": "; ".join(record.positive_sources_hit),
        "positive_sources_missing": "; ".join(record.positive_sources_missing),
    }


def _positive_source_summary_row(record: ProteinRecord) -> dict[str, Any]:
    """Convert one ProteinRecord into a positive source summary row."""
    row = _record_to_row(record)
    return {
        "protein_id": row["protein_id"],
        "description": row["description"],
        "old_locus_tag": row["old_locus_tag"],
        "negative_hit": bool(record.negative_hits),
        "positive_source_count": row["positive_source_count"],
        "positive_sources_hit": row["positive_sources_hit"],
        "positive_sources_missing": row["positive_sources_missing"],
        "positive_hit_count": row["positive_hit_count"],
        "negative_hit_count": row["negative_hit_count"],
        "blast_status": row["blast_status"],
        "best_positive_hit": row["best_positive_hit"],
        "best_positive_bitscore": row["best_positive_bitscore"],
        "best_positive_evalue": row["best_positive_evalue"],
        "best_negative_hit": row["best_negative_hit"],
        "best_negative_bitscore": row["best_negative_bitscore"],
        "best_negative_evalue": row["best_negative_evalue"],
    }


def _score_components(record: ProteinRecord) -> str:
    """Return score components as readable name=value pairs."""
    if record.score is None:
        return ""

    return "; ".join(
        f"{name}={value}" for name, value in record.score.components.items()
    )


def _unique_domain_count(record: ProteinRecord) -> int:
    """Return the number of unique domain accessions in first-seen order."""
    return len(
        {
            str(domain.accession)
            for domain in record.domains
            if str(domain.accession)
        }
    )


def _unique_domain_names(record: ProteinRecord) -> str:
    """Return unique readable domain names while skipping internal numeric ids."""
    return _unique_join(
        domain.name
        for domain in record.domains
        if not _looks_like_internal_domain_name(domain.name)
    )


def _looks_like_internal_domain_name(value: object) -> bool:
    """Return True for numeric/internal-looking domain names."""
    return bool(re.fullmatch(r"\d{6,}", str(value).strip()))


def _best_hit(hits: list[BlastHit]) -> BlastHit | None:
    """Return the best BLAST hit by highest bitscore, then lowest e-value."""
    if not hits:
        return None

    return max(hits, key=lambda hit: (hit.bitscore, -hit.evalue))


def _join(values: Iterable[object]) -> str:
    """Join non-empty string values with a semicolon separator."""
    return "; ".join(str(value) for value in values if str(value))


def _unique_join(values: Iterable[object]) -> str:
    """Join unique non-empty string values while preserving input order."""
    seen: set[str] = set()
    unique_values: list[str] = []

    for value in values:
        text = str(value)
        if not text or text in seen:
            continue

        seen.add(text)
        unique_values.append(text)

    return "; ".join(unique_values)


def _blast_status(record: ProteinRecord) -> str:
    """Return a compact BLAST hit status for Excel output."""
    has_positive = bool(record.positive_hits)
    has_negative = bool(record.negative_hits)

    if has_positive and has_negative:
        return "positive_and_negative"

    if has_positive:
        return "positive_only"

    if has_negative:
        return "negative_only"

    return "no_hits"


__all__: tuple[str, ...] = (
    "EXCEL_COLUMNS",
    "INDEX_ROWS",
    "INTERACTION_SCORE_EXPLANATIONS",
    "INTERACTION_SCORE_NOTES",
    "POSITIVE_SOURCE_SUMMARY_COLUMNS",
    "positive_source_summary_dataframe",
    "records_to_dataframe",
    "write_classification_workbook",
    "write_records_to_excel",
)
