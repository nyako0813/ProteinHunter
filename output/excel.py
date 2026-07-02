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
) -> Path:
    """Write candidate records and BLAST classification sheets to Excel."""
    resolved_output = Path(output_path).expanduser().resolve()
    resolved_output.parent.mkdir(parents=True, exist_ok=True)
    sheets = {
        "Candidates": records_to_dataframe(candidates),
        "Positive_all_sources": records_to_dataframe(positive_all_sources or {}),
        "Positive_source_summary": positive_source_summary_dataframe(
            positive_source_summary or {}
        ),
        "Negative_unmatched": records_to_dataframe(negative_unmatched or {}),
        "No_hit": records_to_dataframe(no_hit or {}),
        "Negative_hit": records_to_dataframe(negative_hit or {}),
    }

    try:
        with pd.ExcelWriter(resolved_output, engine="openpyxl") as writer:
            for sheet_name, dataframe in sheets.items():
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


def positive_source_summary_dataframe(
    records: dict[str, ProteinRecord],
) -> pd.DataFrame:
    """Return the compact positive-source summary DataFrame."""
    rows = [_positive_source_summary_row(record) for record in records.values()]
    return pd.DataFrame(rows, columns=POSITIVE_SOURCE_SUMMARY_COLUMNS)


def _format_worksheet(worksheet: Worksheet, dataframe: pd.DataFrame) -> None:
    """Apply simple readability formatting to an Excel worksheet."""
    worksheet.freeze_panes = "A2"
    worksheet.auto_filter.ref = worksheet.dimensions

    for header_cell in worksheet[1]:
        header_cell.font = Font(bold=True)

    wrap_columns = {
        "description",
        "domain_descriptions",
        "score_reasons",
        "notes",
        "positive_sources_hit",
        "positive_sources_missing",
    }

    for column_index, column_name in enumerate(dataframe.columns, start=1):
        column_letter = get_column_letter(column_index)
        width = _column_width(column_name, dataframe[column_name])

        if column_name in {"description", "domain_descriptions", "score_reasons"}:
            width = max(width, 35)
        elif column_name == "notes":
            width = max(width, 45)

        worksheet.column_dimensions[column_letter].width = width

        if column_name in wrap_columns:
            for cell in worksheet[column_letter]:
                cell.alignment = Alignment(wrap_text=True, vertical="top")


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
    "POSITIVE_SOURCE_SUMMARY_COLUMNS",
    "positive_source_summary_dataframe",
    "records_to_dataframe",
    "write_classification_workbook",
    "write_records_to_excel",
)
