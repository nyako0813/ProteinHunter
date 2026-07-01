"""Minimal Excel output helpers for ProteinHunter candidate records."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd

from core.exceptions import ExcelOutputError
from core.models import BlastHit, ProteinRecord


EXCEL_COLUMNS: tuple[str, ...] = (
    "protein_id",
    "description",
    "total_score",
    "score_components",
    "score_reasons",
    "domain_sources",
    "domain_names",
    "domain_accessions",
    "domain_descriptions",
    "domain_count",
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
    except Exception as exc:
        message = (
            f"ProteinHunter could not write the Excel file: {resolved_output}. "
            "Please check that the folder is writable and the file is not open."
        )
        raise ExcelOutputError(message) from exc

    return resolved_output


def _record_to_row(record: ProteinRecord) -> dict[str, Any]:
    """Convert one ProteinRecord into one Excel row dictionary."""
    best_positive = _best_hit(record.positive_hits)
    best_negative = _best_hit(record.negative_hits)

    return {
        "protein_id": record.protein_id,
        "description": record.description,
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
    }


def _score_components(record: ProteinRecord) -> str:
    """Return score components as readable name=value pairs."""
    if record.score is None:
        return ""

    return "; ".join(
        f"{name}={value}" for name, value in record.score.components.items()
    )


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
    "records_to_dataframe",
    "write_records_to_excel",
)
