"""Input FASTA summary helpers for ProteinHunter."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from core.exceptions import FileValidationError
from core.fasta import validate_fasta_file


def count_fasta_records(path: str | Path) -> int:
    """Validate a FASTA file and return its record count."""
    fasta_path = validate_fasta_file(path)

    try:
        count = sum(1 for _record in SeqIO.parse(fasta_path, "fasta"))
    except ValueError as exc:
        raise FileValidationError(
            f"The FASTA file contains no readable records: {fasta_path}"
        ) from exc

    if count == 0:
        raise FileValidationError(
            f"The FASTA file contains no readable records: {fasta_path}"
        )

    return count


def summarize_input_fastas(
    target_fasta: str | Path,
    positive_fasta: str | Path,
    negative_fasta: str | Path,
    source_counts: dict[str, int] | None = None,
) -> dict[str, int]:
    """Return record counts for target, positive, and negative FASTA files."""
    summary = {
        "target": count_fasta_records(target_fasta),
        "positive": count_fasta_records(positive_fasta),
        "negative": count_fasta_records(negative_fasta),
    }
    if source_counts:
        summary.update(source_counts)

    return summary


def format_input_summary(summary: dict[str, int]) -> list[str]:
    """Return beginner-friendly input summary lines."""
    lines = [
        f"Target proteins: {summary.get('target', 0)}",
        f"Positive references: {summary.get('positive', 0)}",
        f"Negative references: {summary.get('negative', 0)}",
    ]
    if "target_sources" in summary:
        lines.append(f"Target source folders: {summary.get('target_sources', 0)}")
    if "positive_sources" in summary:
        lines.append(f"Positive source folders: {summary.get('positive_sources', 0)}")
    if "negative_sources" in summary:
        lines.append(f"Negative source folders: {summary.get('negative_sources', 0)}")

    return lines


__all__: tuple[str, ...] = (
    "count_fasta_records",
    "format_input_summary",
    "summarize_input_fastas",
)
