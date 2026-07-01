"""Input FASTA summary helpers for ProteinHunter."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from core.exceptions import FileValidationError
from core.fasta import validate_fasta_file


def count_fasta_records(path: str | Path) -> int:
    """Validate a FASTA file and return its record count."""
    fasta_path = validate_fasta_file(path)
    count = sum(1 for _record in SeqIO.parse(fasta_path, "fasta"))

    if count == 0:
        raise FileValidationError(
            f"The FASTA file contains no readable records: {fasta_path}"
        )

    return count


def summarize_input_fastas(
    target_fasta: str | Path,
    positive_fasta: str | Path,
    negative_fasta: str | Path,
) -> dict[str, int]:
    """Return record counts for target, positive, and negative FASTA files."""
    return {
        "target": count_fasta_records(target_fasta),
        "positive": count_fasta_records(positive_fasta),
        "negative": count_fasta_records(negative_fasta),
    }


def format_input_summary(summary: dict[str, int]) -> list[str]:
    """Return beginner-friendly input summary lines."""
    return [
        f"Target proteins: {summary.get('target', 0)}",
        f"Positive references: {summary.get('positive', 0)}",
        f"Negative references: {summary.get('negative', 0)}",
    ]


__all__: tuple[str, ...] = (
    "count_fasta_records",
    "format_input_summary",
    "summarize_input_fastas",
)
