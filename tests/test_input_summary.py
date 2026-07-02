"""Tests for input FASTA summary helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from analysis.input_summary import (
    count_fasta_records,
    format_input_summary,
    summarize_input_fastas,
)
from core.exceptions import FileValidationError


def write_fasta(path: Path, records: int) -> Path:
    """Write a small FASTA file with the requested number of records."""
    lines: list[str] = []

    for index in range(1, records + 1):
        lines.append(f">protein_{index}")
        lines.append("MSTNPKPQR")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_count_fasta_records_with_valid_fasta(tmp_path: Path) -> None:
    """A valid FASTA should return the number of records."""
    fasta = write_fasta(tmp_path / "proteins.faa", records=2)

    assert count_fasta_records(fasta) == 2


def test_count_fasta_records_missing_file_raises(tmp_path: Path) -> None:
    """A missing FASTA should raise FileValidationError."""
    with pytest.raises(FileValidationError, match="not found"):
        count_fasta_records(tmp_path / "missing.faa")


def test_count_fasta_records_empty_file_raises(tmp_path: Path) -> None:
    """An empty FASTA should raise FileValidationError."""
    fasta = tmp_path / "empty.faa"
    fasta.write_text("", encoding="utf-8")

    with pytest.raises(FileValidationError, match="empty"):
        count_fasta_records(fasta)


def test_count_fasta_records_no_records_raises(tmp_path: Path) -> None:
    """A non-empty file with no FASTA records should raise FileValidationError."""
    fasta = tmp_path / "not_fasta.faa"
    fasta.write_text("this is not a fasta record\n", encoding="utf-8")

    with pytest.raises(FileValidationError, match="no readable records"):
        count_fasta_records(fasta)


def test_summarize_input_fastas_returns_correct_counts(tmp_path: Path) -> None:
    """Input summary should count target, positive, and negative FASTA files."""
    target = write_fasta(tmp_path / "target.faa", records=3)
    positive = write_fasta(tmp_path / "positive.faa", records=1)
    negative = write_fasta(tmp_path / "negative.faa", records=2)

    assert summarize_input_fastas(target, positive, negative) == {
        "target": 3,
        "positive": 1,
        "negative": 2,
    }


def test_format_input_summary_returns_readable_lines() -> None:
    """Formatted input summary should be easy to read."""
    lines = format_input_summary(
        {
            "target": 1854,
            "positive": 3,
            "negative": 120,
        }
    )

    assert lines == [
        "Target proteins: 1854",
        "Positive references: 3",
        "Negative references: 120",
    ]


def test_format_input_summary_includes_directory_source_counts() -> None:
    """Directory mode summaries should include source folder counts."""
    lines = format_input_summary(
        {
            "target": 1854,
            "positive": 3,
            "negative": 120,
            "target_sources": 2,
            "positive_sources": 4,
            "negative_sources": 1,
        }
    )

    assert lines == [
        "Target proteins: 1854",
        "Positive references: 3",
        "Negative references: 120",
        "Target source folders: 2",
        "Positive source folders: 4",
        "Negative source folders: 1",
    ]
