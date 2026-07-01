"""Tests for FASTA parsing helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from core.exceptions import FileValidationError
from core.fasta import (
    read_fasta_as_components,
    read_fasta_descriptions,
    read_fasta_ids,
    read_fasta_sequences,
    validate_fasta_file,
)


def write_two_record_fasta(path: Path) -> Path:
    """Write a small two-record FASTA file."""
    path.write_text(
        ">protein_1 first candidate protein\nMSTNPKPQR\n"
        ">protein_2 second candidate protein\nAAAA\nCCCC\n",
        encoding="utf-8",
    )
    return path


def test_validate_fasta_file_accepts_valid_fasta(tmp_path: Path) -> None:
    """A present non-empty FASTA file should validate."""
    fasta = write_two_record_fasta(tmp_path / "proteins.faa")

    assert validate_fasta_file(fasta) == fasta.resolve()


def test_read_fasta_ids_preserves_order(tmp_path: Path) -> None:
    """FASTA IDs should be returned in the same order as the file."""
    fasta = write_two_record_fasta(tmp_path / "proteins.faa")

    assert read_fasta_ids(fasta) == ["protein_1", "protein_2"]


def test_read_fasta_sequences(tmp_path: Path) -> None:
    """FASTA sequences should be keyed by record ID."""
    fasta = write_two_record_fasta(tmp_path / "proteins.faa")

    assert read_fasta_sequences(fasta) == {
        "protein_1": "MSTNPKPQR",
        "protein_2": "AAAACCCC",
    }


def test_read_fasta_descriptions(tmp_path: Path) -> None:
    """FASTA descriptions should be keyed by record ID."""
    fasta = write_two_record_fasta(tmp_path / "proteins.faa")

    assert read_fasta_descriptions(fasta) == {
        "protein_1": "protein_1 first candidate protein",
        "protein_2": "protein_2 second candidate protein",
    }


def test_read_fasta_as_components(tmp_path: Path) -> None:
    """FASTA components should include IDs, descriptions, and sequences."""
    fasta = write_two_record_fasta(tmp_path / "proteins.faa")

    protein_ids, descriptions, sequences = read_fasta_as_components(fasta)

    assert protein_ids == ["protein_1", "protein_2"]
    assert descriptions == {
        "protein_1": "protein_1 first candidate protein",
        "protein_2": "protein_2 second candidate protein",
    }
    assert sequences == {
        "protein_1": "MSTNPKPQR",
        "protein_2": "AAAACCCC",
    }


def test_missing_file_raises_file_validation_error(tmp_path: Path) -> None:
    """A missing FASTA file should raise FileValidationError."""
    with pytest.raises(FileValidationError, match="not found"):
        validate_fasta_file(tmp_path / "missing.faa")


def test_empty_file_raises_file_validation_error(tmp_path: Path) -> None:
    """An empty FASTA file should raise FileValidationError."""
    fasta = tmp_path / "empty.faa"
    fasta.write_text("", encoding="utf-8")

    with pytest.raises(FileValidationError, match="empty"):
        read_fasta_ids(fasta)


def test_fasta_with_no_records_raises_file_validation_error(tmp_path: Path) -> None:
    """A non-empty file with no FASTA records should raise FileValidationError."""
    fasta = tmp_path / "no_records.faa"
    fasta.write_text("this is not a fasta record\n", encoding="utf-8")

    with pytest.raises(FileValidationError, match="no readable records"):
        read_fasta_sequences(fasta)
