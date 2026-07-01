"""FASTA file parsing helpers for ProteinHunter."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from core.exceptions import FileValidationError


def validate_fasta_file(path: str | Path) -> Path:
    """Return a resolved FASTA path if it exists and is not empty."""
    fasta_path = Path(path).expanduser().resolve()

    if not fasta_path.exists():
        raise FileValidationError(f"The FASTA file was not found: {fasta_path}")

    if not fasta_path.is_file():
        raise FileValidationError(f"The FASTA path is not a file: {fasta_path}")

    if fasta_path.stat().st_size == 0:
        raise FileValidationError(f"The FASTA file is empty: {fasta_path}")

    return fasta_path


def read_fasta_sequences(path: str | Path) -> dict[str, str]:
    """Read a FASTA file and return sequences keyed by record ID."""
    records = _read_records(path)
    return {record.id: str(record.seq) for record in records}


def read_fasta_descriptions(path: str | Path) -> dict[str, str]:
    """Read a FASTA file and return descriptions keyed by record ID."""
    records = _read_records(path)
    return {record.id: record.description for record in records}


def read_fasta_ids(path: str | Path) -> list[str]:
    """Read a FASTA file and return record IDs in input order."""
    records = _read_records(path)
    return [record.id for record in records]


def read_fasta_as_components(
    path: str | Path,
) -> tuple[list[str], dict[str, str], dict[str, str]]:
    """Read a FASTA file as IDs, descriptions, and sequences."""
    records = _read_records(path)
    protein_ids = [record.id for record in records]
    descriptions = {record.id: record.description for record in records}
    sequences = {record.id: str(record.seq) for record in records}

    return protein_ids, descriptions, sequences


def _read_records(path: str | Path) -> list[SeqIO.SeqRecord]:
    """Validate and read FASTA records, raising a friendly error if empty."""
    fasta_path = validate_fasta_file(path)
    records = list(SeqIO.parse(fasta_path, "fasta"))

    if not records:
        raise FileValidationError(
            f"The FASTA file contains no readable records: {fasta_path}"
        )

    return records


__all__: tuple[str, ...] = (
    "read_fasta_as_components",
    "read_fasta_descriptions",
    "read_fasta_ids",
    "read_fasta_sequences",
    "validate_fasta_file",
)
