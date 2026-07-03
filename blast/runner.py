"""Utility functions for creating and running BLAST searches.

This module wraps the external BLAST+ commands with small, typed helpers and
parses BLAST tabular output into ProteinHunter data models.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from core.exceptions import (
    BlastDatabaseError,
    BlastError,
    BlastExecutionError,
    BlastParseError,
)
from core.models import BlastHit


BLAST_OUTFMT_COLUMNS: tuple[str, ...] = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "evalue",
    "bitscore",
)


def validate_fasta(path: str | Path) -> Path:
    """Return a resolved FASTA path if it exists and contains data."""
    fasta_path = Path(path).expanduser().resolve()

    if not fasta_path.exists():
        raise BlastError(f"The FASTA file was not found: {fasta_path}")

    if not fasta_path.is_file():
        raise BlastError(f"The FASTA path is not a file: {fasta_path}")

    if fasta_path.stat().st_size == 0:
        raise BlastError(f"The FASTA file is empty: {fasta_path}")

    return fasta_path


def make_blast_db(
    fasta_path: str | Path,
    db_dir: str | Path,
    db_name: str,
    protein: bool = True,
) -> Path:
    """Create a BLAST database and return its prefix path."""
    validated_fasta = validate_fasta(fasta_path)
    database_dir = Path(db_dir).expanduser().resolve()
    database_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = database_dir / db_name
    dbtype = "prot" if protein else "nucl"

    command = [
        "makeblastdb",
        "-in",
        str(validated_fasta),
        "-dbtype",
        dbtype,
        "-out",
        str(db_prefix),
    ]

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise BlastDatabaseError(
            "The makeblastdb command was not found. Please install BLAST+."
        ) from exc
    except subprocess.CalledProcessError as exc:
        details = _subprocess_error_text(exc)
        raise BlastDatabaseError(
            f"BLAST database creation failed for '{validated_fasta}'. {details}"
        ) from exc

    return db_prefix


def run_blastp(
    query_fasta: str | Path,
    db_prefix: str | Path,
    output_path: str | Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 10,
    threads: int = 1,
) -> Path:
    """Run blastp and return the tabular output path."""
    validated_query = validate_fasta(query_fasta)
    database_prefix = Path(db_prefix).expanduser().resolve()
    blast_output = Path(output_path).expanduser().resolve()
    blast_output.parent.mkdir(parents=True, exist_ok=True)
    outfmt = "6 " + " ".join(BLAST_OUTFMT_COLUMNS)

    command = [
        "blastp",
        "-query",
        str(validated_query),
        "-db",
        str(database_prefix),
        "-out",
        str(blast_output),
        "-outfmt",
        outfmt,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(threads),
    ]

    try:
        subprocess.run(command, check=True, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise BlastExecutionError(
            "The blastp command was not found. Please install BLAST+."
        ) from exc
    except subprocess.CalledProcessError as exc:
        details = _subprocess_error_text(exc)
        raise BlastExecutionError(
            f"BLAST search failed for query '{validated_query}'. {details}"
        ) from exc

    return blast_output


def parse_blast_tabular(path: str | Path, source: str = "blast") -> list[BlastHit]:
    """Parse BLAST outfmt 6 output into a list of BlastHit objects."""
    blast_output = Path(path).expanduser().resolve()

    if not blast_output.exists():
        raise BlastParseError(f"The BLAST output file was not found: {blast_output}")

    hits: list[BlastHit] = []

    with blast_output.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line:
                continue

            columns = line.split("\t")
            if len(columns) not in {6, len(BLAST_OUTFMT_COLUMNS)}:
                raise BlastParseError(
                    f"Malformed BLAST output on line {line_number}: "
                    f"expected 6 or {len(BLAST_OUTFMT_COLUMNS)} columns, "
                    f"got {len(columns)}."
                )

            try:
                query_length = int(columns[4]) if len(columns) == 7 else None
                evalue_index = 5 if len(columns) == 7 else 4
                bitscore_index = 6 if len(columns) == 7 else 5
                hit = BlastHit(
                    query_id=columns[0],
                    subject_id=columns[1],
                    percent_identity=float(columns[2]),
                    alignment_length=int(columns[3]),
                    evalue=float(columns[evalue_index]),
                    bitscore=float(columns[bitscore_index]),
                    source=source,
                    query_length=query_length,
                )
            except ValueError as exc:
                raise BlastParseError(
                    f"Malformed BLAST output on line {line_number}: "
                    "numeric columns could not be read."
                ) from exc

            hits.append(hit)

    return hits


def run_blast_pipeline(
    query_fasta: str | Path,
    subject_fasta: str | Path,
    work_dir: str | Path,
    db_name: str,
    source: str,
    evalue: float = 1e-5,
    max_target_seqs: int = 10,
    threads: int = 1,
) -> list[BlastHit]:
    """Create a subject database, run blastp, and parse the tabular hits."""
    working_dir = Path(work_dir).expanduser().resolve()
    db_dir = working_dir / "db"
    output_path = working_dir / f"{db_name}_blast.tsv"

    db_prefix = make_blast_db(subject_fasta, db_dir, db_name, protein=True)
    blast_output = run_blastp(
        query_fasta=query_fasta,
        db_prefix=db_prefix,
        output_path=output_path,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        threads=threads,
    )

    return parse_blast_tabular(blast_output, source=source)


def _subprocess_error_text(error: subprocess.CalledProcessError) -> str:
    """Return a readable error message from a failed subprocess call."""
    stderr = (error.stderr or "").strip()
    stdout = (error.stdout or "").strip()

    if stderr:
        return stderr

    if stdout:
        return stdout

    return "No extra error message was provided."


__all__: tuple[str, ...] = (
    "BLAST_OUTFMT_COLUMNS",
    "make_blast_db",
    "parse_blast_tabular",
    "run_blast_pipeline",
    "run_blastp",
    "validate_fasta",
)
