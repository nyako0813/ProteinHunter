"""Tests for BLAST runner utilities."""

from __future__ import annotations

import subprocess
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from blast.runner import (
    make_blast_db,
    parse_blast_tabular,
    run_blast_pipeline,
    run_blastp,
    validate_fasta,
)
from core.exceptions import BlastDatabaseError, BlastError, BlastExecutionError, BlastParseError
from core.models import BlastHit


def write_fasta(path: Path) -> Path:
    """Write a tiny protein FASTA file for tests."""
    path.write_text(">protein_1\nMSTNPKPQR\n", encoding="utf-8")
    return path


def test_validate_fasta_success(tmp_path: Path) -> None:
    """A present non-empty FASTA should resolve successfully."""
    fasta = write_fasta(tmp_path / "query.faa")

    assert validate_fasta(fasta) == fasta.resolve()


def test_validate_fasta_missing(tmp_path: Path) -> None:
    """Missing FASTA files should raise BlastError."""
    with pytest.raises(BlastError, match="not found"):
        validate_fasta(tmp_path / "missing.faa")


def test_validate_fasta_empty(tmp_path: Path) -> None:
    """Empty FASTA files should raise BlastError."""
    fasta = tmp_path / "empty.faa"
    fasta.write_text("", encoding="utf-8")

    with pytest.raises(BlastError, match="empty"):
        validate_fasta(fasta)


def test_parse_blast_tabular_valid_output(tmp_path: Path) -> None:
    """Valid BLAST tabular output should become BlastHit objects."""
    output = tmp_path / "blast.tsv"
    output.write_text(
        "query_1\tsubject_1\t95.5\t120\t1e-20\t88.0\n"
        "query_2\tsubject_2\t42.0\t75\t0.001\t31.5\n",
        encoding="utf-8",
    )

    hits = parse_blast_tabular(output, source="positive")

    assert hits == [
        BlastHit("query_1", "subject_1", 95.5, 120, 1e-20, 88.0, "positive"),
        BlastHit("query_2", "subject_2", 42.0, 75, 0.001, 31.5, "positive"),
    ]


def test_parse_blast_tabular_empty_output(tmp_path: Path) -> None:
    """Empty BLAST output should parse as no hits."""
    output = tmp_path / "empty.tsv"
    output.write_text("", encoding="utf-8")

    assert parse_blast_tabular(output) == []


def test_parse_blast_tabular_malformed_line(tmp_path: Path) -> None:
    """Malformed BLAST output should include the line number."""
    output = tmp_path / "bad.tsv"
    output.write_text("query_1\tsubject_1\t95.5\n", encoding="utf-8")

    with pytest.raises(BlastParseError, match="line 1"):
        parse_blast_tabular(output)


def test_parse_blast_tabular_bad_number(tmp_path: Path) -> None:
    """Invalid numeric columns should raise BlastParseError."""
    output = tmp_path / "bad_number.tsv"
    output.write_text("query_1\tsubject_1\tbad\t120\t1e-20\t88.0\n", encoding="utf-8")

    with pytest.raises(BlastParseError, match="line 1"):
        parse_blast_tabular(output)


def test_make_blast_db_runs_makeblastdb(tmp_path: Path) -> None:
    """make_blast_db should call makeblastdb and return the database prefix."""
    fasta = write_fasta(tmp_path / "subject.faa")

    with patch("blast.runner.subprocess.run") as run_mock:
        db_prefix = make_blast_db(fasta, tmp_path / "db", "subject_db")

    assert db_prefix == (tmp_path / "db" / "subject_db").resolve()
    run_mock.assert_called_once()
    command = run_mock.call_args.args[0]
    assert command[:2] == ["makeblastdb", "-in"]
    assert "-dbtype" in command
    assert "prot" in command
    assert run_mock.call_args.kwargs == {
        "check": True,
        "capture_output": True,
        "text": True,
    }


def test_make_blast_db_failure_raises_database_error(tmp_path: Path) -> None:
    """makeblastdb command failures should raise BlastDatabaseError."""
    fasta = write_fasta(tmp_path / "subject.faa")
    error = subprocess.CalledProcessError(
        returncode=1,
        cmd=["makeblastdb"],
        stderr="database failed",
    )

    with patch("blast.runner.subprocess.run", side_effect=error):
        with pytest.raises(BlastDatabaseError, match="database failed"):
            make_blast_db(fasta, tmp_path / "db", "subject_db")


def test_run_blastp_runs_blastp(tmp_path: Path) -> None:
    """run_blastp should call blastp and return the output path."""
    query = write_fasta(tmp_path / "query.faa")
    output = tmp_path / "out" / "hits.tsv"

    with patch("blast.runner.subprocess.run") as run_mock:
        result = run_blastp(query, tmp_path / "db" / "subject_db", output)

    assert result == output.resolve()
    run_mock.assert_called_once()
    command = run_mock.call_args.args[0]
    assert command[0] == "blastp"
    assert "-outfmt" in command
    assert "6 qseqid sseqid pident length evalue bitscore" in command
    assert output.parent.exists()


def test_run_blastp_failure_raises_execution_error(tmp_path: Path) -> None:
    """blastp command failures should raise BlastExecutionError."""
    query = write_fasta(tmp_path / "query.faa")
    error = subprocess.CalledProcessError(
        returncode=1,
        cmd=["blastp"],
        stderr="search failed",
    )

    with patch("blast.runner.subprocess.run", side_effect=error):
        with pytest.raises(BlastExecutionError, match="search failed"):
            run_blastp(query, tmp_path / "db" / "subject_db", tmp_path / "hits.tsv")


def test_run_blast_pipeline_uses_helpers(tmp_path: Path) -> None:
    """The pipeline should create a database, run BLAST, and parse output."""
    query = tmp_path / "query.faa"
    subject = tmp_path / "subject.faa"
    output = tmp_path / "work" / "db_name_blast.tsv"
    hits = [BlastHit("q", "s", 90.0, 10, 1e-5, 50.0, "positive")]

    with (
        patch("blast.runner.make_blast_db", return_value=tmp_path / "db" / "db_name") as db_mock,
        patch("blast.runner.run_blastp", return_value=output) as blastp_mock,
        patch("blast.runner.parse_blast_tabular", return_value=hits) as parse_mock,
    ):
        result = run_blast_pipeline(
            query_fasta=query,
            subject_fasta=subject,
            work_dir=tmp_path / "work",
            db_name="db_name",
            source="positive",
        )

    assert result == hits
    db_mock.assert_called_once_with(subject, (tmp_path / "work").resolve() / "db", "db_name", protein=True)
    blastp_mock.assert_called_once()
    parse_mock.assert_called_once_with(output, source="positive")
