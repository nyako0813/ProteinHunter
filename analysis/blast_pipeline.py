"""High-level BLAST candidate assembly pipeline."""

from __future__ import annotations

from pathlib import Path

from analysis.candidates import (
    build_candidate_records,
    filter_positive_without_negative,
)
from blast.runner import run_blast_pipeline
from core.fasta import read_fasta_as_components
from core.models import ProteinRecord


def run_blast_candidate_pipeline(
    target_fasta: str | Path,
    positive_fasta: str | Path,
    negative_fasta: str | Path,
    work_dir: str | Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 10,
    threads: int = 1,
) -> dict[str, ProteinRecord]:
    """Build candidate records from positive and negative BLAST searches."""
    protein_ids, descriptions, sequences = read_fasta_as_components(target_fasta)
    working_dir = Path(work_dir).expanduser().resolve()

    positive_hits = run_blast_pipeline(
        query_fasta=target_fasta,
        subject_fasta=positive_fasta,
        work_dir=working_dir / "positive",
        db_name="positive",
        source="positive",
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        threads=threads,
    )
    negative_hits = run_blast_pipeline(
        query_fasta=target_fasta,
        subject_fasta=negative_fasta,
        work_dir=working_dir / "negative",
        db_name="negative",
        source="negative",
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        threads=threads,
    )

    records = build_candidate_records(
        protein_ids=protein_ids,
        descriptions=descriptions,
        sequences=sequences,
        positive_hits=positive_hits,
        negative_hits=negative_hits,
    )

    return filter_positive_without_negative(records)


__all__: tuple[str, ...] = ("run_blast_candidate_pipeline",)
