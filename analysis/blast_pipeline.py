"""High-level BLAST candidate assembly pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from analysis.candidates import (
    build_candidate_records,
    filter_positive_without_negative,
)
from blast.runner import run_blast_pipeline
from core.fasta import read_fasta_as_components
from core.models import ProteinRecord


@dataclass(frozen=True)
class BlastClassificationResult:
    """BLAST classification groups for all target proteins."""

    all_records: dict[str, ProteinRecord]
    positive_only_records: dict[str, ProteinRecord]
    negative_unmatched_records: dict[str, ProteinRecord]
    no_hit_records: dict[str, ProteinRecord]
    negative_hit_records: dict[str, ProteinRecord]


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
    return run_blast_classification_pipeline(
        target_fasta=target_fasta,
        positive_fasta=positive_fasta,
        negative_fasta=negative_fasta,
        work_dir=work_dir,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        threads=threads,
    ).positive_only_records


def run_blast_classification_pipeline(
    target_fasta: str | Path,
    positive_fasta: str | Path,
    negative_fasta: str | Path,
    work_dir: str | Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 10,
    threads: int = 1,
) -> BlastClassificationResult:
    """Build BLAST classification groups for all target proteins."""
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

    positive_only_records = filter_positive_without_negative(records)
    negative_unmatched_records = {
        protein_id: record
        for protein_id, record in records.items()
        if not record.negative_hits
    }
    no_hit_records = {
        protein_id: record
        for protein_id, record in records.items()
        if not record.positive_hits and not record.negative_hits
    }
    negative_hit_records = {
        protein_id: record
        for protein_id, record in records.items()
        if record.negative_hits
    }

    return BlastClassificationResult(
        all_records=records,
        positive_only_records=positive_only_records,
        negative_unmatched_records=negative_unmatched_records,
        no_hit_records=no_hit_records,
        negative_hit_records=negative_hit_records,
    )


__all__: tuple[str, ...] = (
    "BlastClassificationResult",
    "run_blast_candidate_pipeline",
    "run_blast_classification_pipeline",
)
