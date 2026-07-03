"""High-level BLAST candidate assembly pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

from analysis.candidates import (
    build_candidate_records,
    filter_positive_without_negative,
)
from analysis.ortholog_filter import (
    is_excluded_by_negative_mode,
    populate_negative_hit_evidence,
)
from blast.runner import run_blast_pipeline
from Bio import SeqIO
from core.fasta import read_fasta_as_components
from core.models import ProteinRecord


SOURCE_LABEL_PATTERN = re.compile(r"\[source=([^\]]+)\]")


@dataclass(frozen=True)
class BlastClassificationResult:
    """BLAST classification groups for all target proteins."""

    all_records: dict[str, ProteinRecord]
    positive_only_records: dict[str, ProteinRecord]
    negative_unmatched_records: dict[str, ProteinRecord]
    no_hit_records: dict[str, ProteinRecord]
    negative_hit_records: dict[str, ProteinRecord]
    candidates_relaxed_records: dict[str, ProteinRecord]
    negative_strong_hit_records: dict[str, ProteinRecord]
    negative_medium_hit_records: dict[str, ProteinRecord]
    negative_weak_hit_records: dict[str, ProteinRecord]
    positive_all_sources_records: dict[str, ProteinRecord]
    positive_source_labels: tuple[str, ...]


def run_blast_candidate_pipeline(
    target_fasta: str | Path,
    positive_fasta: str | Path,
    negative_fasta: str | Path,
    work_dir: str | Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 10,
    threads: int = 1,
    positive_source_labels: tuple[str, ...] | None = None,
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
        positive_source_labels=positive_source_labels,
    ).positive_only_records


def run_blast_classification_pipeline(
    target_fasta: str | Path,
    positive_fasta: str | Path,
    negative_fasta: str | Path,
    work_dir: str | Path,
    evalue: float = 1e-5,
    max_target_seqs: int = 10,
    threads: int = 1,
    positive_source_labels: tuple[str, ...] | None = None,
    ortholog_filter: object | None = None,
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
    _populate_hit_query_lengths(positive_hits, sequences)
    _populate_hit_query_lengths(negative_hits, sequences)

    records = build_candidate_records(
        protein_ids=protein_ids,
        descriptions=descriptions,
        sequences=sequences,
        positive_hits=positive_hits,
        negative_hits=negative_hits,
    )
    source_labels = _resolve_positive_source_labels(
        positive_source_labels=positive_source_labels,
        positive_fasta=positive_fasta,
    )
    subject_sources = _positive_subject_sources(
        positive_fasta=positive_fasta,
        fallback_source=source_labels[0] if len(source_labels) == 1 else None,
    )
    negative_subject_sources = _subject_sources(
        fasta_path=negative_fasta,
        fallback_source="negative_fasta",
    )
    _populate_hit_sources(negative_hits, negative_subject_sources)
    _populate_positive_source_summary(
        records=records,
        positive_source_labels=source_labels,
        subject_sources=subject_sources,
    )
    if ortholog_filter is not None:
        populate_negative_hit_evidence(records, ortholog_filter)

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
    if ortholog_filter is None:
        candidates_relaxed_records = dict(positive_only_records)
    else:
        candidates_relaxed_records = {
            protein_id: record
            for protein_id, record in records.items()
            if record.positive_hits
            and not is_excluded_by_negative_mode(
                record,
                ortholog_filter.negative_exclusion_mode,
            )
        }
    negative_strong_hit_records = {
        protein_id: record
        for protein_id, record in records.items()
        if record.negative_strong_hit_count > 0
    }
    negative_medium_hit_records = {
        protein_id: record
        for protein_id, record in records.items()
        if record.negative_medium_hit_count > 0
        and record.negative_strong_hit_count == 0
    }
    negative_weak_hit_records = {
        protein_id: record
        for protein_id, record in records.items()
        if record.negative_weak_hit_count > 0
        and record.negative_strong_hit_count == 0
        and record.negative_medium_hit_count == 0
    }
    positive_all_sources_records = {
        protein_id: record
        for protein_id, record in records.items()
        if not record.negative_hits
        and set(record.positive_sources_hit) == set(source_labels)
    }

    return BlastClassificationResult(
        all_records=records,
        positive_only_records=positive_only_records,
        negative_unmatched_records=negative_unmatched_records,
        no_hit_records=no_hit_records,
        negative_hit_records=negative_hit_records,
        candidates_relaxed_records=candidates_relaxed_records,
        negative_strong_hit_records=negative_strong_hit_records,
        negative_medium_hit_records=negative_medium_hit_records,
        negative_weak_hit_records=negative_weak_hit_records,
        positive_all_sources_records=positive_all_sources_records,
        positive_source_labels=source_labels,
    )


def _resolve_positive_source_labels(
    positive_source_labels: tuple[str, ...] | None,
    positive_fasta: str | Path,
) -> tuple[str, ...]:
    """Return configured or discovered positive source labels."""
    if positive_source_labels:
        return tuple(positive_source_labels)

    discovered = _discover_positive_source_labels(positive_fasta)
    if discovered:
        return discovered

    return ("positive_fasta",)


def _discover_positive_source_labels(positive_fasta: str | Path) -> tuple[str, ...]:
    """Discover source labels from combined FASTA descriptions."""
    fasta_path = Path(positive_fasta)
    if not fasta_path.exists():
        return ()

    labels: list[str] = []
    seen: set[str] = set()
    for record in SeqIO.parse(fasta_path, "fasta"):
        label = _source_label_from_description(record.description)
        if label is None or label in seen:
            continue

        seen.add(label)
        labels.append(label)

    return tuple(labels)


def _positive_subject_sources(
    positive_fasta: str | Path,
    fallback_source: str | None,
) -> dict[str, set[str]]:
    """Map positive FASTA subject IDs to source labels."""
    return _subject_sources(positive_fasta, fallback_source)


def _subject_sources(
    fasta_path: str | Path,
    fallback_source: str | None,
) -> dict[str, set[str]]:
    """Map FASTA subject IDs to source labels."""
    fasta = Path(fasta_path)
    if not fasta.exists():
        return {}

    mapping: dict[str, set[str]] = {}
    for record in SeqIO.parse(fasta, "fasta"):
        source_label = _source_label_from_description(record.description)
        if source_label is None:
            source_label = fallback_source
        if source_label is None:
            continue

        for key in _protein_id_lookup_keys(record.id):
            mapping.setdefault(key, set()).add(source_label)

    return mapping


def _populate_hit_sources(
    hits: list,
    subject_sources: dict[str, set[str]],
) -> None:
    """Replace generic hit source with source labels when they are known."""
    for hit in hits:
        sources = subject_sources.get(hit.subject_id)
        if sources is None:
            sources = subject_sources.get(_without_version(hit.subject_id))
        if sources:
            hit.source = sorted(sources)[0]


def _populate_hit_query_lengths(
    hits: list,
    sequences: dict[str, str],
) -> None:
    """Fill missing query lengths from target FASTA sequences."""
    for hit in hits:
        if hit.query_length is None:
            hit.query_length = len(sequences.get(hit.query_id, ""))


def _populate_positive_source_summary(
    records: dict[str, ProteinRecord],
    positive_source_labels: tuple[str, ...],
    subject_sources: dict[str, set[str]],
) -> None:
    """Populate positive source hit and missing summaries on each record."""
    all_sources = list(positive_source_labels)
    all_source_set = set(all_sources)
    fallback_source = all_sources[0] if len(all_sources) == 1 else None

    for record in records.values():
        hit_sources: set[str] = set()
        for hit in record.positive_hits:
            sources = subject_sources.get(hit.subject_id)
            if sources is None:
                sources = subject_sources.get(_without_version(hit.subject_id))
            if sources is None and fallback_source is not None:
                sources = {fallback_source}
            if sources:
                hit_sources.update(sources)

        ordered_hits = [label for label in all_sources if label in hit_sources]
        missing = [label for label in all_sources if label not in hit_sources]
        record.positive_source_count = len(hit_sources & all_source_set)
        record.positive_sources_hit = ordered_hits
        record.positive_sources_missing = missing


def _source_label_from_description(description: str) -> str | None:
    """Extract a [source=...] label from a FASTA description."""
    match = SOURCE_LABEL_PATTERN.search(description)
    if match:
        return match.group(1)

    return None


def _protein_id_lookup_keys(protein_id: str) -> list[str]:
    """Return versioned and unversioned lookup keys for a protein ID."""
    keys = [protein_id]
    unversioned = _without_version(protein_id)
    if unversioned != protein_id:
        keys.append(unversioned)

    return keys


def _without_version(protein_id: str) -> str:
    """Return a protein ID without a trailing numeric version."""
    return re.sub(r"\.\d+$", "", protein_id)


__all__: tuple[str, ...] = (
    "BlastClassificationResult",
    "SOURCE_LABEL_PATTERN",
    "run_blast_candidate_pipeline",
    "run_blast_classification_pipeline",
)
