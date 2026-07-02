"""Apply UniProt and AlphaFold annotations to ProteinRecord objects."""

from __future__ import annotations

from annotation.alphafold import get_alphafold_url_if_exists
from annotation.uniprot import (
    extract_uniprot_accession,
    extract_uniprot_old_locus_tag,
    search_uniprot_by_protein_id,
)
from core.cache import JsonCache
from core.models import AnnotationResult, ProteinRecord


def annotate_uniprot_and_alphafold(
    record: ProteinRecord,
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> ProteinRecord:
    """Annotate one protein record with UniProt and AlphaFold metadata."""
    annotate_uniprot(record, cache=cache, timeout=timeout)
    annotate_alphafold(record, cache=cache, timeout=timeout)
    return record


def annotate_uniprot(
    record: ProteinRecord,
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> ProteinRecord:
    """Annotate one protein record with UniProt metadata."""
    try:
        metadata = search_uniprot_by_protein_id(
            record.protein_id,
            cache=cache,
            timeout=timeout,
        )
        accession = extract_uniprot_accession(metadata)
        old_locus_tag = extract_uniprot_old_locus_tag(metadata)
        record.uniprot_accession = accession
        if old_locus_tag is not None:
            record.old_locus_tag = old_locus_tag
        record.annotations["uniprot"] = AnnotationResult(
            protein_id=record.protein_id,
            source="uniprot",
            success=True,
            metadata=metadata,
        )
    except Exception as exc:
        message = f"UniProt annotation failed for {record.protein_id}: {exc}"
        record.notes.append(message)
        record.annotations["uniprot"] = AnnotationResult(
            protein_id=record.protein_id,
            source="uniprot",
            success=False,
            error=message,
        )
        return record

    return record


def annotate_records_uniprot(
    records: dict[str, ProteinRecord],
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> dict[str, ProteinRecord]:
    """Annotate all records with UniProt metadata and return the same dictionary."""
    for record in records.values():
        annotate_uniprot(record, cache=cache, timeout=timeout)

    return records


def annotate_alphafold(
    record: ProteinRecord,
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> ProteinRecord:
    """Annotate one protein record with an AlphaFold URL when possible."""
    accession = record.uniprot_accession

    if accession is None or not accession.strip():
        message = "AlphaFold annotation skipped because UniProt accession is missing."
        record.notes.append(message)
        record.annotations["alphafold"] = AnnotationResult(
            protein_id=record.protein_id,
            source="alphafold",
            success=True,
            metadata={"accession": None, "url": None, "exists": False},
        )
        return record

    try:
        alphafold_url = get_alphafold_url_if_exists(
            accession,
            cache=cache,
            timeout=timeout,
        )
        record.alphafold_url = alphafold_url
        record.annotations["alphafold"] = AnnotationResult(
            protein_id=record.protein_id,
            source="alphafold",
            success=True,
            metadata={
                "accession": accession,
                "url": alphafold_url,
                "exists": alphafold_url is not None,
            },
        )
    except Exception as exc:
        message = f"AlphaFold annotation failed for {record.protein_id}: {exc}"
        record.notes.append(message)
        record.annotations["alphafold"] = AnnotationResult(
            protein_id=record.protein_id,
            source="alphafold",
            success=False,
            error=message,
            metadata={"accession": accession, "url": None, "exists": False},
        )

    return record


def annotate_records_alphafold(
    records: dict[str, ProteinRecord],
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> dict[str, ProteinRecord]:
    """Annotate all records with AlphaFold URLs and return the same dictionary."""
    for record in records.values():
        annotate_alphafold(record, cache=cache, timeout=timeout)

    return records


def annotate_records_uniprot_and_alphafold(
    records: dict[str, ProteinRecord],
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> dict[str, ProteinRecord]:
    """Annotate all records and return the same dictionary."""
    annotate_records_uniprot(records, cache=cache, timeout=timeout)
    annotate_records_alphafold(records, cache=cache, timeout=timeout)
    return records


__all__: tuple[str, ...] = (
    "annotate_alphafold",
    "annotate_records_alphafold",
    "annotate_records_uniprot",
    "annotate_records_uniprot_and_alphafold",
    "annotate_uniprot",
    "annotate_uniprot_and_alphafold",
)
