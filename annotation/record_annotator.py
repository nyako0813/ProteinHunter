"""Apply UniProt and AlphaFold annotations to ProteinRecord objects."""

from __future__ import annotations

from annotation.alphafold import get_alphafold_url_if_exists
from annotation.uniprot import extract_uniprot_accession, search_uniprot_by_protein_id
from core.cache import JsonCache
from core.models import AnnotationResult, ProteinRecord


def annotate_uniprot_and_alphafold(
    record: ProteinRecord,
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> ProteinRecord:
    """Annotate one protein record with UniProt and AlphaFold metadata."""
    try:
        metadata = search_uniprot_by_protein_id(
            record.protein_id,
            cache=cache,
            timeout=timeout,
        )
        accession = extract_uniprot_accession(metadata)
        record.uniprot_accession = accession
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

    if accession is None:
        message = "No UniProt accession was found, so AlphaFold was skipped."
        record.annotations["alphafold"] = AnnotationResult(
            protein_id=record.protein_id,
            source="alphafold",
            success=False,
            error=message,
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


def annotate_records_uniprot_and_alphafold(
    records: dict[str, ProteinRecord],
    cache: JsonCache | None = None,
    timeout: int = 30,
) -> dict[str, ProteinRecord]:
    """Annotate all records and return the same dictionary."""
    for record in records.values():
        annotate_uniprot_and_alphafold(record, cache=cache, timeout=timeout)

    return records


__all__: tuple[str, ...] = (
    "annotate_records_uniprot_and_alphafold",
    "annotate_uniprot_and_alphafold",
)
