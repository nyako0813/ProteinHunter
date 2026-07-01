"""Apply domain annotations to ProteinRecord objects."""

from __future__ import annotations

from annotation.cdd import search_cdd_by_sequence
from core.cache import JsonCache
from core.models import AnnotationResult, ProteinRecord


def annotate_cdd_domains(
    record: ProteinRecord,
    cache: JsonCache | None = None,
    timeout: int = 60,
) -> ProteinRecord:
    """Annotate one protein record with CDD domain hits."""
    try:
        domains = search_cdd_by_sequence(
            protein_id=record.protein_id,
            sequence=record.sequence,
            cache=cache,
            timeout=timeout,
        )
    except Exception as exc:
        message = f"CDD annotation failed for {record.protein_id}: {exc}"
        record.notes.append(message)
        record.annotations["cdd"] = AnnotationResult(
            protein_id=record.protein_id,
            source="cdd",
            success=False,
            error=message,
        )
        return record

    record.domains.extend(domains)
    record.annotations["cdd"] = AnnotationResult(
        protein_id=record.protein_id,
        source="cdd",
        success=True,
        domains=list(domains),
        metadata={"domain_count": len(domains)},
    )

    return record


def annotate_records_cdd(
    records: dict[str, ProteinRecord],
    cache: JsonCache | None = None,
    timeout: int = 60,
) -> dict[str, ProteinRecord]:
    """Annotate all records with CDD domains and return the same dictionary."""
    for record in records.values():
        annotate_cdd_domains(record, cache=cache, timeout=timeout)

    return records


__all__: tuple[str, ...] = (
    "annotate_cdd_domains",
    "annotate_records_cdd",
)
