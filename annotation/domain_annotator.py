"""Apply domain annotations to ProteinRecord objects."""

from __future__ import annotations

from annotation.cdd import search_cdd_by_sequence
from annotation.pfam import search_pfam_by_sequence
from core.cache import JsonCache
from core.models import AnnotationResult, DomainHit, ProteinRecord


DEFAULT_PFAM_EVALUE_THRESHOLD = 1e-5


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


def annotate_pfam_domains(
    record: ProteinRecord,
    cache: JsonCache | None = None,
    timeout: int = 60,
    evalue_threshold: float = DEFAULT_PFAM_EVALUE_THRESHOLD,
) -> ProteinRecord:
    """Annotate one protein record with Pfam domain hits."""
    try:
        domains = search_pfam_by_sequence(
            protein_id=record.protein_id,
            sequence=record.sequence,
            cache=cache,
            timeout=timeout,
        )
    except Exception as exc:
        message = f"Pfam annotation failed for {record.protein_id}: {exc}"
        record.notes.append(message)
        record.annotations["pfam"] = AnnotationResult(
            protein_id=record.protein_id,
            source="pfam",
            success=False,
            error=message,
        )
        return record

    domains = filter_pfam_domains(domains, evalue_threshold=evalue_threshold)
    record.domains.extend(domains)
    record.annotations["pfam"] = AnnotationResult(
        protein_id=record.protein_id,
        source="pfam",
        success=True,
        domains=list(domains),
        metadata={"domain_count": len(domains)},
    )

    return record


def annotate_records_pfam(
    records: dict[str, ProteinRecord],
    cache: JsonCache | None = None,
    timeout: int = 60,
    evalue_threshold: float = DEFAULT_PFAM_EVALUE_THRESHOLD,
) -> dict[str, ProteinRecord]:
    """Annotate all records with Pfam domains and return the same dictionary."""
    for record in records.values():
        annotate_pfam_domains(
            record,
            cache=cache,
            timeout=timeout,
            evalue_threshold=evalue_threshold,
        )

    return records


def filter_pfam_domains(
    domains: list[DomainHit],
    evalue_threshold: float = DEFAULT_PFAM_EVALUE_THRESHOLD,
) -> list[DomainHit]:
    """Keep Pfam hits at or below the e-value threshold.

    Pfam hits with missing e-values are kept for now because some HMMER response
    shapes may omit e-values even when the accession/name is useful.
    """
    filtered: list[DomainHit] = []

    for domain in domains:
        if domain.evalue is None or domain.evalue <= evalue_threshold:
            filtered.append(domain)

    return filtered


__all__: tuple[str, ...] = (
    "annotate_cdd_domains",
    "annotate_pfam_domains",
    "annotate_records_cdd",
    "annotate_records_pfam",
    "filter_pfam_domains",
)
