"""Apply domain annotations to ProteinRecord objects."""

from __future__ import annotations

from annotation.cdd import (
    CDD_MAX_POLL_SECONDS,
    CDD_MAX_QUERIES_PER_BATCH,
    CDD_POLL_INTERVAL_SECONDS,
    domain_hit_from_dict,
    domain_hit_to_dict,
    fetch_cdd_batch_results,
    parse_cdd_batch_response,
    poll_cdd_batch,
    search_cdd_by_sequence,
    submit_cdd_batch,
)
from annotation.pfam import enrich_pfam_domains_with_metadata, search_pfam_by_sequence
from core.cache import JsonCache
from core.models import AnnotationResult, DomainHit, ProteinRecord


DEFAULT_PFAM_EVALUE_THRESHOLD = 1e-5


def annotate_cdd_domains(
    record: ProteinRecord,
    cache: JsonCache | None = None,
    timeout: int = 60,
) -> ProteinRecord:
    """Annotate one protein record with CDD domain hits.

    Submits a single-sequence Batch CD-Search job. For annotating many
    records at once, prefer :func:`annotate_records_cdd`, which batches
    cache-miss records into shared jobs instead of one job per protein.
    """
    try:
        domains = search_cdd_by_sequence(
            protein_id=record.protein_id,
            sequence=record.sequence,
            cache=cache,
            timeout=timeout,
        )
    except Exception as exc:
        _apply_cdd_failure(record, exc)
        return record

    _apply_cdd_success(record, domains)
    return record


def annotate_records_cdd(
    records: dict[str, ProteinRecord],
    cache: JsonCache | None = None,
    timeout: int = 60,
    poll_interval: float = CDD_POLL_INTERVAL_SECONDS,
    max_wait: float = CDD_MAX_POLL_SECONDS,
) -> dict[str, ProteinRecord]:
    """Annotate all records with CDD domains and return the same dictionary.

    Cache-hit records are filled in without a network call. Cache-miss
    records are submitted together as one or more Batch CD-Search jobs
    (chunked at CDD_MAX_QUERIES_PER_BATCH sequences per job, NCBI's
    documented per-request limit) instead of one job per protein, and each
    job is polled to completion (see analysis/cdd.py::poll_cdd_batch)
    before its results are distributed back to the matching records. If a
    chunk's job fails outright (submission error, terminal status, or
    polling timeout), every record in that chunk gets the same "could not
    evaluate" failure recorded (note + AnnotationResult(success=False)) as
    a single failed record would via annotate_cdd_domains -- it is never
    silently treated as zero domain hits.
    """
    pending: list[ProteinRecord] = []
    for record in records.values():
        if cache is not None and cache.has("cdd", record.protein_id):
            _apply_cdd_success(record, _cached_cdd_domains(record.protein_id, cache))
            continue
        if not record.sequence.strip():
            _apply_cdd_success(record, [])
            continue
        pending.append(record)

    for chunk in _chunked(pending, CDD_MAX_QUERIES_PER_BATCH):
        _run_cdd_batch(
            chunk,
            cache=cache,
            timeout=timeout,
            poll_interval=poll_interval,
            max_wait=max_wait,
        )

    return records


def _run_cdd_batch(
    chunk: list[ProteinRecord],
    cache: JsonCache | None,
    timeout: int,
    poll_interval: float,
    max_wait: float,
) -> None:
    """Submit, poll, and distribute the results of one CDD batch job."""
    queries = [(record.protein_id, record.sequence) for record in chunk]

    try:
        cdsid = submit_cdd_batch(queries, timeout=timeout)
        poll_cdd_batch(cdsid, timeout=timeout, poll_interval=poll_interval, max_wait=max_wait)
        hits_by_id = parse_cdd_batch_response(fetch_cdd_batch_results(cdsid, timeout=timeout))
    except Exception as exc:
        for record in chunk:
            _apply_cdd_failure(record, exc)
        return

    for record in chunk:
        domains = hits_by_id.get(record.protein_id, [])
        _apply_cdd_success(record, domains)
        if cache is not None:
            cache.set("cdd", record.protein_id, [domain_hit_to_dict(hit) for hit in domains])


def _cached_cdd_domains(protein_id: str, cache: JsonCache) -> list[DomainHit]:
    """Return previously cached CDD domain hits for one protein, if any."""
    cached = cache.get("cdd", protein_id)
    if isinstance(cached, list):
        return [domain_hit_from_dict(item) for item in cached if isinstance(item, dict)]
    return []


def _apply_cdd_success(record: ProteinRecord, domains: list[DomainHit]) -> None:
    """Record a successful CDD annotation outcome, empty or not."""
    record.domains.extend(domains)
    record.annotations["cdd"] = AnnotationResult(
        protein_id=record.protein_id,
        source="cdd",
        success=True,
        domains=list(domains),
        metadata={"domain_count": len(domains)},
    )


def _apply_cdd_failure(record: ProteinRecord, exc: Exception) -> None:
    """Record a CDD annotation that could not be evaluated (never a phantom zero)."""
    message = f"CDD annotation failed for {record.protein_id}: {exc}"
    record.notes.append(message)
    record.annotations["cdd"] = AnnotationResult(
        protein_id=record.protein_id,
        source="cdd",
        success=False,
        error=message,
    )


def _chunked(items: list[ProteinRecord], size: int) -> list[list[ProteinRecord]]:
    """Split records into chunks of at most `size`, preserving order."""
    return [items[index : index + size] for index in range(0, len(items), size)]


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
    domains = enrich_pfam_domains_with_metadata(
        domains,
        cache=cache,
        timeout=timeout,
    )
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
