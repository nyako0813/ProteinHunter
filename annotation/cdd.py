"""CDD annotation helpers for ProteinHunter."""

from __future__ import annotations

import csv
import re
import time
from io import StringIO
from typing import Callable

import requests

from core.cache import JsonCache
from core.exceptions import CDDAnnotationError
from core.models import DomainHit


CDD_SEARCH_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"

#: Maximum protein sequences/identifiers submitted in one Batch CD-Search
#: request. NCBI's documentation (cdd_help.shtml, "Batch CD-Search Help >
#: Input > Maximum input") states "1,000 protein sequences and/or
#: identifiers", but a live submission of exactly 1000 queries was
#: empirically rejected with "#status 2 msg Too many queries. Please
#: submit 1000 or less queries per request" (observed while extending CDD
#: annotation to interaction_scoring's larger candidate pools -- see
#: docs/implementation_plan_sequence_evidence.md's CDD investigation
#: notes). NCBI's real, enforced limit is stricter than its own documented
#: one, so this stays one below the documented number as a safety margin
#: rather than trusting "1,000 or less" literally.
CDD_MAX_QUERIES_PER_BATCH = 999

#: Seconds to wait between job-status checks. This matches NCBI's own
#: sample client (bwrpsb.pl, "checking for completion, wait 5 seconds
#: between checks" -> sleep(5)) -- not a value we invented.
CDD_POLL_INTERVAL_SECONDS = 5.0

#: Maximum total seconds to keep polling before giving up. NCBI's own
#: sample client polls indefinitely and documents no maximum wait time, so
#: this bound is our own safety margin (not an NCBI-recommended value) to
#: keep a stuck job from hanging the pipeline forever.
CDD_MAX_POLL_SECONDS = 600.0

#: How many times to retry one status check after a transport-level
#: failure (connection error, timeout, ...) before giving up on it.
#: Observed transient failure rates (~a few percent of individual requests
#: in this environment, see the single-request-per-protein era of CDD
#: annotation before batching) make a single unretried attempt too fragile
#: once one status check runs dozens of times over a large batch's
#: lifetime -- one blip used to fail the entire (up to
#: CDD_MAX_QUERIES_PER_BATCH-sequence) batch. Not an NCBI-documented value.
CDD_NETWORK_RETRY_ATTEMPTS = 3

#: Seconds to wait between retry attempts after a transport-level failure.
#: Short and separate from CDD_POLL_INTERVAL_SECONDS -- this is recovering
#: from a connectivity blip, not waiting for the job itself to progress.
CDD_NETWORK_RETRY_INTERVAL_SECONDS = 2.0

#: Batch CD-Search job status codes, verbatim from cdd_help.shtml
#: ("Batch CD-Search Help > Scripted Data Downloads (Web API) > Check
#: status > job status codes"). "0" (success) and "3" (still running) are
#: handled separately by the polling loop below.
CDD_STATUS_MESSAGES: dict[str, str] = {
    "1": "Invalid search ID",
    "2": "No effective input (usually no query proteins or search ID specified)",
    "4": "Queue manager (qman) service error",
    "5": "Data is corrupted or no longer available (cache cleaned, etc)",
}

_CDSID_PATTERN = re.compile(r"^#cdsid\s+(\S+)", re.MULTILINE)
_STATUS_PATTERN = re.compile(r"^#status\s+(\d)", re.MULTILINE)
_QUERY_COLUMN_PATTERN = re.compile(r"^Q#\d+\s*-\s*>(\S+)")


def parse_cdd_response(text: str) -> list[DomainHit]:
    """Parse a tolerant CDD-like text response into domain hits."""
    hits: list[DomainHit] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if not line or line.startswith("#") or line.startswith(">"):
            continue

        hit = _parse_cdd_line(line)
        if hit is not None:
            hits.append(hit)

    return hits


def search_cdd_by_sequence(
    protein_id: str,
    sequence: str,
    cache: JsonCache | None = None,
    timeout: int = 60,
) -> list[DomainHit]:
    """Search CDD for one protein sequence and return domain hits."""
    if cache is not None and cache.has("cdd", protein_id):
        cached = cache.get("cdd", protein_id)
        if isinstance(cached, list):
            return [
                domain_hit_from_dict(item)
                for item in cached
                if isinstance(item, dict)
            ]

    if not sequence.strip():
        return []

    fasta_text = f">{protein_id}\n{sequence.strip()}\n"
    data = {
        "queries": fasta_text,
        "tdata": "hits",
        "dmode": "rep",
    }

    try:
        response = requests.post(CDD_SEARCH_URL, data=data, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise CDDAnnotationError(
            f"CDD search failed for '{protein_id}'. Please check the network connection."
        ) from exc

    try:
        hits = parse_cdd_response(response.text)
    except Exception as exc:
        raise CDDAnnotationError(
            f"CDD returned an unexpected response for '{protein_id}'."
        ) from exc

    if cache is not None:
        cache.set("cdd", protein_id, [domain_hit_to_dict(hit) for hit in hits])

    return hits


def submit_cdd_batch(queries: list[tuple[str, str]], timeout: int = 60) -> str:
    """Submit up to CDD_MAX_QUERIES_PER_BATCH (protein_id, sequence) pairs as one job.

    Returns the job's search ID (``cdsid``), used to poll for completion and
    retrieve results. Raises :class:`CDDAnnotationError` if the request
    fails or the response does not contain a search ID.
    """
    if not queries:
        raise CDDAnnotationError("CDD batch submission requires at least one query.")
    if len(queries) > CDD_MAX_QUERIES_PER_BATCH:
        raise CDDAnnotationError(
            f"CDD batch submission exceeds the {CDD_MAX_QUERIES_PER_BATCH}-sequence "
            f"limit ({len(queries)} given)."
        )

    fasta_text = "".join(
        f">{protein_id}\n{sequence.strip()}\n" for protein_id, sequence in queries
    )
    data = {"queries": fasta_text, "tdata": "hits", "dmode": "rep"}

    try:
        response = requests.post(CDD_SEARCH_URL, data=data, timeout=timeout)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise CDDAnnotationError(
            "CDD batch submission failed. Please check the network connection."
        ) from exc

    match = _CDSID_PATTERN.search(response.text)
    if match is None:
        raise CDDAnnotationError(
            "CDD batch submission did not return a search ID. "
            f"Response: {response.text[:200]!r}"
        )
    return match.group(1)


def poll_cdd_batch(
    cdsid: str,
    timeout: int = 60,
    poll_interval: float = CDD_POLL_INTERVAL_SECONDS,
    max_wait: float = CDD_MAX_POLL_SECONDS,
    sleep_fn: Callable[[float], None] = time.sleep,
    network_retry_attempts: int = CDD_NETWORK_RETRY_ATTEMPTS,
    network_retry_interval: float = CDD_NETWORK_RETRY_INTERVAL_SECONDS,
) -> None:
    """Poll a submitted batch job until it completes (status 0).

    Follows the status codes documented in cdd_help.shtml and NCBI's own
    sample client's polling loop (bwrpsb.pl). Raises
    :class:`CDDAnnotationError` on a terminal error status (1/2/4/5), a
    malformed status response, or exceeding ``max_wait`` while still
    "running" (status 3) -- NCBI's own client has no such bound, so this is
    our own safety timeout, not an NCBI-mandated value.

    A transport-level failure (connection error, timeout, ...) on one
    status check is retried up to ``network_retry_attempts`` times,
    ``network_retry_interval`` seconds apart, before being raised as a
    hard failure. This is deliberately narrower than the job-status
    handling above: NCBI's own status codes 1/2/4/5 mean the server
    responded and reported a real failure, so retrying them would not
    help and they are never retried. Retry backoff sleeps are tracked
    separately from ``elapsed`` and never count against ``max_wait`` --
    ``max_wait`` bounds how long we wait for the job itself to finish, not
    how long we spend working around flaky connectivity.
    """
    elapsed = 0.0
    while True:
        sleep_fn(poll_interval)
        elapsed += poll_interval

        response = _post_with_network_retries(
            data={"tdata": "hits", "cdsid": cdsid},
            timeout=timeout,
            cdsid=cdsid,
            attempts=network_retry_attempts,
            retry_interval=network_retry_interval,
            sleep_fn=sleep_fn,
        )

        match = _STATUS_PATTERN.search(response.text)
        if match is None:
            raise CDDAnnotationError(
                f"CDD status check returned an unexpected response for search "
                f"'{cdsid}': {response.text[:200]!r}"
            )

        status = match.group(1)
        if status == "0":
            return
        if status in CDD_STATUS_MESSAGES:
            raise CDDAnnotationError(
                f"CDD batch search '{cdsid}' failed: {CDD_STATUS_MESSAGES[status]} "
                f"(status={status})."
            )
        # status == "3": still running/waiting -- keep polling.
        if elapsed >= max_wait:
            raise CDDAnnotationError(
                f"CDD batch search '{cdsid}' did not complete within "
                f"{max_wait:.0f} seconds."
            )


def _post_with_network_retries(
    *,
    data: dict[str, str],
    timeout: int,
    cdsid: str,
    attempts: int,
    retry_interval: float,
    sleep_fn: Callable[[float], None],
) -> requests.Response:
    """POST to the CDD endpoint, retrying only transport-level failures.

    A response that comes back at all -- even one reporting a job failure
    via NCBI's own status codes -- is returned as-is on the first try;
    only ``requests.RequestException`` (connection errors, timeouts, ...)
    triggers a retry, since those mean no response was received at all.
    """
    last_exc: requests.RequestException | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = requests.post(CDD_SEARCH_URL, data=data, timeout=timeout)
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            last_exc = exc
            if attempt < attempts:
                sleep_fn(retry_interval)

    raise CDDAnnotationError(
        f"CDD status check failed for search '{cdsid}' after {attempts} attempts. "
        "Please check the network connection."
    ) from last_exc


def fetch_cdd_batch_results(cdsid: str, timeout: int = 60) -> str:
    """Retrieve the raw domain-hit results text for a completed batch job."""
    try:
        response = requests.post(
            CDD_SEARCH_URL, data={"tdata": "hits", "cdsid": cdsid}, timeout=timeout
        )
        response.raise_for_status()
    except requests.RequestException as exc:
        raise CDDAnnotationError(
            f"CDD result retrieval failed for search '{cdsid}'. "
            "Please check the network connection."
        ) from exc
    return response.text


def parse_cdd_batch_response(text: str) -> dict[str, list[DomainHit]]:
    """Parse a completed Batch CD-Search response into per-query domain hits.

    Each data row's ``Query`` column has the form ``Q#<index> -
    ><protein_id>`` (observed directly from a live Batch CD-Search job;
    see docs/implementation_plan_sequence_evidence.md's CDD investigation
    notes -- this is NOT the same line shape ``parse_cdd_response`` expects,
    which never reflected real Batch CD-Search "Concise Results" output).
    Column positions are read from the response's own header row rather
    than hard-coded, so a future harmless column reorder does not silently
    break parsing.
    """
    hits_by_query: dict[str, list[DomainHit]] = {}
    column_index: dict[str, int] | None = None

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        columns = line.split("\t")
        if columns[0] == "Query":
            column_index = {name.strip(): index for index, name in enumerate(columns)}
            continue
        if column_index is None:
            continue

        query_match = _QUERY_COLUMN_PATTERN.match(columns[0])
        if query_match is None:
            continue
        protein_id = query_match.group(1)

        hit = _domain_hit_from_columns(columns, column_index)
        if hit is not None:
            hits_by_query.setdefault(protein_id, []).append(hit)

    return hits_by_query


def _domain_hit_from_columns(
    columns: list[str], column_index: dict[str, int]
) -> DomainHit | None:
    """Build one DomainHit from a Concise Results data row and its header map."""
    accession = _column_value(columns, column_index, "Accession")
    if not accession:
        return None

    return DomainHit(
        source="CDD",
        accession=accession,
        name=_column_value(columns, column_index, "Short name") or accession,
        description=_column_value(columns, column_index, "Superfamily"),
        evalue=_optional_float(_column_value(columns, column_index, "E-Value")),
        bitscore=_optional_float(_column_value(columns, column_index, "Bitscore")),
        start=_optional_int(_column_value(columns, column_index, "From")),
        end=_optional_int(_column_value(columns, column_index, "To")),
    )


def _column_value(
    columns: list[str], column_index: dict[str, int], name: str
) -> str | None:
    index = column_index.get(name)
    if index is None or index >= len(columns):
        return None
    value = columns[index].strip()
    return value or None


def domain_hit_to_dict(hit: DomainHit) -> dict[str, str | int | float | None]:
    """Convert a DomainHit to a JSON-serializable dictionary."""
    return {
        "source": hit.source,
        "accession": hit.accession,
        "name": hit.name,
        "description": hit.description,
        "evalue": hit.evalue,
        "bitscore": hit.bitscore,
        "start": hit.start,
        "end": hit.end,
    }


def domain_hit_from_dict(data: dict[str, object]) -> DomainHit:
    """Convert cached domain data back into a DomainHit."""
    return DomainHit(
        source=_string_or_default(data.get("source"), "CDD"),
        accession=_string_or_default(data.get("accession"), ""),
        name=_string_or_default(data.get("name"), ""),
        description=_string_or_default(data.get("description"), ""),
        evalue=_optional_float(data.get("evalue")),
        bitscore=_optional_float(data.get("bitscore")),
        start=_optional_int(data.get("start")),
        end=_optional_int(data.get("end")),
    )


def _parse_cdd_line(line: str) -> DomainHit | None:
    """Parse one non-comment CDD-like line."""
    parts = _split_line(line)
    if len(parts) < 2:
        return None

    accession_index = _find_accession_index(parts)
    if accession_index is None:
        return None

    accession = parts[accession_index]
    name = _field_after(parts, accession_index, default=accession)
    evalue = _find_labeled_float(parts, ("evalue", "e-value", "eval"))
    bitscore = _find_labeled_float(parts, ("bitscore", "bit_score", "score"))
    start = _find_labeled_int(parts, ("start", "from"))
    end = _find_labeled_int(parts, ("end", "to"))

    if evalue is None:
        evalue = _first_float_after(parts, accession_index)

    if bitscore is None:
        bitscore = _second_float_after(parts, accession_index)

    if start is None or end is None:
        coordinates = _find_coordinate_pair(parts)
        if coordinates is not None:
            start, end = coordinates

    description = _description_from_parts(parts, accession_index)

    return DomainHit(
        source="CDD",
        accession=accession,
        name=name,
        description=description,
        evalue=evalue,
        bitscore=bitscore,
        start=start,
        end=end,
    )


def _split_line(line: str) -> list[str]:
    """Split CDD-like text using tabs, commas, or repeated whitespace."""
    if "\t" in line:
        return [part.strip() for part in line.split("\t") if part.strip()]

    try:
        csv_parts = next(csv.reader(StringIO(line)))
    except csv.Error:
        csv_parts = []

    if len(csv_parts) > 1:
        return [part.strip() for part in csv_parts if part.strip()]

    return [part.strip() for part in re.split(r"\s{2,}|\s+", line) if part.strip()]


def _find_accession_index(parts: list[str]) -> int | None:
    """Find the first CDD-like accession field."""
    for index, part in enumerate(parts):
        if re.search(r"\b(?:cd|CDD|pfam|smart|COG)\d+\b", part, re.IGNORECASE):
            return index

    return None


def _field_after(parts: list[str], index: int, default: str) -> str:
    """Return the next field after an index, or a fallback value."""
    if index + 1 < len(parts):
        return parts[index + 1]

    return default


def _description_from_parts(parts: list[str], accession_index: int) -> str:
    """Return a readable description when extra fields are available."""
    start_index = accession_index + 2
    description_parts = [
        part
        for part in parts[start_index:]
        if not _looks_like_labeled_value(part)
        and _optional_float(part) is None
        and _optional_int(part) is None
        and _find_coordinate_pair([part]) is None
    ]
    return " ".join(description_parts)


def _find_labeled_float(parts: list[str], labels: tuple[str, ...]) -> float | None:
    """Find a float from fields like evalue=1e-5 or score:42.0."""
    for part in parts:
        key, value = _split_labeled_value(part)
        if key in labels:
            return _optional_float(value)

    return None


def _find_labeled_int(parts: list[str], labels: tuple[str, ...]) -> int | None:
    """Find an int from fields like start=10 or end:80."""
    for part in parts:
        key, value = _split_labeled_value(part)
        if key in labels:
            return _optional_int(value)

    return None


def _split_labeled_value(part: str) -> tuple[str, str]:
    """Split a key-value field into lowercase key and raw value."""
    if "=" in part:
        key, value = part.split("=", 1)
    elif ":" in part:
        key, value = part.split(":", 1)
    else:
        return "", part

    return key.strip().lower(), value.strip()


def _looks_like_labeled_value(part: str) -> bool:
    """Return True for simple key=value or key:value fields."""
    return "=" in part or ":" in part


def _first_float_after(parts: list[str], index: int) -> float | None:
    """Return the first float after an index."""
    for part in parts[index + 1 :]:
        value = _optional_float(part)
        if value is not None:
            return value

    return None


def _second_float_after(parts: list[str], index: int) -> float | None:
    """Return the second float after an index."""
    seen = 0

    for part in parts[index + 1 :]:
        value = _optional_float(part)
        if value is None:
            continue

        seen += 1
        if seen == 2:
            return value

    return None


def _find_coordinate_pair(parts: list[str]) -> tuple[int, int] | None:
    """Find a simple coordinate field such as 10-80."""
    for part in parts:
        match = re.fullmatch(r"(\d+)\s*-\s*(\d+)", part)
        if match:
            return int(match.group(1)), int(match.group(2))

    return None


def _string_or_default(value: object, default: str) -> str:
    """Return a string value or a fallback."""
    if isinstance(value, str):
        return value

    return default


def _optional_float(value: object) -> float | None:
    """Return a float when the value can be parsed."""
    if value is None or isinstance(value, bool):
        return None

    try:
        return float(str(value).strip())
    except ValueError:
        return None


def _optional_int(value: object) -> int | None:
    """Return an int when the value can be parsed."""
    if value is None or isinstance(value, bool):
        return None

    try:
        return int(str(value).strip())
    except ValueError:
        return None


__all__: tuple[str, ...] = (
    "CDD_MAX_POLL_SECONDS",
    "CDD_MAX_QUERIES_PER_BATCH",
    "CDD_NETWORK_RETRY_ATTEMPTS",
    "CDD_NETWORK_RETRY_INTERVAL_SECONDS",
    "CDD_POLL_INTERVAL_SECONDS",
    "CDD_SEARCH_URL",
    "CDD_STATUS_MESSAGES",
    "domain_hit_from_dict",
    "domain_hit_to_dict",
    "fetch_cdd_batch_results",
    "parse_cdd_batch_response",
    "parse_cdd_response",
    "poll_cdd_batch",
    "search_cdd_by_sequence",
    "submit_cdd_batch",
)
