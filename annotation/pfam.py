"""Pfam annotation helpers for ProteinHunter."""

from __future__ import annotations

import csv
import json
import re
import time
from collections.abc import Iterable
from io import StringIO

import requests

from core.cache import JsonCache
from core.exceptions import PfamAnnotationError
from core.models import DomainHit


PFAM_SEARCH_URL = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan"
PFAM_RESULT_URL = "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result"
PFAM_REQUEST_HEADERS = {"Accept": "application/json"}
PFAM_JSON_HEADERS = {
    "Accept": "application/json",
    "Content-Type": "application/json",
}
PFAM_POLL_ATTEMPTS = 5


def parse_pfam_response(text: str) -> list[DomainHit]:
    """Parse a tolerant Pfam or HMMER-like response into domain hits."""
    hits: list[DomainHit] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()

        if not line or line.startswith("#") or line.startswith(">"):
            continue

        hit = _parse_pfam_line(line)
        if hit is not None:
            hits.append(hit)

    return hits


def search_pfam_by_sequence(
    protein_id: str,
    sequence: str,
    cache: JsonCache | None = None,
    timeout: int = 60,
) -> list[DomainHit]:
    """Search Pfam/HMMER for one protein sequence and return domain hits."""
    if cache is not None and cache.has("pfam", protein_id):
        cached = cache.get("pfam", protein_id)
        if isinstance(cached, list):
            return [
                domain_hit_from_dict(item)
                for item in cached
                if isinstance(item, dict)
            ]

    if not sequence.strip():
        return []

    payload = {
        "input": sequence.strip(),
        "database": "pfam",
    }

    response, body_label = _post_pfam_search(
        protein_id=protein_id,
        payload=payload,
        timeout=timeout,
    )

    status_code = _response_status_code(response)
    response_text = _response_text(response)

    if status_code is not None and status_code >= 400:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="request",
                detail=f"Pfam returned an HTTP error response after {body_label}.",
                status_code=status_code,
                response_text=response_text,
                timeout=timeout,
            )
        )

    try:
        response.raise_for_status()
    except requests.RequestException as exc:
        error_response = getattr(exc, "response", None) or response
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="request",
                detail=f"Pfam request failed after {body_label}: {exc}",
                status_code=_response_status_code(error_response),
                response_text=_response_text(error_response),
                timeout=timeout,
            )
        ) from exc

    try:
        hits = _parse_pfam_search_or_fetch_result(
            protein_id=protein_id,
            text=response_text,
            timeout=timeout,
        )
    except PfamAnnotationError:
        raise
    except Exception as exc:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="parse",
                detail=f"Pfam response could not be parsed: {exc}",
                status_code=status_code,
                response_text=response_text,
                timeout=timeout,
            )
        ) from exc

    if cache is not None:
        cache.set("pfam", protein_id, [domain_hit_to_dict(hit) for hit in hits])

    return hits


def _post_pfam_search(
    protein_id: str,
    payload: dict[str, str],
    timeout: int,
) -> tuple[requests.Response, str]:
    """Submit a Pfam search, trying JSON first and form data only when needed."""
    try:
        response = requests.post(
            PFAM_SEARCH_URL,
            json=payload,
            headers=PFAM_JSON_HEADERS,
            timeout=timeout,
        )
    except requests.Timeout as exc:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="request",
                detail=f"The JSON body request timed out after {timeout} seconds.",
                timeout=timeout,
            )
        ) from exc
    except requests.RequestException as exc:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="request",
                detail=f"The JSON body request could not be completed: {exc}",
                timeout=timeout,
            )
        ) from exc

    if not _should_retry_with_form_body(response):
        return response, "JSON body"

    try:
        form_response = requests.post(
            PFAM_SEARCH_URL,
            data=payload,
            headers=PFAM_REQUEST_HEADERS,
            timeout=timeout,
        )
    except requests.Timeout as exc:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="request",
                detail=(
                    "The JSON body was rejected as unparseable, then the form body "
                    f"request timed out after {timeout} seconds."
                ),
                status_code=_response_status_code(response),
                response_text=_response_text(response),
                timeout=timeout,
            )
        ) from exc
    except requests.RequestException as exc:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="request",
                detail=(
                    "The JSON body was rejected as unparseable, then the form body "
                    f"request could not be completed: {exc}"
                ),
                status_code=_response_status_code(response),
                response_text=_response_text(response),
                timeout=timeout,
            )
        ) from exc

    return form_response, "form body after JSON body retry"


def _should_retry_with_form_body(response: object) -> bool:
    """Return True when HMMER says the JSON request body could not be parsed."""
    status_code = _response_status_code(response)
    response_text = _response_text(response).lower()
    return status_code == 400 and "cannot parse request body" in response_text


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
        source=_string_or_default(data.get("source"), "Pfam"),
        accession=_string_or_default(data.get("accession"), ""),
        name=_string_or_default(data.get("name"), ""),
        description=_string_or_default(data.get("description"), ""),
        evalue=_optional_float(data.get("evalue")),
        bitscore=_optional_float(data.get("bitscore")),
        start=_optional_int(data.get("start")),
        end=_optional_int(data.get("end")),
    )


def _format_pfam_error(
    protein_id: str,
    phase: str,
    detail: str,
    status_code: int | None = None,
    response_text: str | None = None,
    timeout: int | None = None,
    endpoint: str = PFAM_SEARCH_URL,
) -> str:
    """Build a beginner-readable Pfam diagnostic message."""
    parts = [
        f"Pfam search failed for '{protein_id}' during the {phase} phase.",
        f"Endpoint: {endpoint}",
    ]

    if status_code is not None:
        parts.append(f"HTTP status: {status_code}")

    if timeout is not None:
        parts.append(f"Timeout setting: {timeout} seconds")

    if response_text is not None:
        parts.append(f"Response preview: {_response_preview(response_text)}")

    parts.append(f"Details: {detail}")
    return " ".join(parts)


def _response_status_code(response: object) -> int | None:
    """Return an HTTP status code when it is available."""
    status_code = getattr(response, "status_code", None)

    if isinstance(status_code, int):
        return status_code

    return None


def _response_text(response: object) -> str:
    """Return response text as a safe string."""
    text = getattr(response, "text", "")

    if isinstance(text, str):
        return text

    return str(text)


def _response_preview(text: str, limit: int = 300) -> str:
    """Return a short one-line response preview."""
    preview = " ".join(text.split())

    if len(preview) > limit:
        return preview[:limit] + "..."

    return preview


def _parse_pfam_search_or_fetch_result(
    protein_id: str,
    text: str,
    timeout: int,
) -> list[DomainHit]:
    """Parse direct search results or fetch an asynchronous HMMER result."""
    stripped = text.lstrip()

    if not stripped:
        return []

    if not stripped.startswith(("{", "[")):
        return parse_pfam_response(text)

    try:
        data = json.loads(stripped)
    except json.JSONDecodeError as exc:
        raise ValueError("Pfam returned invalid JSON instead of text output.") from exc

    try:
        return _parse_pfam_json_response(data)
    except ValueError as exc:
        job_id = _find_result_job_id(data)
        if job_id is None:
            raise exc

    return _poll_pfam_result(
        protein_id=protein_id,
        job_id=job_id,
        timeout=timeout,
    )


def _poll_pfam_result(
    protein_id: str,
    job_id: str,
    timeout: int,
) -> list[DomainHit]:
    """Poll the HMMER result endpoint until the job succeeds or fails."""
    endpoint = _result_endpoint(job_id)
    last_status = "UNKNOWN"
    last_response_text = ""
    last_status_code: int | None = None

    for attempt in range(PFAM_POLL_ATTEMPTS):
        response = _get_pfam_result(
            protein_id=protein_id,
            endpoint=endpoint,
            timeout=timeout,
        )
        last_status_code = _response_status_code(response)
        last_response_text = _response_text(response)

        if last_status_code is not None and last_status_code >= 400:
            raise PfamAnnotationError(
                _format_pfam_error(
                    protein_id=protein_id,
                    phase="result",
                    detail="Pfam result retrieval returned an HTTP error response.",
                    status_code=last_status_code,
                    response_text=last_response_text,
                    timeout=timeout,
                    endpoint=endpoint,
                )
            )

        try:
            response.raise_for_status()
        except requests.RequestException as exc:
            error_response = getattr(exc, "response", None) or response
            raise PfamAnnotationError(
                _format_pfam_error(
                    protein_id=protein_id,
                    phase="result",
                    detail=f"Pfam result retrieval failed: {exc}",
                    status_code=_response_status_code(error_response),
                    response_text=_response_text(error_response),
                    timeout=timeout,
                    endpoint=endpoint,
                )
            ) from exc

        stripped = last_response_text.lstrip()
        if not stripped.startswith(("{", "[")):
            return parse_pfam_response(last_response_text)

        try:
            data = json.loads(stripped)
        except json.JSONDecodeError as exc:
            raise PfamAnnotationError(
                _format_pfam_error(
                    protein_id=protein_id,
                    phase="parse",
                    detail=f"Pfam result JSON could not be parsed: {exc}",
                    status_code=last_status_code,
                    response_text=last_response_text,
                    timeout=timeout,
                    endpoint=endpoint,
                )
            ) from exc

        status = _find_json_status(data)
        if status is not None:
            last_status = status

        if _is_failure_status(status):
            raise PfamAnnotationError(
                _format_pfam_error(
                    protein_id=protein_id,
                    phase="result",
                    detail=f"Pfam result job ended with status {status}.",
                    status_code=last_status_code,
                    response_text=last_response_text,
                    timeout=timeout,
                    endpoint=endpoint,
                )
            )

        if _is_pending_status(status):
            if attempt < PFAM_POLL_ATTEMPTS - 1:
                time.sleep(_poll_sleep_seconds(timeout))
                continue

            break

        try:
            return _parse_pfam_json_response(data)
        except ValueError as exc:
            raise PfamAnnotationError(
                _format_pfam_error(
                    protein_id=protein_id,
                    phase="parse",
                    detail=f"Pfam result response could not be parsed: {exc}",
                    status_code=last_status_code,
                    response_text=last_response_text,
                    timeout=timeout,
                    endpoint=endpoint,
                )
            ) from exc

    raise PfamAnnotationError(
        _format_pfam_error(
            protein_id=protein_id,
            phase="result",
            detail=(
                "Pfam result job did not finish before polling stopped. "
                f"Last status: {last_status}."
            ),
            status_code=last_status_code,
            response_text=last_response_text,
            timeout=timeout,
            endpoint=endpoint,
        )
    )


def _get_pfam_result(
    protein_id: str,
    endpoint: str,
    timeout: int,
) -> requests.Response:
    """Fetch one HMMER result response."""
    try:
        return requests.get(
            endpoint,
            headers=PFAM_REQUEST_HEADERS,
            timeout=timeout,
        )
    except requests.Timeout as exc:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="result",
                detail=f"The result request timed out after {timeout} seconds.",
                timeout=timeout,
                endpoint=endpoint,
            )
        ) from exc
    except requests.RequestException as exc:
        raise PfamAnnotationError(
            _format_pfam_error(
                protein_id=protein_id,
                phase="result",
                detail=f"The result request could not be completed: {exc}",
                timeout=timeout,
                endpoint=endpoint,
            )
        ) from exc


def _result_endpoint(job_id: str) -> str:
    """Build the HMMER result endpoint for a job id."""
    return f"{PFAM_RESULT_URL}/{job_id}"


def _poll_sleep_seconds(timeout: int) -> float:
    """Return a short pause between result polling attempts."""
    if timeout <= 0:
        return 0.0

    return min(1.0, timeout / PFAM_POLL_ATTEMPTS)


def _parse_pfam_search_response(text: str) -> list[DomainHit]:
    """Parse either HMMER API JSON or older text/tabular output."""
    stripped = text.lstrip()

    if not stripped:
        return []

    if stripped.startswith(("{", "[")):
        try:
            data = json.loads(stripped)
        except json.JSONDecodeError as exc:
            raise ValueError("Pfam returned invalid JSON instead of text output.") from exc

        return _parse_pfam_json_response(data)

    return parse_pfam_response(text)


def _parse_pfam_json_response(data: object) -> list[DomainHit]:
    """Parse common HMMER API JSON shapes into domain hits."""
    hits = _collect_json_domain_hits(data)

    if hits:
        return hits

    if _json_clearly_has_no_hits(data):
        return []

    keys = _json_top_level_keys(data)
    raise ValueError(
        "Pfam returned JSON, but no Pfam domain fields were recognized. "
        f"Available top-level keys: {keys}."
    )


def _collect_json_domain_hits(
    value: object,
    parent: dict[str, object] | None = None,
) -> list[DomainHit]:
    """Recursively collect Pfam hits from common HMMER JSON containers."""
    hits: list[DomainHit] = []

    if isinstance(value, dict):
        merged = _merge_hit_context(parent, value)
        hit = _domain_hit_from_json_dict(merged)
        child_hits = _collect_json_children(value, merged)

        if child_hits:
            hits.extend(child_hits)
        elif hit is not None:
            hits.append(hit)

        return hits

    if isinstance(value, list):
        for item in value:
            hits.extend(_collect_json_domain_hits(item, parent))

    return hits


def _collect_json_children(
    data: dict[str, object],
    parent: dict[str, object],
) -> list[DomainHit]:
    """Collect hits from nested result containers without assuming one schema."""
    hits: list[DomainHit] = []

    for key, value in data.items():
        if key.lower() in {"results", "hits", "domains", "matches"}:
            hits.extend(_collect_json_domain_hits(value, parent))

    return hits


def _merge_hit_context(
    parent: dict[str, object] | None,
    child: dict[str, object],
) -> dict[str, object]:
    """Merge parent hit metadata with child domain metadata."""
    if parent is None:
        return dict(child)

    merged = dict(parent)
    merged.update(child)
    return merged


def _domain_hit_from_json_dict(data: dict[str, object]) -> DomainHit | None:
    """Build a DomainHit when a JSON object contains Pfam-like fields."""
    accession = _find_json_accession(data)
    if accession is None:
        return None

    return DomainHit(
        source="Pfam",
        accession=accession,
        name=_find_json_string(
            data,
            ("name", "target_name", "hmm_name", "model_name", "id"),
            accession,
        ),
        description=_find_json_string(data, ("description", "desc", "summary"), ""),
        evalue=_find_json_float(
            data,
            (
                "evalue",
                "e_value",
                "i_evalue",
                "ievalue",
                "independent_evalue",
                "conditional_evalue",
            ),
        ),
        bitscore=_find_json_float(data, ("bitscore", "bit_score", "score", "domain_score")),
        start=_find_json_int(data, ("start", "from", "ali_from", "env_from", "hmm_from")),
        end=_find_json_int(data, ("end", "to", "ali_to", "env_to", "hmm_to")),
    )


def _find_json_accession(data: dict[str, object]) -> str | None:
    """Find a Pfam accession in common accession fields or string values."""
    for key in ("accession", "acc", "target_accession", "hmm_acc", "model_acc"):
        value = data.get(key)
        if isinstance(value, str):
            accession = _clean_accession(value)
            if accession.upper().startswith("PF"):
                return accession

    for value in data.values():
        if isinstance(value, str):
            match = re.search(r"PF[0-9]{5}(?:[.][0-9]+)?", value, re.IGNORECASE)
            if match:
                return match.group(0)

    return None


def _find_json_string(
    data: dict[str, object],
    keys: Iterable[str],
    default: str,
) -> str:
    """Return the first useful string from a JSON object."""
    for key in keys:
        value = data.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()

    return default


def _find_json_float(data: dict[str, object], keys: Iterable[str]) -> float | None:
    """Return the first useful float from a JSON object."""
    for key in keys:
        value = _optional_float(data.get(key))
        if value is not None:
            return value

    return None


def _find_json_int(data: dict[str, object], keys: Iterable[str]) -> int | None:
    """Return the first useful int from a JSON object."""
    for key in keys:
        value = _optional_int(data.get(key))
        if value is not None:
            return value

    return None


def _json_clearly_has_no_hits(data: object) -> bool:
    """Return True only when a JSON response clearly reports empty results."""
    if isinstance(data, list):
        return len(data) == 0

    if not isinstance(data, dict):
        return False

    for key in ("results", "hits", "domains", "matches"):
        value = data.get(key)
        if value == []:
            return True
        if isinstance(value, dict) and _json_clearly_has_no_hits(value):
            return True

    return False


def _json_top_level_keys(data: object) -> str:
    """Return readable top-level JSON keys for parse diagnostics."""
    if isinstance(data, dict):
        keys = sorted(str(key) for key in data.keys())
        return ", ".join(keys) if keys else "(none)"

    if isinstance(data, list):
        return "(JSON list)"

    return f"({type(data).__name__})"


def _find_result_job_id(data: object) -> str | None:
    """Find a HMMER asynchronous result id in common JSON fields."""
    if not isinstance(data, dict):
        return None

    for key in ("id", "job_id", "jobId", "uuid", "result_id", "resultId"):
        value = data.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip()

    for key in ("job", "result", "location"):
        value = data.get(key)
        if isinstance(value, str) and value.strip():
            return value.rstrip("/").split("/")[-1]

    return None


def _find_json_status(data: object) -> str | None:
    """Find a normalized job status in a JSON response."""
    if not isinstance(data, dict):
        return None

    for key in ("status", "state", "job_status", "jobStatus"):
        value = data.get(key)
        if isinstance(value, str) and value.strip():
            return value.strip().upper()

    return None


def _is_pending_status(status: str | None) -> bool:
    """Return True for statuses that mean the result is not ready yet."""
    return status in {"PENDING", "RUNNING", "QUEUED", "SUBMITTED"}


def _is_failure_status(status: str | None) -> bool:
    """Return True for statuses that mean the job failed."""
    return status in {"FAILURE", "FAILED", "ERROR"}


def _parse_pfam_line(line: str) -> DomainHit | None:
    """Parse one non-comment Pfam/HMMER-like line."""
    parts = _split_line(line)
    if len(parts) < 2:
        return None

    accession_index = _find_accession_index(parts)
    if accession_index is None:
        return None

    accession = _clean_accession(parts[accession_index])
    name = _find_name(parts, accession_index, accession)
    evalue = _find_labeled_float(parts, ("evalue", "e-value", "eval", "i-evalue"))
    bitscore = _find_labeled_float(parts, ("bitscore", "bit_score", "score"))
    start = _find_labeled_int(parts, ("start", "from", "ali_from", "env_from"))
    end = _find_labeled_int(parts, ("end", "to", "ali_to", "env_to"))

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
        source="Pfam",
        accession=accession,
        name=name,
        description=description,
        evalue=evalue,
        bitscore=bitscore,
        start=start,
        end=end,
    )


def _split_line(line: str) -> list[str]:
    """Split Pfam-like text using tabs, commas, or whitespace."""
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
    """Find the first Pfam-like accession field."""
    for index, part in enumerate(parts):
        if re.search(r"\bPF\d{5}(?:\.\d+)?\b", part, re.IGNORECASE):
            return index

    return None


def _clean_accession(value: str) -> str:
    """Return the Pfam accession without surrounding punctuation."""
    match = re.search(r"\bPF\d{5}(?:\.\d+)?\b", value, re.IGNORECASE)
    if match:
        return match.group(0)

    return value.strip()


def _find_name(parts: list[str], accession_index: int, default: str) -> str:
    """Find a readable Pfam name near the accession."""
    if accession_index > 0 and not _looks_like_value(parts[accession_index - 1]):
        return parts[accession_index - 1]

    if accession_index + 1 < len(parts):
        return parts[accession_index + 1]

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


def _looks_like_value(part: str) -> bool:
    """Return True for fields that look numeric, labeled, or coordinate-like."""
    return (
        _looks_like_labeled_value(part)
        or _optional_float(part) is not None
        or _find_coordinate_pair([part]) is not None
    )


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
    "PFAM_SEARCH_URL",
    "domain_hit_from_dict",
    "domain_hit_to_dict",
    "parse_pfam_response",
    "search_pfam_by_sequence",
)
