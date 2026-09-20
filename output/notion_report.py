"""Notion export: the report as Notion pages, a third output beside Excel and Word.

Renders the same renderer-neutral ``NarrativeSection`` list as the Word report
(``output/report_sections.py``), so wording, language (``report_language``) and
candidate selection are decided once. Nothing here changes any computed result.

Layout in the workspace (claude/notion_export_design.md):

    parent page (notion_export.parent_page_id, shared with the integration)
      └─ Run page          one per pipeline run: title-page lines, table of
           │               contents, the per-query ranking tables, and links to:
           ├─ evidence pages   one child page per evidence category (5.1-5.7)
           └─ Candidates       an inline database, one row per candidate
                └─ page        the candidate's "why it ranks highly",
                               "biological interpretation" and domain
                               information, opened from the row

Database properties carry the data behind each row (query, rank, final score,
tier, candidate source, candidate id) so the database can be filtered and
grouped in Notion; their values are the same data values the Excel workbook
has and are never translated.

The Notion API is called through ``notion-client``'s low-level ``request`` (not
its per-endpoint helpers, some of which drop fields depending on SDK version --
3.x silently omitted ``properties`` from ``databases.create``), pinned to API
version 2022-06-28, with the SDK's own automatic retries switched off: this module
throttles to Notion's ~3 requests/second average and retries 429/5xx/network
errors with exponential backoff itself (honouring ``Retry-After``), so the
behaviour is the same whichever SDK version is installed and is testable.
Requests are split to Notion's limits: at most 100 blocks per request (and a
conservative 900 including table rows), 2000 characters per rich-text object.
"""

from __future__ import annotations

import logging
import random
import re
import time
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Callable, Iterator, Sequence

from core.exceptions import NotionExportError
from output.report_i18n import DEFAULT_LANGUAGE, normalize_language, t
from output.word_narrative import NarrativeSection

NOTION_VERSION = "2022-06-28"
RICH_TEXT_LIMIT = 2000
MAX_CHILDREN_PER_REQUEST = 100
MAX_BLOCKS_PER_REQUEST = 900  # Notion's cap is 1000 block elements per request, nested rows included
MAX_TABLE_ROWS = 100
SELECT_NAME_LIMIT = 100
MAX_CONSECUTIVE_CANDIDATE_FAILURES = 3

_RETRYABLE_STATUSES = frozenset({409, 429, 500, 502, 503, 504})
_TRANSIENT_ERROR_NAMES = frozenset(
    {
        "RequestTimeoutError",
        "TimeoutException",
        "ConnectTimeout",
        "ReadTimeout",
        "WriteTimeout",
        "PoolTimeout",
        "ConnectError",
        "ReadError",
        "WriteError",
        "NetworkError",
        "RemoteProtocolError",
        "TransportError",
    }
)
_MAX_RETRY_AFTER_SECONDS = 120.0


# ---------------------------------------------------------------------------
# Sections -> Notion blocks
# ---------------------------------------------------------------------------


def rich_text(text: str, *, bold: bool = False) -> list[dict[str, Any]]:
    """Notion rich-text objects for ``text``, split so none exceeds 2000 characters."""
    items: list[dict[str, Any]] = []
    for start in range(0, len(text), RICH_TEXT_LIMIT):
        item: dict[str, Any] = {"type": "text", "text": {"content": text[start : start + RICH_TEXT_LIMIT]}}
        if bold:
            item["annotations"] = {"bold": True}
        items.append(item)
    return items


def _block(block_type: str, rich: list[dict[str, Any]]) -> dict[str, Any]:
    return {"object": "block", "type": block_type, block_type: {"rich_text": rich}}


def _table_blocks(header: Sequence[str], rows: Sequence[Sequence[str]]) -> list[dict[str, Any]]:
    """One table block per up to MAX_TABLE_ROWS rows (header repeated), with rows as nested children."""

    def row_block(cells: Sequence[str]) -> dict[str, Any]:
        return {"object": "block", "type": "table_row", "table_row": {"cells": [rich_text(cell) for cell in cells]}}

    body_capacity = MAX_TABLE_ROWS - 1
    chunks = [rows[i : i + body_capacity] for i in range(0, len(rows), body_capacity)] or [[]]
    return [
        {
            "object": "block",
            "type": "table",
            "table": {
                "table_width": len(header),
                "has_column_header": True,
                "has_row_header": False,
                "children": [row_block(header), *(row_block(row) for row in chunk)],
            },
        }
        for chunk in chunks
    ]


def section_to_blocks(section: NarrativeSection) -> list[dict[str, Any]]:
    """Notion blocks for one section. The document title (heading level 0) becomes the page title, not a block."""
    if section.kind == "heading":
        if section.level == 0:
            return []
        block_type = f"heading_{min(max(section.level, 1), 3)}"  # Notion has three heading levels
        return [_block(block_type, rich_text(section.text))]
    if section.kind == "paragraph":
        segments = (rich_text(section.label, bold=True) if section.label else []) + rich_text(section.text)
        return [_block("paragraph", segments)]
    if section.kind == "bullet_list":
        return [_block("bulleted_list_item", rich_text(item)) for item in section.items]
    if section.kind == "table":
        return _table_blocks(section.header, section.rows)
    if section.kind == "toc":
        return [{"object": "block", "type": "table_of_contents", "table_of_contents": {}}]
    raise ValueError(f"unknown report section kind: {section.kind!r}")


def count_blocks(block: dict[str, Any]) -> int:
    """The block plus everything nested under it (table rows), for Notion's per-request cap."""
    children = block.get(block.get("type", ""), {}).get("children", [])
    return 1 + sum(count_blocks(child) for child in children)


def batch_blocks(blocks: Sequence[dict[str, Any]]) -> Iterator[list[dict[str, Any]]]:
    """Group top-level blocks into requests within Notion's limits (100 children, 900 block elements)."""
    batch: list[dict[str, Any]] = []
    weight = 0
    for block in blocks:
        block_weight = count_blocks(block)
        if batch and (len(batch) >= MAX_CHILDREN_PER_REQUEST or weight + block_weight > MAX_BLOCKS_PER_REQUEST):
            yield batch
            batch, weight = [], 0
        batch.append(block)
        weight += block_weight
    if batch:
        yield batch


# ---------------------------------------------------------------------------
# Splitting the section list into the Run page and candidate pages
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CandidatePage:
    title: str
    meta: dict[str, str]
    blocks: list[dict[str, Any]]


@dataclass(frozen=True)
class EvidencePage:
    title: str
    blocks: list[dict[str, Any]]


@dataclass(frozen=True)
class ReportLayout:
    """The report split into Notion pages.

    ``run_head`` and ``run_tail`` are the Run page's blocks before and after
    the evidence child pages (which Notion lists at the point they are created,
    so the export writes head, then the evidence pages, then tail, then the
    database). Without evidence sections everything is in ``run_head``.
    """

    title: str
    run_head: list[dict[str, Any]]
    evidence_pages: list[EvidencePage]
    run_tail: list[dict[str, Any]]
    candidates: list[CandidatePage]

    @property
    def run_blocks(self) -> list[dict[str, Any]]:
        """All Run-page blocks in reading order (head then tail), ignoring the child pages between them."""
        return [*self.run_head, *self.run_tail]


def layout_report(sections: Sequence[NarrativeSection], language: str = DEFAULT_LANGUAGE) -> ReportLayout:
    """Split the report into the Run page, evidence child pages and one page per candidate.

    * Before the ``evidence_root`` heading: Run page (head).
    * The ``evidence_root`` heading and its intro paragraph: Run page (head),
      followed by a note that each category has its own page; every ``evidence``
      heading then opens a child page collecting the paragraphs after it.
    * Everything up to the ``details`` heading, that heading, and (when there
      are candidates) a note pointing at the database: Run page (tail). The
      ranking tables stay here as the overview.
    * Each ``candidate`` heading opens a candidate page that collects the
      sections after it (narrative, domain information, Excel reference);
      per-query headings inside the details are dropped (the query is a
      database property instead).

    Without those role markers the whole report stays on the Run page.
    """
    language = normalize_language(language)
    title = next((s.text for s in sections if s.kind == "heading" and s.level == 0), "ProteinHunter")
    head: list[dict[str, Any]] = []
    tail: list[dict[str, Any]] = []
    evidence_pages: list[EvidencePage] = []
    candidates: list[CandidatePage] = []
    current_evidence: tuple[str, list[dict[str, Any]]] | None = None
    current_candidate: tuple[str, dict[str, str], list[dict[str, Any]]] | None = None
    state = "head"  # head -> evidence (after the evidence_root heading) -> tail -> details
    details_note_index: int | None = None

    def close_evidence() -> None:
        nonlocal current_evidence
        if current_evidence is not None:
            evidence_pages.append(EvidencePage(title=current_evidence[0], blocks=current_evidence[1]))
            current_evidence = None

    def close_candidate() -> None:
        nonlocal current_candidate
        if current_candidate is not None:
            candidates.append(CandidatePage(title=current_candidate[0], meta=current_candidate[1], blocks=current_candidate[2]))
            current_candidate = None

    for section in sections:
        blocks = section_to_blocks(section)
        is_heading = section.kind == "heading"
        if is_heading and section.role == "evidence_root":
            head.extend(blocks)
            state = "evidence"
            continue
        if state == "evidence":
            if is_heading and section.role == "evidence":
                close_evidence()
                current_evidence = (section.text, [])
                continue
            if not is_heading and (current_evidence is not None or section.kind == "paragraph"):
                # Paragraphs after an evidence heading belong to its page; the paragraph right after the
                # root heading (the intro) stays on the Run page.
                (current_evidence[1] if current_evidence is not None else head).extend(blocks)
                continue
            close_evidence()
            if evidence_pages:
                head.append(_block("paragraph", rich_text(t("notion.evidence_note", language))))
            state = "tail"
        if is_heading and section.role == "details":
            close_evidence()
            state = "details"
            tail.extend(blocks)
            details_note_index = len(tail)
            continue
        if state == "details" and is_heading and section.role == "query":
            close_candidate()
            continue
        if state == "details" and is_heading and section.role == "candidate":
            close_candidate()
            current_candidate = (section.text, section.meta_dict(), [])
            continue
        if state == "details" and current_candidate is not None:
            current_candidate[2].extend(blocks)
        elif state == "head":
            head.extend(blocks)
        else:
            tail.extend(blocks)
    close_evidence()
    close_candidate()

    if state == "evidence" and evidence_pages:  # report ended inside the evidence section
        head.append(_block("paragraph", rich_text(t("notion.evidence_note", language))))
    if candidates and details_note_index is not None:
        tail.insert(details_note_index, _block("paragraph", rich_text(t("notion.details_note", language))))
    return ReportLayout(title=title, run_head=head, evidence_pages=evidence_pages, run_tail=tail, candidates=candidates)


# ---------------------------------------------------------------------------
# Database schema and row properties
# ---------------------------------------------------------------------------


def property_names(language: str) -> dict[str, str]:
    """Property (column) names, in the report language; the same labels as the ranking table."""
    return {
        "title": t("ranking.col.candidate", language),
        "query": t("notion.prop.query", language),
        "rank": t("ranking.col.rank", language),
        "final_score": t("ranking.col.final_score", language),
        "tier": t("ranking.col.tier", language),
        "source": t("ranking.col.source", language),
        "candidate_id": t("notion.prop.candidate_id", language),
    }


def database_properties(language: str) -> dict[str, Any]:
    names = property_names(language)
    return {
        names["title"]: {"title": {}},
        names["query"]: {"select": {}},
        names["rank"]: {"number": {"format": "number"}},
        names["final_score"]: {"number": {"format": "number"}},
        names["tier"]: {"select": {}},
        names["source"]: {"select": {}},
        names["candidate_id"]: {"rich_text": {}},
    }


def _select_name(value: str) -> str:
    """A legal select option name: Notion forbids commas and caps the length at 100."""
    return value.replace(",", " ").strip()[:SELECT_NAME_LIMIT]


def _number(value: str) -> float | int | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if number.is_integer() and "." not in value else number


def candidate_properties(page: CandidatePage, language: str) -> dict[str, Any]:
    """Row properties for one candidate page; empty or non-numeric data values are simply omitted."""
    names = property_names(language)
    meta = page.meta
    properties: dict[str, Any] = {names["title"]: {"title": rich_text(page.title)}}
    if meta.get("query_id"):
        properties[names["query"]] = {"select": {"name": _select_name(meta["query_id"])}}
    rank = _number(meta.get("rank", ""))
    if rank is not None:
        properties[names["rank"]] = {"number": rank}
    score = _number(meta.get("final_score", ""))
    if score is not None:
        properties[names["final_score"]] = {"number": score}
    if meta.get("tier"):
        properties[names["tier"]] = {"select": {"name": _select_name(meta["tier"])}}
    if meta.get("source"):
        properties[names["source"]] = {"select": {"name": _select_name(meta["source"])}}
    if meta.get("candidate_id"):
        properties[names["candidate_id"]] = {"rich_text": rich_text(meta["candidate_id"])}
    return properties


# ---------------------------------------------------------------------------
# API access: throttling and retries
# ---------------------------------------------------------------------------


def _describe(exc: Exception) -> str:
    status = getattr(exc, "status", None)
    code = getattr(exc, "code", None)
    detail = str(exc).strip().splitlines()[0][:200] if str(exc).strip() else ""
    parts = [type(exc).__name__]
    if status is not None:
        parts.append(f"status={status}")
    if code is not None:
        parts.append(f"code={getattr(code, 'value', code)}")
    return " ".join(parts) + (f": {detail}" if detail else "")


def classify_error(exc: Exception) -> tuple[bool, float | None]:
    """``(retryable, retry_after_seconds)`` for an exception raised by the Notion client."""
    status = getattr(exc, "status", None)
    if isinstance(status, int):
        if status not in _RETRYABLE_STATUSES:
            return False, None
        retry_after: float | None = None
        headers = getattr(exc, "headers", None)
        raw = headers.get("retry-after") if hasattr(headers, "get") else None
        if raw is not None:
            try:
                retry_after = min(max(float(raw), 0.0), _MAX_RETRY_AFTER_SECONDS)
            except ValueError:
                retry_after = None
        return True, retry_after
    if isinstance(exc, (TimeoutError, ConnectionError)):
        return True, None
    if any(cls.__name__ in _TRANSIENT_ERROR_NAMES for cls in type(exc).__mro__):
        return True, None
    return False, None


class NotionExporter:
    """Thin wrapper around a ``notion_client.Client`` that throttles and retries every call."""

    def __init__(
        self,
        client: Any,
        *,
        sleep: Callable[[float], None] = time.sleep,
        monotonic: Callable[[], float] = time.monotonic,
        jitter: Callable[[float], float] | None = None,
        max_retries: int = 6,
        base_delay: float = 1.0,
        max_delay: float = 30.0,
        min_interval: float = 0.34,
        log: Callable[[str], None] | None = None,
    ) -> None:
        self.client = client
        self._sleep = sleep
        self._monotonic = monotonic
        self._jitter = jitter if jitter is not None else (lambda delay: random.uniform(0.0, 0.25 * delay))
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.min_interval = min_interval
        self.log = log or (lambda message: None)
        self._last_call: float | None = None

    def _throttle(self) -> None:
        if self._last_call is not None:
            wait = self.min_interval - (self._monotonic() - self._last_call)
            if wait > 0:
                self._sleep(wait)
        self._last_call = self._monotonic()

    def backoff_delay(self, attempt: int) -> float:
        delay = min(self.max_delay, self.base_delay * (2**attempt))
        return delay + self._jitter(delay)

    def call(self, description: str, func: Callable[..., Any], **kwargs: Any) -> Any:
        """Run one API call, retrying rate-limit / server / network errors with exponential backoff."""
        for attempt in range(self.max_retries + 1):
            self._throttle()
            try:
                return func(**kwargs)
            except Exception as exc:  # noqa: BLE001 -- classified below; anything unexpected is wrapped
                retryable, retry_after = classify_error(exc)
                if not retryable or attempt == self.max_retries:
                    raise NotionExportError(f"Notion request failed ({description}): {_describe(exc)}") from exc
                delay = retry_after if retry_after is not None else self.backoff_delay(attempt)
                self.log(
                    f"Notion request '{description}' hit {_describe(exc)}; retrying in {delay:.1f}s "
                    f"(attempt {attempt + 1} of {self.max_retries})"
                )
                self._sleep(delay)
        raise AssertionError("unreachable")  # pragma: no cover

    # -- raw API calls: the body sent is exactly the body built here ----------
    # (``client.request`` rather than the SDK's per-endpoint helpers, which pick
    # a fixed set of fields per SDK version and dropped ``properties`` from
    # databases.create in notion-client 3.x.)

    def api(self, description: str, method: str, path: str, body: dict[str, Any]) -> Any:
        return self.call(description, self.client.request, path=path, method=method, body=body)

    def create_page(
        self,
        parent: dict[str, Any],
        properties: dict[str, Any],
        blocks: Sequence[dict[str, Any]],
        description: str,
    ) -> dict[str, Any]:
        """Create a page; blocks beyond the first request's limits are appended afterwards."""
        batches = list(batch_blocks(blocks))
        body: dict[str, Any] = {"parent": parent, "properties": properties}
        if batches:
            body["children"] = batches[0]
        page = self.api(description, "POST", "pages", body)
        for batch in batches[1:]:
            self.api(f"{description} (more blocks)", "PATCH", f"blocks/{page['id']}/children", {"children": batch})
        return page

    def append_blocks(self, block_id: str, blocks: Sequence[dict[str, Any]], description: str) -> None:
        for batch in batch_blocks(blocks):
            self.api(description, "PATCH", f"blocks/{block_id}/children", {"children": batch})

    def create_database(self, parent_page_id: str, title: str, properties: dict[str, Any], description: str) -> dict[str, Any]:
        return self.api(
            description,
            "POST",
            "databases",
            {
                "parent": {"type": "page_id", "page_id": parent_page_id},
                "title": rich_text(title),
                "is_inline": True,
                "properties": properties,
            },
        )


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NotionExportResult:
    run_page_id: str
    run_page_url: str
    database_id: str
    candidates_total: int
    candidates_created: int
    failed: tuple[str, ...] = field(default=())
    aborted: bool = False
    evidence_pages_total: int = 0
    evidence_pages_created: int = 0


def export_sections_to_notion(
    sections: Sequence[NarrativeSection],
    *,
    exporter: NotionExporter,
    parent_page_id: str,
    language: str = DEFAULT_LANGUAGE,
    run_name: str = "",
    generated_at: datetime | None = None,
) -> NotionExportResult:
    """Create the Run page, its evidence pages, the candidate database and one page per candidate.

    Order matters because Notion lists a child page where it was created:
    Run page (head blocks) -> evidence child pages -> remaining Run-page blocks
    -> candidate database. Raises NotionExportError if the Run page or the
    database cannot be created (the message then names the Run page that
    already exists). A failed evidence or candidate page is recorded in
    ``failed`` and skipped; after MAX_CONSECUTIVE_CANDIDATE_FAILURES candidate
    failures in a row the rest are skipped (``aborted``) rather than waiting
    out retries for every remaining page.
    """
    language = normalize_language(language)
    layout = layout_report(sections, language)
    timestamp = f"{generated_at or datetime.now():%Y-%m-%d %H:%M}"
    run_title = t("notion.run_title", language, title=layout.title, name=run_name, timestamp=timestamp)

    run_page = exporter.create_page(
        {"page_id": parent_page_id},
        {"title": {"title": rich_text(run_title)}},
        layout.run_head,
        "create Run page",
    )
    run_page_id = run_page["id"]
    run_page_url = str(run_page.get("url", ""))
    where = f" (the Run page was already created: {run_page_url or run_page_id})"

    failed: list[str] = []
    evidence_created = 0
    try:
        for evidence in layout.evidence_pages:
            try:
                exporter.create_page(
                    {"page_id": run_page_id},
                    {"title": {"title": rich_text(evidence.title)}},
                    evidence.blocks,
                    f"create evidence page '{evidence.title[:60]}'",
                )
                evidence_created += 1
            except NotionExportError as exc:
                failed.append(evidence.title)
                exporter.log(str(exc))
        exporter.append_blocks(run_page_id, layout.run_tail, "add Run page blocks")

        if not layout.candidates:
            return NotionExportResult(
                run_page_id, run_page_url, "", 0, 0, tuple(failed), False, len(layout.evidence_pages), evidence_created
            )

        database = exporter.create_database(
            run_page_id, t("notion.database_title", language), database_properties(language), "create candidate database"
        )
    except NotionExportError as exc:
        raise NotionExportError(f"{exc}{where}") from exc
    database_id = database["id"]

    created = 0
    consecutive_failures = 0
    aborted = False
    for index, page in enumerate(layout.candidates):
        if consecutive_failures >= MAX_CONSECUTIVE_CANDIDATE_FAILURES:
            aborted = True
            failed.extend(remaining.title for remaining in layout.candidates[index:])
            exporter.log(
                f"Notion export stopped after {MAX_CONSECUTIVE_CANDIDATE_FAILURES} consecutive failures; "
                f"{len(layout.candidates) - index} candidate page(s) were not created."
            )
            break
        if not page.blocks:
            exporter.log(f"Candidate '{page.title[:60]}' has no narrative text to export.")
        try:
            exporter.create_page(
                {"database_id": database_id},
                candidate_properties(page, language),
                page.blocks,
                f"create candidate page '{page.title[:60]}'",
            )
            created += 1
            consecutive_failures = 0
        except NotionExportError as exc:
            failed.append(page.title)
            consecutive_failures += 1
            exporter.log(str(exc))

    return NotionExportResult(
        run_page_id,
        run_page_url,
        database_id,
        len(layout.candidates),
        created,
        tuple(failed),
        aborted,
        len(layout.evidence_pages),
        evidence_created,
    )


def create_client(token: str, http_client: Any | None = None) -> Any:
    """A ``notion_client.Client`` on API version 2022-06-28 with the SDK's own retries and chatter off.

    ``http_client`` (an ``httpx.Client``) lets tests inspect the real HTTP requests.
    """
    from notion_client import Client

    options: dict[str, Any] = {"auth": token, "notion_version": NOTION_VERSION, "log_level": logging.ERROR}
    if http_client is not None:
        options["client"] = http_client
    try:
        return Client(retry=False, **options)
    except TypeError:  # older SDKs have no built-in retry option to switch off
        return Client(**options)


def export_run_to_notion(
    config: Any,
    sections: Sequence[NarrativeSection],
    *,
    token: str,
    run_name: str,
    generated_at: datetime | None = None,
    client: Any | None = None,
    log: Callable[[str], None] | None = None,
    **exporter_options: Any,
) -> NotionExportResult:
    """Export ``sections`` to the page configured in ``config.notion_export``; the token is passed in, never stored."""
    settings = config.notion_export
    exporter = NotionExporter(client if client is not None else create_client(token), log=log, **exporter_options)
    return export_sections_to_notion(
        sections,
        exporter=exporter,
        parent_page_id=settings.parent_page_id,
        language=getattr(config, "report_language", DEFAULT_LANGUAGE),
        run_name=run_name,
        generated_at=generated_at,
    )


__all__: tuple[str, ...] = (
    "NOTION_VERSION",
    "NotionExportResult",
    "NotionExporter",
    "batch_blocks",
    "candidate_properties",
    "classify_error",
    "create_client",
    "database_properties",
    "export_run_to_notion",
    "export_sections_to_notion",
    "layout_report",
    "rich_text",
    "section_to_blocks",
)
