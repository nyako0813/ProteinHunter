"""Notion export: the report as Notion pages, a third output beside Excel and Word.

Renders the same renderer-neutral ``NarrativeSection`` list as the Word report
(``output/report_sections.py``), so wording, language (``report_language``) and
candidate selection are decided once. Nothing here changes any computed result.

Layout in the workspace (claude/notion_export_design.md):

    parent page (notion_export.parent_page_id, shared with the integration)
      └─ Run page          one per pipeline run: title-page lines, table of
                           contents, "5. Evidence Architecture", the per-query
                           ranking tables
           └─ Candidates   an inline database, one row per candidate
                └─ page    the candidate's "why it ranks highly" /
                           "biological interpretation" text, opened from the row

Database properties carry the data behind each row (query, rank, final score,
tier, candidate source, candidate id) so the database can be filtered and
grouped in Notion; their values are the same data values the Excel workbook
has and are never translated.

The Notion API is called through ``notion-client`` pinned to API version
2022-06-28, with the SDK's own automatic retries switched off: this module
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
class ReportLayout:
    title: str
    run_blocks: list[dict[str, Any]]
    candidates: list[CandidatePage]


def layout_report(sections: Sequence[NarrativeSection], language: str = DEFAULT_LANGUAGE) -> ReportLayout:
    """Split the report into the Run page's blocks and one page per candidate.

    Everything before the ``details`` heading, and that heading itself, is the
    Run page; each ``candidate`` heading opens a candidate page that collects
    the paragraphs after it. Per-query headings inside the details are dropped
    (the query is a database property instead). Without a ``details`` heading
    the whole report stays on the Run page.
    """
    language = normalize_language(language)
    title = next((s.text for s in sections if s.kind == "heading" and s.level == 0), "ProteinHunter")
    run_blocks: list[dict[str, Any]] = []
    candidates: list[CandidatePage] = []
    in_details = False
    current: tuple[str, dict[str, str], list[dict[str, Any]]] | None = None

    def close_current() -> None:
        nonlocal current
        if current is not None:
            candidates.append(CandidatePage(title=current[0], meta=current[1], blocks=current[2]))
            current = None

    details_note_index: int | None = None
    for section in sections:
        if section.kind == "heading" and section.role == "details":
            in_details = True
            run_blocks.extend(section_to_blocks(section))
            details_note_index = len(run_blocks)
            continue
        if in_details and section.kind == "heading" and section.role == "query":
            close_current()
            continue
        if in_details and section.kind == "heading" and section.role == "candidate":
            close_current()
            current = (section.text, section.meta_dict(), [])
            continue
        blocks = section_to_blocks(section)
        if in_details and current is not None:
            current[2].extend(blocks)
        else:
            run_blocks.extend(blocks)
    close_current()

    if candidates and details_note_index is not None:
        run_blocks.insert(details_note_index, _block("paragraph", rich_text(t("notion.details_note", language))))
    return ReportLayout(title=title, run_blocks=run_blocks, candidates=candidates)


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

    # -- page/database helpers, all sending blocks in limit-respecting batches ----

    def create_page(self, parent: dict[str, Any], properties: dict[str, Any], blocks: Sequence[dict[str, Any]], description: str) -> dict[str, Any]:
        batches = list(batch_blocks(blocks))
        first = batches[0] if batches else []
        page = self.call(
            description,
            self.client.pages.create,
            parent=parent,
            properties=properties,
            **({"children": first} if first else {}),
        )
        for batch in batches[1:]:
            self.call(f"{description} (more blocks)", self.client.blocks.children.append, block_id=page["id"], children=batch)
        return page


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


def export_sections_to_notion(
    sections: Sequence[NarrativeSection],
    *,
    exporter: NotionExporter,
    parent_page_id: str,
    language: str = DEFAULT_LANGUAGE,
    run_name: str = "",
    generated_at: datetime | None = None,
) -> NotionExportResult:
    """Create the Run page, the candidate database and one page per candidate under ``parent_page_id``.

    Raises NotionExportError if the Run page or the database cannot be created
    (nothing useful exists yet). Once the database exists, a candidate whose
    page cannot be created is recorded in ``failed`` and skipped; after
    MAX_CONSECUTIVE_CANDIDATE_FAILURES in a row the rest are skipped
    (``aborted``) rather than waiting out retries for every remaining page.
    """
    language = normalize_language(language)
    layout = layout_report(sections, language)
    timestamp = f"{generated_at or datetime.now():%Y-%m-%d %H:%M}"
    run_title = t("notion.run_title", language, title=layout.title, name=run_name, timestamp=timestamp)

    run_page = exporter.create_page(
        {"page_id": parent_page_id},
        {"title": {"title": rich_text(run_title)}},
        layout.run_blocks,
        "create Run page",
    )
    run_page_id = run_page["id"]
    run_page_url = str(run_page.get("url", ""))

    if not layout.candidates:
        return NotionExportResult(run_page_id, run_page_url, "", 0, 0)

    database = exporter.call(
        "create candidate database",
        exporter.client.databases.create,
        parent={"type": "page_id", "page_id": run_page_id},
        title=rich_text(t("notion.database_title", language)),
        is_inline=True,
        properties=database_properties(language),
    )
    database_id = database["id"]

    created = 0
    failed: list[str] = []
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

    return NotionExportResult(run_page_id, run_page_url, database_id, len(layout.candidates), created, tuple(failed), aborted)


def create_client(token: str) -> Any:
    """A ``notion_client.Client`` on API version 2022-06-28 with the SDK's own retries and chatter off."""
    from notion_client import Client

    try:
        return Client(auth=token, notion_version=NOTION_VERSION, retry=False, log_level=logging.ERROR)
    except TypeError:  # older SDKs have no built-in retry option to switch off
        return Client(auth=token, notion_version=NOTION_VERSION, log_level=logging.ERROR)


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
