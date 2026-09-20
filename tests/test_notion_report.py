"""Tests for the Notion export (output/notion_report.py) against a fake Notion client."""

from __future__ import annotations

from datetime import datetime
from types import SimpleNamespace
from typing import Any, Callable

import pytest

from core.exceptions import NotionExportError
from output import notion_report as nr
from output.report_sections import build_report_sections
from output.word_narrative import CategoryRef, NarrativeSection, bullet_list, heading, paragraph, table, toc

NOW = datetime(2026, 9, 20, 12, 30)
REFS = (
    CategoryRef("candidate_priority_score", "sequence", "Sequence/Source Classification", 30.0),
    CategoryRef("same_gene_neighborhood_score", "genomic_context", "Genomic Context", 25.0),
)


def row(query: str, candidate: str, rank: int, tier: str = "Tier2_Strong", source: str = "Candidates", **extra: Any) -> dict:
    base = {
        "query_id": query,
        "candidate_protein_id": candidate,
        "candidate_description": "a hypothetical enzyme",
        "candidate_rank": rank,
        "candidate_source": source,
        "negative_hit_strength": "none",
        "final_score": 61.25,
        "final_score_tier": tier,
        "evidence_category_count": 3,
        "candidate_priority_score": 30.0,
        "same_gene_neighborhood_score": 12.5,
    }
    base.update(extra)
    return base


GROUPED = [
    ("MA_4115", [row("MA_4115", "MA_0363", 1), row("MA_4115", "MA_4110", 2, tier="Tier3_Moderate", source="No_hit", final_score=33.0)]),
    ("MA_0688", [row("MA_0688", "MA_0687", 1, tier="Tier1_VeryStrong", final_score=88.5)]),
]


def sections(language: str = "en", grouped=GROUPED) -> list[NarrativeSection]:
    return build_report_sections(
        language=language,
        generated_at=NOW,
        scoring_model="v2_evidence_based",
        category_caps={"source_classification": 30.0, "genomic_context": 25.0},
        pih_bundle_configured=False,
        grouped=grouped,
        category_refs=REFS,
        max_per_query=15,
        excel_filename="run.xlsx",
    )


# ---------------------------------------------------------------------------
# Fake Notion client
# ---------------------------------------------------------------------------


class FakeClient:
    """Records every call; ``fail(name, kwargs, n)`` may return an exception to raise instead."""

    def __init__(self, fail: Callable[[str, dict, int], Exception | None] | None = None) -> None:
        self.calls: list[tuple[str, dict]] = []
        self._fail = fail
        self._ids = 0
        self.pages = SimpleNamespace(create=lambda **kw: self._do("pages.create", kw, self._page))
        self.databases = SimpleNamespace(create=lambda **kw: self._do("databases.create", kw, self._database))
        self.blocks = SimpleNamespace(children=SimpleNamespace(append=lambda **kw: self._do("blocks.children.append", kw, lambda kw: {"results": []})))

    def _do(self, name: str, kwargs: dict, make: Callable[[dict], dict]) -> dict:
        self.calls.append((name, kwargs))
        error = self._fail(name, kwargs, len(self.calls)) if self._fail else None
        if error is not None:
            raise error
        return make(kwargs)

    def _next(self, prefix: str) -> str:
        self._ids += 1
        return f"{prefix}-{self._ids}"

    def _page(self, kwargs: dict) -> dict:
        page_id = self._next("page")
        return {"id": page_id, "url": f"https://www.notion.so/{page_id}"}

    def _database(self, kwargs: dict) -> dict:
        return {"id": self._next("db")}

    def named(self, name: str) -> list[dict]:
        return [kwargs for call, kwargs in self.calls if call == name]


class HttpError(Exception):
    """Duck-typed stand-in for notion_client.errors.APIResponseError."""

    def __init__(self, status: int, message: str = "boom", headers: dict | None = None, code: str = "err") -> None:
        super().__init__(message)
        self.status = status
        self.code = code
        self.headers = headers or {}


def make_exporter(client: FakeClient, **overrides: Any) -> tuple[nr.NotionExporter, list[float], list[str]]:
    sleeps: list[float] = []
    logs: list[str] = []
    options: dict[str, Any] = dict(sleep=sleeps.append, jitter=lambda delay: 0.0, min_interval=0.0, log=logs.append)
    options.update(overrides)
    return nr.NotionExporter(client, **options), sleeps, logs


def export(client: FakeClient, language: str = "en", grouped=GROUPED, **overrides: Any):
    exporter, sleeps, logs = make_exporter(client, **overrides)
    result = nr.export_sections_to_notion(
        sections(language, grouped),
        exporter=exporter,
        parent_page_id="3e1ef408bf298005a303f93efcff9871",
        language=language,
        run_name="MA_4115_2",
        generated_at=NOW,
    )
    return result, exporter, sleeps, logs


# ---------------------------------------------------------------------------
# Blocks
# ---------------------------------------------------------------------------


def test_rich_text_is_split_at_the_2000_character_limit() -> None:
    text = "あ" * 4500

    items = nr.rich_text(text)

    assert [len(item["text"]["content"]) for item in items] == [2000, 2000, 500]
    assert "".join(item["text"]["content"] for item in items) == text
    assert nr.rich_text("") == []


def test_rich_text_bold_annotation_only_when_requested() -> None:
    assert nr.rich_text("x", bold=True)[0]["annotations"] == {"bold": True}
    assert "annotations" not in nr.rich_text("x")[0]


def test_headings_map_to_notion_levels_and_the_title_becomes_no_block() -> None:
    assert nr.section_to_blocks(heading("Title", 0)) == []
    assert [nr.section_to_blocks(heading("h", level))[0]["type"] for level in (1, 2, 3, 4)] == ["heading_1", "heading_2", "heading_3", "heading_3"]


def test_labelled_paragraph_has_a_bold_lead_in_then_the_text() -> None:
    block = nr.section_to_blocks(paragraph("body", label="Why: "))[0]

    rich = block["paragraph"]["rich_text"]
    assert [item["text"]["content"] for item in rich] == ["Why: ", "body"]
    assert rich[0]["annotations"] == {"bold": True} and "annotations" not in rich[1]


def test_long_paragraph_stays_one_block_with_several_rich_text_objects() -> None:
    blocks = nr.section_to_blocks(paragraph("x" * 4100, label="L: "))

    assert len(blocks) == 1
    lengths = [len(item["text"]["content"]) for item in blocks[0]["paragraph"]["rich_text"]]
    assert lengths == [3, 2000, 2000, 100] and max(lengths) <= nr.RICH_TEXT_LIMIT


def test_bullet_list_and_toc_blocks() -> None:
    bullets = nr.section_to_blocks(bullet_list(["a", "b"]))
    assert [b["type"] for b in bullets] == ["bulleted_list_item"] * 2
    assert bullets[1]["bulleted_list_item"]["rich_text"][0]["text"]["content"] == "b"
    assert nr.section_to_blocks(toc("placeholder")) == [{"object": "block", "type": "table_of_contents", "table_of_contents": {}}]


def test_table_block_has_header_flag_width_and_nested_rows_with_empty_cells_allowed() -> None:
    block = nr.section_to_blocks(table(["Rank", "Candidate"], [("1", "MA_0363"), ("2", "")]))[0]

    assert block["type"] == "table"
    assert block["table"]["table_width"] == 2 and block["table"]["has_column_header"] is True
    rows = block["table"]["children"]
    assert len(rows) == 3 and all(r["type"] == "table_row" for r in rows)
    assert rows[1]["table_row"]["cells"][1][0]["text"]["content"] == "MA_0363"
    assert rows[2]["table_row"]["cells"][1] == []


def test_long_tables_are_split_with_the_header_repeated() -> None:
    blocks = nr.section_to_blocks(table(["h"], [(str(i),) for i in range(250)]))

    assert len(blocks) == 3
    assert [len(b["table"]["children"]) for b in blocks] == [100, 100, 53]  # 99 + 99 + 52 body rows, each with the header row
    assert all(b["table"]["children"][0]["table_row"]["cells"][0][0]["text"]["content"] == "h" for b in blocks)


def test_unknown_section_kind_is_rejected() -> None:
    with pytest.raises(ValueError, match="unknown report section kind"):
        nr.section_to_blocks(NarrativeSection("carousel"))


def test_batches_respect_the_100_block_limit() -> None:
    blocks = [nr._block("paragraph", nr.rich_text(str(i))) for i in range(250)]

    assert [len(batch) for batch in nr.batch_blocks(blocks)] == [100, 100, 50]


def test_batches_also_respect_the_total_block_cap_including_table_rows() -> None:
    tables = [nr.section_to_blocks(table(["h"], [(str(i),) for i in range(59)]))[0] for _ in range(20)]  # 1 + 60 rows each

    batches = list(nr.batch_blocks(tables))

    assert all(sum(nr.count_blocks(b) for b in batch) <= nr.MAX_BLOCKS_PER_REQUEST for batch in batches)
    assert sum(len(batch) for batch in batches) == 20 and len(batches) > 1


# ---------------------------------------------------------------------------
# Layout and properties
# ---------------------------------------------------------------------------


def test_layout_splits_run_page_from_candidate_pages() -> None:
    layout = nr.layout_report(sections("en"), "en")

    assert layout.title == "ProteinHunter Candidate Report"
    assert [c.title for c in layout.candidates] == [
        "MA_0363 — a hypothetical enzyme",
        "MA_4110 — a hypothetical enzyme",
        "MA_0687 — a hypothetical enzyme",
    ]
    assert [c.meta["query_id"] for c in layout.candidates] == ["MA_4115", "MA_4115", "MA_0688"]
    types = [b["type"] for b in layout.run_blocks]
    assert "table_of_contents" in types and types.count("table") == 2  # the per-query ranking tables stay on the Run page
    assert "heading_1" in types
    run_text = " ".join(item["text"]["content"] for b in layout.run_blocks for item in b.get(b["type"], {}).get("rich_text", []))
    assert "Candidate Details" in run_text and "One page per candidate" in run_text
    assert "Why this candidate ranks highly" not in run_text  # that text lives on the candidate pages


def test_candidate_page_body_is_the_three_narrative_paragraphs() -> None:
    page = nr.layout_report(sections("en"), "en").candidates[0]

    assert [b["type"] for b in page.blocks] == ["paragraph"] * 3
    first = page.blocks[0]["paragraph"]["rich_text"]
    assert first[0]["text"]["content"] == "Why this candidate ranks highly: " and first[0]["annotations"] == {"bold": True}
    assert "MA_0363" in page.blocks[2]["paragraph"]["rich_text"][0]["text"]["content"]  # the Excel cross-reference


def test_layout_without_candidates_keeps_everything_on_the_run_page() -> None:
    layout = nr.layout_report(sections("en", grouped=[]), "en")

    assert layout.candidates == []
    text = " ".join(item["text"]["content"] for b in layout.run_blocks for item in b.get(b["type"], {}).get("rich_text", []))
    assert "No query-specific candidates were produced by this run." in text
    assert "One page per candidate" not in text


def test_layout_is_localized_but_carries_the_same_data() -> None:
    en, ja = nr.layout_report(sections("en"), "en"), nr.layout_report(sections("ja"), "ja")

    assert [c.meta for c in en.candidates] == [c.meta for c in ja.candidates]
    assert ja.title == "ProteinHunter 候補レポート"
    ja_text = " ".join(item["text"]["content"] for b in ja.run_blocks for item in b.get(b["type"], {}).get("rich_text", []))
    assert "候補ごとのページは、下のデータベースにある" in ja_text


def test_database_schema_and_property_names_follow_the_language() -> None:
    en = nr.database_properties("en")
    ja = nr.database_properties("ja")

    assert list(en) == ["Candidate", "Query", "Rank", "Final Score", "Tier", "Candidate Source", "Candidate ID"]
    assert list(ja) == ["候補", "クエリ", "順位", "最終スコア", "Tier(信頼度階層)", "候補ソース", "候補ID"]
    assert en["Candidate"] == {"title": {}} and en["Rank"] == {"number": {"format": "number"}}
    assert en["Tier"] == {"select": {}} and en["Candidate ID"] == {"rich_text": {}}
    assert [list(v)[0] for v in en.values()].count("title") == 1  # Notion needs exactly one title property


def test_candidate_properties_carry_data_values_unchanged() -> None:
    page = nr.layout_report(sections("en"), "en").candidates[1]

    props = nr.candidate_properties(page, "en")

    assert props["Candidate"]["title"][0]["text"]["content"] == "MA_4110 — a hypothetical enzyme"
    assert props["Query"] == {"select": {"name": "MA_4115"}}
    assert props["Rank"] == {"number": 2} and props["Final Score"] == {"number": 33.0}
    assert props["Tier"] == {"select": {"name": "Tier3_Moderate"}}
    assert props["Candidate Source"] == {"select": {"name": "No_hit"}}
    assert props["Candidate ID"]["rich_text"][0]["text"]["content"] == "MA_4110"
    assert nr.candidate_properties(page, "ja")["クエリ"] == {"select": {"name": "MA_4115"}}  # same value, Japanese column


def test_select_values_drop_commas_and_are_capped_and_missing_numbers_are_omitted() -> None:
    page = nr.CandidatePage("t", {"query_id": "a,b" + "x" * 200, "rank": "", "final_score": "n/a", "tier": "", "source": "S"}, [])

    props = nr.candidate_properties(page, "en")

    name = props["Query"]["select"]["name"]
    assert "," not in name and len(name) <= nr.SELECT_NAME_LIMIT
    assert "Rank" not in props and "Final Score" not in props and "Tier" not in props
    assert props["Candidate Source"] == {"select": {"name": "S"}}


# ---------------------------------------------------------------------------
# Export flow
# ---------------------------------------------------------------------------


def test_export_builds_run_page_then_database_then_one_page_per_candidate() -> None:
    client = FakeClient()

    result, _, _, _ = export(client)

    names = [name for name, _ in client.calls]
    assert names == ["pages.create", "databases.create", "pages.create", "pages.create", "pages.create"]
    run_call, db_call, *candidate_calls = [kwargs for _, kwargs in client.calls]
    assert run_call["parent"] == {"page_id": "3e1ef408bf298005a303f93efcff9871"}
    assert run_call["properties"]["title"]["title"][0]["text"]["content"] == "ProteinHunter Candidate Report — MA_4115_2 (2026-09-20 12:30)"
    assert db_call["parent"] == {"type": "page_id", "page_id": "page-1"} and db_call["is_inline"] is True
    assert db_call["title"][0]["text"]["content"] == "Candidates"
    assert all(call["parent"] == {"database_id": "db-2"} for call in candidate_calls)
    assert result.run_page_url == "https://www.notion.so/page-1" and result.database_id == "db-2"
    assert (result.candidates_total, result.candidates_created, result.failed, result.aborted) == (3, 3, (), False)


def test_japanese_export_uses_japanese_titles_and_columns_with_unchanged_data() -> None:
    client = FakeClient()

    export(client, "ja")

    run_call, db_call, first = [kwargs for _, kwargs in client.calls][:3]
    assert run_call["properties"]["title"]["title"][0]["text"]["content"] == "ProteinHunter 候補レポート — MA_4115_2(2026-09-20 12:30)"
    assert db_call["title"][0]["text"]["content"] == "候補" and "クエリ" in db_call["properties"]
    assert first["properties"]["クエリ"] == {"select": {"name": "MA_4115"}}
    assert first["properties"]["候補ID"]["rich_text"][0]["text"]["content"] == "MA_0363"


def test_run_page_blocks_beyond_100_are_appended_in_further_requests() -> None:
    many = [(f"Q{i}", [row(f"Q{i}", f"C{i}", 1)]) for i in range(60)]  # 60 queries -> 60 headings + 60 tables on the Run page
    client = FakeClient()

    export(client, grouped=many)

    creates = client.named("pages.create")
    assert len(creates[0]["children"]) <= nr.MAX_CHILDREN_PER_REQUEST
    appends = client.named("blocks.children.append")
    assert appends and appends[0]["block_id"] == "page-1"
    for call in (creates[0], *appends):
        blocks = call["children"]
        assert len(blocks) <= nr.MAX_CHILDREN_PER_REQUEST
        assert sum(nr.count_blocks(b) for b in blocks) <= nr.MAX_BLOCKS_PER_REQUEST


def test_a_failing_candidate_page_is_recorded_and_the_rest_continue() -> None:
    def fail(name: str, kwargs: dict, n: int) -> Exception | None:
        if name == "pages.create" and kwargs["parent"].get("database_id") and "MA_4110" in kwargs["properties"]["Candidate"]["title"][0]["text"]["content"]:
            return HttpError(400, "validation_error")
        return None

    client = FakeClient(fail)

    result, _, _, logs = export(client)

    assert result.candidates_created == 2 and result.failed == ("MA_4110 — a hypothetical enzyme",) and not result.aborted
    assert any("MA_4110" in message and "status=400" in message for message in logs)


def test_export_stops_after_repeated_candidate_failures() -> None:
    def fail(name: str, kwargs: dict, n: int) -> Exception | None:
        return HttpError(400) if name == "pages.create" and "database_id" in kwargs["parent"] else None

    grouped = [("Q", [row("Q", f"C{i}", i + 1) for i in range(6)])]

    result, _, _, logs = export(FakeClient(fail), grouped=grouped)

    assert result.aborted and result.candidates_created == 0
    assert len(result.failed) == 6  # 3 tried, 3 skipped
    assert any("stopped after 3 consecutive failures" in message for message in logs)


def test_run_page_failure_is_fatal() -> None:
    client = FakeClient(lambda name, kwargs, n: HttpError(404, "object_not_found", code="object_not_found"))

    with pytest.raises(NotionExportError, match="create Run page.*status=404"):
        export(client)


def test_database_failure_is_fatal() -> None:
    client = FakeClient(lambda name, kwargs, n: HttpError(400) if name == "databases.create" else None)

    with pytest.raises(NotionExportError, match="create candidate database"):
        export(client)


def test_no_database_is_created_when_there_are_no_candidates() -> None:
    client = FakeClient()

    result, _, _, _ = export(client, grouped=[])

    assert [name for name, _ in client.calls] == ["pages.create"]
    assert (result.database_id, result.candidates_total) == ("", 0)


# ---------------------------------------------------------------------------
# Retries and throttling
# ---------------------------------------------------------------------------


def flaky(errors: list[Exception]) -> Callable[[str, dict, int], Exception | None]:
    remaining = list(errors)
    return lambda name, kwargs, n: remaining.pop(0) if remaining else None


def test_rate_limited_requests_honour_retry_after() -> None:
    client = FakeClient(flaky([HttpError(429, "rate_limited", {"retry-after": "3"})]))
    exporter, sleeps, logs = make_exporter(client)

    exporter.call("x", client.pages.create, parent={}, properties={})

    assert sleeps == [3.0] and len(client.calls) == 2
    assert "rate_limited" in logs[0] or "status=429" in logs[0]


def test_server_errors_back_off_exponentially_up_to_the_cap() -> None:
    client = FakeClient(flaky([HttpError(503)] * 5))
    exporter, sleeps, _ = make_exporter(client, base_delay=1.0, max_delay=4.0)

    exporter.call("x", client.pages.create, parent={}, properties={})

    assert sleeps == [1.0, 2.0, 4.0, 4.0, 4.0]


def test_gives_up_after_max_retries_with_a_clear_error() -> None:
    client = FakeClient(lambda name, kwargs, n: HttpError(503, "unavailable"))
    exporter, sleeps, _ = make_exporter(client, max_retries=2)

    with pytest.raises(NotionExportError, match=r"\(x\).*status=503"):
        exporter.call("x", client.pages.create, parent={}, properties={})

    assert len(client.calls) == 3 and len(sleeps) == 2  # first try + 2 retries


@pytest.mark.parametrize("status", [400, 401, 403, 404])
def test_client_errors_are_not_retried(status: int) -> None:
    client = FakeClient(lambda name, kwargs, n: HttpError(status))
    exporter, sleeps, _ = make_exporter(client)

    with pytest.raises(NotionExportError):
        exporter.call("x", client.pages.create, parent={}, properties={})

    assert len(client.calls) == 1 and sleeps == []


@pytest.mark.parametrize("error", [TimeoutError("t"), ConnectionError("c"), type("ReadTimeout", (Exception,), {})("r"), type("RequestTimeoutError", (Exception,), {})("q")])
def test_network_errors_are_retried(error: Exception) -> None:
    client = FakeClient(flaky([error]))
    exporter, sleeps, _ = make_exporter(client)

    exporter.call("x", client.pages.create, parent={}, properties={})

    assert len(client.calls) == 2 and len(sleeps) == 1


def test_unexpected_exceptions_are_wrapped_not_retried() -> None:
    client = FakeClient(lambda name, kwargs, n: ValueError("bad payload"))
    exporter, sleeps, _ = make_exporter(client)

    with pytest.raises(NotionExportError, match="ValueError"):
        exporter.call("x", client.pages.create, parent={}, properties={})

    assert sleeps == []


def test_requests_are_throttled_to_the_minimum_interval() -> None:
    clock = [0.0]
    sleeps: list[float] = []

    def sleep(seconds: float) -> None:
        sleeps.append(seconds)
        clock[0] += seconds

    client = FakeClient()
    exporter = nr.NotionExporter(client, sleep=sleep, monotonic=lambda: clock[0], min_interval=0.34, jitter=lambda d: 0.0)

    for _ in range(3):
        exporter.call("x", client.pages.create, parent={}, properties={})

    assert sleeps == pytest.approx([0.34, 0.34])  # none before the first call, then one gap each


def test_default_jitter_stays_within_a_quarter_of_the_delay() -> None:
    exporter = nr.NotionExporter(FakeClient(), base_delay=4.0)

    delays = [exporter.backoff_delay(0) for _ in range(50)]

    assert all(4.0 <= d <= 5.0 for d in delays)


# ---------------------------------------------------------------------------
# Real SDK objects
# ---------------------------------------------------------------------------


def test_classification_of_the_real_sdk_errors() -> None:
    errors = pytest.importorskip("notion_client.errors")
    httpx = pytest.importorskip("httpx")

    def api_error(status: int, headers: dict[str, str] | None = None):
        return errors.APIResponseError("rate_limited", status, "msg", httpx.Headers(headers or {}), "{}")

    assert nr.classify_error(api_error(429, {"retry-after": "3"})) == (True, 3.0)
    assert nr.classify_error(api_error(429)) == (True, None)
    assert nr.classify_error(api_error(503)) == (True, None)
    assert nr.classify_error(api_error(401)) == (False, None)
    assert nr.classify_error(errors.RequestTimeoutError()) == (True, None)


def test_create_client_pins_the_api_version_and_turns_sdk_retries_off() -> None:
    pytest.importorskip("notion_client")

    client = nr.create_client("secret-token")

    assert client.options.notion_version == nr.NOTION_VERSION == "2022-06-28"
    assert client.options.retry is False
