"""Tests for main._run_notion_export: opt-in, token from the environment, never fails the run."""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from main import _run_notion_export
from output import notion_report, report_sections

PARENT = "3e1ef408bf298005a303f93efcff9871"
TOKEN = "secret_do-not-log-me"


class Logger:
    def __init__(self) -> None:
        self.infos: list[str] = []
        self.warnings: list[str] = []

    def info(self, message: str) -> None:
        self.infos.append(message)

    def warning(self, message: str) -> None:
        self.warnings.append(message)

    @contextmanager
    def timer(self, name: str):
        yield

    def everything(self) -> str:
        return "\n".join(self.infos + self.warnings)


def config(enabled: bool) -> Any:
    return SimpleNamespace(notion_export=SimpleNamespace(enabled=enabled, parent_page_id=PARENT), report_language="en")


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> dict[str, Any]:
    recorded: dict[str, Any] = {"sections": 0, "exports": []}

    def fake_sections(*args: Any, **kwargs: Any) -> list:
        recorded["sections"] += 1
        recorded["excel_filename"] = kwargs.get("excel_filename")
        return ["section"]

    def fake_export(cfg: Any, sections: Any, **kwargs: Any) -> Any:
        recorded["exports"].append((cfg, sections, kwargs))
        return notion_report.NotionExportResult("run-id", "https://www.notion.so/run-id", "db-id", 3, 3)

    monkeypatch.setattr(report_sections, "build_run_sections", fake_sections)
    monkeypatch.setattr(notion_report, "export_run_to_notion", fake_export)
    return recorded


def run(logger: Logger, cfg: Any) -> None:
    _run_notion_export(logger, cfg, "classification", "interaction", Path("/out/MA_4115_2.xlsx"), "provenance")


def test_disabled_by_default_does_nothing(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.setenv("NOTION_TOKEN", TOKEN)
    logger = Logger()

    run(logger, config(False))

    assert calls["exports"] == [] and calls["sections"] == 0
    assert any("disabled" in message for message in logger.infos) and logger.warnings == []


def test_enabled_without_a_token_warns_and_skips(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.delenv("NOTION_TOKEN", raising=False)
    logger = Logger()

    run(logger, config(True))

    assert calls["exports"] == [] and calls["sections"] == 0
    assert len(logger.warnings) == 1 and "NOTION_TOKEN" in logger.warnings[0] and "unaffected" in logger.warnings[0]


def test_a_blank_token_counts_as_missing(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.setenv("NOTION_TOKEN", "   ")
    logger = Logger()

    run(logger, config(True))

    assert calls["exports"] == [] and logger.warnings


def test_enabled_with_a_token_exports_and_logs_the_run_page(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.setenv("NOTION_TOKEN", TOKEN)
    logger = Logger()

    run(logger, config(True))

    (cfg, sections, kwargs), = calls["exports"]
    assert sections == ["section"] and kwargs["token"] == TOKEN and kwargs["run_name"] == "MA_4115_2"
    assert calls["excel_filename"] == "MA_4115_2.xlsx"
    assert any("https://www.notion.so/run-id" in message for message in logger.infos)
    assert any("3 of 3" in message for message in logger.infos)
    assert logger.warnings == []


def test_the_token_is_never_logged(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.setenv("NOTION_TOKEN", TOKEN)
    logger = Logger()

    run(logger, config(True))

    assert TOKEN not in logger.everything()


def test_a_notion_failure_is_a_warning_not_an_error(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.setenv("NOTION_TOKEN", TOKEN)

    def boom(*args: Any, **kwargs: Any) -> Any:
        raise notion_report.NotionExportError("Notion request failed (create Run page): APIResponseError status=404")

    monkeypatch.setattr(notion_report, "export_run_to_notion", boom)
    logger = Logger()

    run(logger, config(True))  # must not raise

    assert len(logger.warnings) == 1 and "status=404" in logger.warnings[0] and "unaffected" in logger.warnings[0]
    assert TOKEN not in logger.everything()


def test_a_missing_notion_client_package_is_a_warning(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.setenv("NOTION_TOKEN", TOKEN)

    def no_package(*args: Any, **kwargs: Any) -> Any:
        raise ModuleNotFoundError("No module named 'notion_client'")

    monkeypatch.setattr(notion_report, "export_run_to_notion", no_package)
    logger = Logger()

    run(logger, config(True))

    assert len(logger.warnings) == 1 and "notion-client" in logger.warnings[0] and "requirements.txt" in logger.warnings[0]


def test_failed_candidate_pages_are_reported(monkeypatch: pytest.MonkeyPatch, calls: dict) -> None:
    monkeypatch.setenv("NOTION_TOKEN", TOKEN)
    monkeypatch.setattr(
        notion_report,
        "export_run_to_notion",
        lambda *a, **k: notion_report.NotionExportResult("r", "u", "d", 5, 2, ("A", "B", "C"), aborted=True),
    )
    logger = Logger()

    run(logger, config(True))

    assert any("2 of 5" in message for message in logger.infos)
    assert any("3 Notion candidate page(s)" in message and "stopped early" in message for message in logger.warnings)
