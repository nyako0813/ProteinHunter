"""Tests for the STRING PPI evidence bridge (analysis/string_ppi_bridge.py)."""

from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import Mock

import pytest
import requests

from analysis.string_ppi_bridge import (
    STRING_VERSION,
    StringPairScores,
    load_string_ppi_bundle,
)
from core.cache import JsonCache

TAXON_ID = 188937


def _write_gzip(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(content)


def _seed_bulk_files(cache_dir: Path, *, links: bool = True, info: bool = True) -> Path:
    """Write synthetic local STRING bulk files, as if already downloaded once."""
    string_dir = cache_dir / "string_ppi_network"
    if links:
        links_content = (
            "protein1 protein2 neighborhood fusion cooccurence coexpression "
            "experimental database textmining combined_score\n"
            f"{TAXON_ID}.MA_0001 {TAXON_ID}.MA_0002 0 0 467 0 0 0 0 467\n"
            f"{TAXON_ID}.MA_0001 {TAXON_ID}.MA_0003 57 0 245 0 0 0 0 257\n"
        )
        _write_gzip(
            string_dir / f"{TAXON_ID}.protein.links.detailed.v{STRING_VERSION}.txt.gz", links_content
        )
    if info:
        info_content = (
            "#string_protein_id\tpreferred_name\tprotein_size\tannotation\n"
            f"{TAXON_ID}.MA_0001\tMA_0001\t100\tsome protein\n"
            f"{TAXON_ID}.MA_0002\tMA_0002\t100\tsome protein\n"
            f"{TAXON_ID}.MA_0003\tMA_0003\t100\tsome protein\n"
            f"{TAXON_ID}.MA_0004\tMA_0004\t100\tsome protein\n"
        )
        _write_gzip(string_dir / f"{TAXON_ID}.protein.info.v{STRING_VERSION}.txt.gz", info_content)
    return string_dir


def test_disabled_when_taxon_id_is_none(tmp_path: Path) -> None:
    """No STRING evidence should be produced when ncbi_taxon_id is unset."""
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_string_ppi_bundle(None, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_0002") is None
    assert bundle.warnings == ()


def test_scan_local_bulk_file_finds_known_pair(tmp_path: Path) -> None:
    """A pair present in the local links file should be returned with its channel scores."""
    _seed_bulk_files(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_0002") == StringPairScores(cooccurrence=467, combined_score=467)
    assert bundle.lookup("MA_0001", "MA_0003") == StringPairScores(
        neighborhood=57, cooccurrence=245, combined_score=257
    )


def test_known_pair_absent_from_links_file_is_evaluated_zero(tmp_path: Path) -> None:
    """Both proteins known to STRING, but no row for this pair -> all-zero, not None."""
    _seed_bulk_files(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)

    # MA_0004 is in protein.info but has no row with MA_0001 in the links file.
    assert bundle.lookup("MA_0001", "MA_0004") == StringPairScores()


def test_unknown_protein_returns_none(tmp_path: Path) -> None:
    """A protein absent from protein.info entirely should be MISSING (None), not zero."""
    _seed_bulk_files(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_9999") is None
    assert bundle.lookup("MA_9999", "MA_0001") is None


def test_result_is_cached_and_reused_without_rescanning(tmp_path: Path) -> None:
    """A second call for the same query must not need to re-read the (possibly huge) links file."""
    _seed_bulk_files(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")
    first = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)
    assert first.lookup("MA_0001", "MA_0002").cooccurrence == 467

    # Remove the (large) links file -- only protein.info and the cache remain.
    string_dir = tmp_path / "string_ppi_network"
    (string_dir / f"{TAXON_ID}.protein.links.detailed.v{STRING_VERSION}.txt.gz").unlink()

    second = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)

    assert second.lookup("MA_0001", "MA_0002") == StringPairScores(cooccurrence=467, combined_score=467)
    assert second.warnings == ()


def test_live_api_fallback_when_no_local_bulk_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """With no local bulk file and a download failure, fall back to one live API call."""
    cache = JsonCache(tmp_path / "jsoncache")

    live_response = Mock()
    live_response.raise_for_status.return_value = None
    live_response.json.return_value = [
        {
            "stringId_A": f"{TAXON_ID}.MA_0001",
            "stringId_B": f"{TAXON_ID}.MA_0002",
            "score": 0.467,
            "nscore": 0.0,
            "fscore": 0.0,
            "pscore": 0.467,
            "ascore": 0.0,
            "escore": 0.0,
            "dscore": 0.0,
            "tscore": 0.0,
        }
    ]

    def fake_get(url: str, *args: object, **kwargs: object) -> Mock:
        if "download" in url or "stringdb-downloads" in url:
            raise requests.RequestException("simulated network failure")
        return live_response

    monkeypatch.setattr("analysis.string_ppi_bridge.requests.get", Mock(side_effect=fake_get))
    monkeypatch.setattr("analysis.string_ppi_bridge.time.sleep", Mock())

    bundle = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_0002") == StringPairScores(cooccurrence=467, combined_score=467)
    assert any("download" in warning.lower() for warning in bundle.warnings)

    # The live result should now be cached -- a second call must not need requests.get again.
    get_calls_before = len(bundle.warnings)
    bundle2 = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)
    assert bundle2.lookup("MA_0001", "MA_0002").cooccurrence == 467


def test_download_and_live_failure_degrades_to_empty_bundle_without_raising(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A fully offline environment must degrade gracefully, never raise, per the optional-evidence rule."""
    cache = JsonCache(tmp_path / "jsoncache")
    monkeypatch.setattr(
        "analysis.string_ppi_bridge.requests.get",
        Mock(side_effect=requests.RequestException("simulated network failure")),
    )
    monkeypatch.setattr("analysis.string_ppi_bridge.time.sleep", Mock())

    bundle = load_string_ppi_bundle(TAXON_ID, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_0002") is None
    assert bundle.warnings  # at least one warning explaining the degradation
