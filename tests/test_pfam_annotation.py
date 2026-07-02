"""Tests for Pfam annotation helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import pytest
import requests

from annotation.pfam import (
    domain_hit_from_dict,
    domain_hit_to_dict,
    parse_pfam_response,
    search_pfam_by_sequence,
)
from core.cache import JsonCache
from core.exceptions import PfamAnnotationError
from core.models import DomainHit


def test_parse_pfam_response_with_simple_pfam_like_response() -> None:
    """A simple Pfam-like line should parse into a DomainHit."""
    text = (
        "# Pfam results\n"
        "PF00001\t7tm_1\tSeven transmembrane receptor\t1e-20\t55.5\t10-80\n"
    )

    hits = parse_pfam_response(text)

    assert hits == [
        DomainHit(
            source="Pfam",
            accession="PF00001",
            name="7tm_1",
            description="Seven transmembrane receptor",
            evalue=1e-20,
            bitscore=55.5,
            start=10,
            end=80,
        )
    ]


def test_parse_pfam_response_with_labeled_values() -> None:
    """Labeled Pfam-like fields should also parse."""
    text = "DomainName PF12345 description evalue=2e-10 bitscore=44 start=5 end=50"

    hits = parse_pfam_response(text)

    assert len(hits) == 1
    assert hits[0].accession == "PF12345"
    assert hits[0].name == "DomainName"
    assert hits[0].evalue == 2e-10
    assert hits[0].bitscore == 44.0
    assert hits[0].start == 5
    assert hits[0].end == 50


def test_parse_pfam_response_with_no_hits_returns_empty_list() -> None:
    """Text without Pfam accessions should return no hits."""
    assert parse_pfam_response("# no hits\nquery complete\n") == []


def test_search_pfam_by_sequence_returns_cached_hits_without_request(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Cached Pfam hits should be returned without calling requests."""
    cache = JsonCache(tmp_path)
    hit = DomainHit(
        source="Pfam",
        accession="PF00001",
        name="7tm_1",
        description="Seven transmembrane receptor",
        evalue=1e-20,
        bitscore=55.5,
        start=10,
        end=80,
    )
    cache.set("pfam", "protein_1", [domain_hit_to_dict(hit)])
    post_mock = Mock()
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    hits = search_pfam_by_sequence("protein_1", "MSTNPKPQR", cache=cache)

    assert hits == [hit]
    post_mock.assert_not_called()


def test_search_pfam_by_sequence_empty_sequence_returns_empty_list(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Blank sequences should not be sent to Pfam."""
    post_mock = Mock()
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    assert search_pfam_by_sequence("protein_1", "   ") == []
    post_mock.assert_not_called()


def test_search_pfam_by_sequence_request_failure_raises_pfam_annotation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Request failures should be wrapped in PfamAnnotationError."""
    monkeypatch.setattr(
        "annotation.pfam.requests.post",
        Mock(side_effect=requests.RequestException("network down")),
    )

    with pytest.raises(PfamAnnotationError, match="Pfam search failed"):
        search_pfam_by_sequence("protein_1", "MSTNPKPQR")


def test_search_pfam_by_sequence_http_error_includes_status_and_preview(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """HTTP errors should include the status code and short response text."""
    response = Mock()
    response.status_code = 405
    response.text = "Method Not Allowed: please use a supported Pfam endpoint."
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    with pytest.raises(PfamAnnotationError) as exc_info:
        search_pfam_by_sequence("protein_1", "MSTNPKPQR", timeout=12)

    message = str(exc_info.value)
    assert "request phase" in message
    assert "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan" in message
    assert "HTTP status: 405" in message
    assert "Method Not Allowed" in message
    assert "Timeout setting: 12 seconds" in message


def test_search_pfam_by_sequence_timeout_includes_timeout_information(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Timeout errors should say that the request timed out."""
    monkeypatch.setattr(
        "annotation.pfam.requests.post",
        Mock(side_effect=requests.Timeout("timed out")),
    )

    with pytest.raises(PfamAnnotationError) as exc_info:
        search_pfam_by_sequence("protein_1", "MSTNPKPQR", timeout=7)

    message = str(exc_info.value)
    assert "request phase" in message
    assert "timed out after 7 seconds" in message
    assert "Timeout setting: 7 seconds" in message
    assert "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan" in message


def test_search_pfam_by_sequence_invalid_json_gives_parse_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Invalid JSON-like responses should produce a clear parse diagnostic."""
    response = Mock()
    response.status_code = 200
    response.text = "{not valid json"
    response.raise_for_status.return_value = None
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    with pytest.raises(PfamAnnotationError) as exc_info:
        search_pfam_by_sequence("protein_1", "MSTNPKPQR")

    message = str(exc_info.value)
    assert "parse phase" in message
    assert "invalid JSON" in message
    assert "Response preview: {not valid json" in message


def test_search_pfam_by_sequence_unexpected_json_format_gives_parse_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unexpected JSON responses should explain that text output was expected."""
    response = Mock()
    response.status_code = 200
    response.text = '{"results": []}'
    response.raise_for_status.return_value = None
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    with pytest.raises(PfamAnnotationError) as exc_info:
        search_pfam_by_sequence("protein_1", "MSTNPKPQR")

    message = str(exc_info.value)
    assert "parse phase" in message
    assert "expects text or tabular output" in message
    assert "Response preview:" in message


def test_search_pfam_by_sequence_parses_and_caches_response(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pfam response text should be parsed and stored in cache."""
    cache = JsonCache(tmp_path)
    response = Mock()
    response.raise_for_status.return_value = None
    response.text = "PF00001\tDomain\tDescription\t1e-5\t30.0\t1-20\n"
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    hits = search_pfam_by_sequence("protein_1", "MSTNPKPQR", cache=cache, timeout=5)

    assert hits[0].accession == "PF00001"
    assert cache.has("pfam", "protein_1") is True
    assert post_mock.call_args.kwargs["timeout"] == 5


def test_domain_hit_to_dict_and_from_dict_round_trip() -> None:
    """DomainHit cache serialization should round trip cleanly."""
    hit = DomainHit(
        source="Pfam",
        accession="PF00001",
        name="7tm_1",
        description="Seven transmembrane receptor",
        evalue=1e-20,
        bitscore=55.5,
        start=10,
        end=80,
    )

    data = domain_hit_to_dict(hit)
    restored = domain_hit_from_dict(data)

    assert restored == hit
