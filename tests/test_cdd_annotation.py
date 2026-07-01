"""Tests for CDD annotation helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import pytest
import requests

from annotation.cdd import (
    domain_hit_from_dict,
    domain_hit_to_dict,
    parse_cdd_response,
    search_cdd_by_sequence,
)
from core.cache import JsonCache
from core.exceptions import CDDAnnotationError
from core.models import DomainHit


def test_parse_cdd_response_with_simple_domain_like_response() -> None:
    """A simple CDD-like line should parse into a DomainHit."""
    text = (
        "# CDD results\n"
        "query1\tcd12345\tThioredoxin_like\tredox domain\t1e-20\t55.5\t10-80\n"
    )

    hits = parse_cdd_response(text)

    assert hits == [
        DomainHit(
            source="CDD",
            accession="cd12345",
            name="Thioredoxin_like",
            description="redox domain",
            evalue=1e-20,
            bitscore=55.5,
            start=10,
            end=80,
        )
    ]


def test_parse_cdd_response_with_labeled_values() -> None:
    """Labeled CDD-like fields should also parse."""
    text = "query1 cd54321 DomainName description evalue=2e-10 bitscore=44 start=5 end=50"

    hits = parse_cdd_response(text)

    assert len(hits) == 1
    assert hits[0].accession == "cd54321"
    assert hits[0].evalue == 2e-10
    assert hits[0].bitscore == 44.0
    assert hits[0].start == 5
    assert hits[0].end == 50


def test_parse_cdd_response_with_no_hits_returns_empty_list() -> None:
    """Text without domain accessions should return no hits."""
    assert parse_cdd_response("# no hits\nquery complete\n") == []


def test_search_cdd_by_sequence_returns_cached_hits_without_request(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Cached CDD hits should be returned without calling requests."""
    cache = JsonCache(tmp_path)
    hit = DomainHit(
        source="CDD",
        accession="cd12345",
        name="Thioredoxin_like",
        description="redox domain",
        evalue=1e-20,
        bitscore=55.5,
        start=10,
        end=80,
    )
    cache.set("cdd", "protein_1", [domain_hit_to_dict(hit)])
    post_mock = Mock()
    monkeypatch.setattr("annotation.cdd.requests.post", post_mock)

    hits = search_cdd_by_sequence("protein_1", "MSTNPKPQR", cache=cache)

    assert hits == [hit]
    post_mock.assert_not_called()


def test_search_cdd_by_sequence_empty_sequence_returns_empty_list(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Blank sequences should not be sent to CDD."""
    post_mock = Mock()
    monkeypatch.setattr("annotation.cdd.requests.post", post_mock)

    assert search_cdd_by_sequence("protein_1", "   ") == []
    post_mock.assert_not_called()


def test_search_cdd_by_sequence_request_failure_raises_cdd_annotation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Request failures should be wrapped in CDDAnnotationError."""
    monkeypatch.setattr(
        "annotation.cdd.requests.post",
        Mock(side_effect=requests.RequestException("network down")),
    )

    with pytest.raises(CDDAnnotationError, match="CDD search failed"):
        search_cdd_by_sequence("protein_1", "MSTNPKPQR")


def test_search_cdd_by_sequence_parses_and_caches_response(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """CDD response text should be parsed and stored in cache."""
    cache = JsonCache(tmp_path)
    response = Mock()
    response.raise_for_status.return_value = None
    response.text = "query1\tcd12345\tDomain\tDescription\t1e-5\t30.0\t1-20\n"
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.cdd.requests.post", post_mock)

    hits = search_cdd_by_sequence("protein_1", "MSTNPKPQR", cache=cache, timeout=5)

    assert hits[0].accession == "cd12345"
    assert cache.has("cdd", "protein_1") is True
    assert post_mock.call_args.kwargs["timeout"] == 5


def test_domain_hit_to_dict_and_from_dict_round_trip() -> None:
    """DomainHit cache serialization should round trip cleanly."""
    hit = DomainHit(
        source="CDD",
        accession="cd12345",
        name="Thioredoxin_like",
        description="redox domain",
        evalue=1e-20,
        bitscore=55.5,
        start=10,
        end=80,
    )

    data = domain_hit_to_dict(hit)
    restored = domain_hit_from_dict(data)

    assert restored == hit
