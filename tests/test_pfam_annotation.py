"""Tests for Pfam annotation helpers."""

from __future__ import annotations

import json
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
    assert "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan" in message
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
    assert "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan" in message


def test_search_pfam_by_sequence_uses_api_endpoint_and_json_body(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pfam searches should use the HMMER API endpoint and JSON body first."""
    response = Mock()
    response.status_code = 200
    response.text = "# no hits\n"
    response.raise_for_status.return_value = None
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    search_pfam_by_sequence("protein_1", "MSTNPKPQR", timeout=9)

    call = post_mock.call_args
    assert call.args[0] == "https://www.ebi.ac.uk/Tools/hmmer/api/v1/search/hmmscan"
    assert call.kwargs["json"] == {"input": "MSTNPKPQR", "database": "pfam"}
    assert call.kwargs["headers"] == {
        "Accept": "application/json",
        "Content-Type": "application/json",
    }
    assert call.kwargs["timeout"] == 9


def test_search_pfam_by_sequence_retries_form_body_when_json_body_is_rejected(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A JSON body parse error from HMMER should trigger one form retry."""
    json_response = Mock()
    json_response.status_code = 400
    json_response.text = '{"detail": "Cannot parse request body"}'
    form_response = Mock()
    form_response.status_code = 200
    form_response.text = "# no hits\n"
    form_response.raise_for_status.return_value = None
    post_mock = Mock(side_effect=[json_response, form_response])
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    assert search_pfam_by_sequence("protein_1", "MSTNPKPQR") == []

    first_call = post_mock.call_args_list[0]
    second_call = post_mock.call_args_list[1]
    assert first_call.kwargs["json"] == {"input": "MSTNPKPQR", "database": "pfam"}
    assert first_call.kwargs["headers"] == {
        "Accept": "application/json",
        "Content-Type": "application/json",
    }
    assert second_call.kwargs["data"] == {"input": "MSTNPKPQR", "database": "pfam"}
    assert second_call.kwargs["headers"] == {"Accept": "application/json"}


def test_search_pfam_by_sequence_parses_json_hit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A JSON HMMER API response with one Pfam hit should become a DomainHit."""
    response = Mock()
    response.status_code = 200
    response.text = json.dumps(
        dict(
            results=[
                dict(
                    accession="PF00001",
                    name="7tm_1",
                    description="Seven transmembrane receptor",
                    evalue="1e-20",
                    score=55.5,
                    start=10,
                    end=80,
                )
            ]
        )
    )
    response.raise_for_status.return_value = None
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    hits = search_pfam_by_sequence("protein_1", "MSTNPKPQR")

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


def test_search_pfam_by_sequence_fetches_result_when_initial_response_has_id(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An asynchronous job id should be resolved through the result endpoint."""
    submit_response = Mock()
    submit_response.status_code = 200
    submit_response.text = json.dumps(dict(id="job-123", status="RUNNING"))
    submit_response.raise_for_status.return_value = None
    result_response = Mock()
    result_response.status_code = 200
    result_response.text = json.dumps(
        dict(
            status="SUCCESS",
            results=[
                dict(
                    accession="PF00002",
                    name="ABC_tran",
                    description="ABC transporter",
                    i_evalue="3e-12",
                    bitscore=72.0,
                    ali_from=4,
                    ali_to=95,
                )
            ],
        )
    )
    result_response.raise_for_status.return_value = None
    post_mock = Mock(return_value=submit_response)
    get_mock = Mock(return_value=result_response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)
    monkeypatch.setattr("annotation.pfam.requests.get", get_mock)

    hits = search_pfam_by_sequence("protein_1", "MSTNPKPQR")

    assert get_mock.call_args.args[0] == (
        "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/job-123"
    )
    assert hits == [
        DomainHit(
            source="Pfam",
            accession="PF00002",
            name="ABC_tran",
            description="ABC transporter",
            evalue=3e-12,
            bitscore=72.0,
            start=4,
            end=95,
        )
    ]


def test_search_pfam_by_sequence_result_error_status_raises_clear_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """FAILURE or ERROR result statuses should be reported as result failures."""
    submit_response = Mock()
    submit_response.status_code = 200
    submit_response.text = json.dumps(dict(id="job-123"))
    submit_response.raise_for_status.return_value = None
    result_response = Mock()
    result_response.status_code = 200
    result_response.text = json.dumps(dict(status="ERROR", detail="Job failed"))
    result_response.raise_for_status.return_value = None
    monkeypatch.setattr("annotation.pfam.requests.post", Mock(return_value=submit_response))
    monkeypatch.setattr("annotation.pfam.requests.get", Mock(return_value=result_response))

    with pytest.raises(PfamAnnotationError) as exc_info:
        search_pfam_by_sequence("protein_1", "MSTNPKPQR")

    message = str(exc_info.value)
    assert "result phase" in message
    assert "status ERROR" in message
    assert "https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/job-123" in message
    assert "Response preview:" in message


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
    """Unexpected JSON responses should include available top-level keys."""
    response = Mock()
    response.status_code = 200
    response.text = json.dumps(dict(message="done", status="DONE"))
    response.raise_for_status.return_value = None
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.pfam.requests.post", post_mock)

    with pytest.raises(PfamAnnotationError) as exc_info:
        search_pfam_by_sequence("protein_1", "MSTNPKPQR")

    message = str(exc_info.value)
    assert "parse phase" in message
    assert "no Pfam domain fields were recognized" in message
    assert "Available top-level keys: message, status" in message
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
