"""Tests for CDD annotation helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import pytest
import requests

from annotation.cdd import (
    CDD_MAX_QUERIES_PER_BATCH,
    CDD_POLL_INTERVAL_SECONDS,
    domain_hit_from_dict,
    domain_hit_to_dict,
    fetch_cdd_batch_results,
    parse_cdd_batch_response,
    parse_cdd_response,
    poll_cdd_batch,
    search_cdd_by_sequence,
    submit_cdd_batch,
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


def test_submit_cdd_batch_returns_search_id(monkeypatch: pytest.MonkeyPatch) -> None:
    """A successful submission should parse and return the #cdsid value."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.text = (
        "#Batch CD-search tool\tNIH/NLM/NCBI\n"
        "#cdsid\tQM3-qcdsearch-ABCDEF\n"
        "#datatype\thitsConcise Results\n"
        "#status\t3\tmsg\tJob is still running\n"
    )
    post_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.cdd.requests.post", post_mock)

    cdsid = submit_cdd_batch([("protein_1", "MSTNPKPQR"), ("protein_2", "MKVSVLF")])

    assert cdsid == "QM3-qcdsearch-ABCDEF"
    sent_fasta = post_mock.call_args.kwargs["data"]["queries"]
    assert sent_fasta == ">protein_1\nMSTNPKPQR\n>protein_2\nMKVSVLF\n"


def test_submit_cdd_batch_rejects_empty_queries() -> None:
    """Submitting zero queries should fail fast instead of contacting NCBI."""
    with pytest.raises(CDDAnnotationError, match="at least one query"):
        submit_cdd_batch([])


def test_submit_cdd_batch_rejects_over_limit() -> None:
    """More than NCBI's documented 1,000-sequence limit should be rejected."""
    queries = [(f"protein_{i}", "MSTNPKPQR") for i in range(CDD_MAX_QUERIES_PER_BATCH + 1)]

    with pytest.raises(CDDAnnotationError, match="exceeds the 1000-sequence limit"):
        submit_cdd_batch(queries)


def test_submit_cdd_batch_missing_cdsid_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """A response without a #cdsid line should raise, not return a bogus ID."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.text = "#status\t2\n"
    monkeypatch.setattr("annotation.cdd.requests.post", Mock(return_value=response))

    with pytest.raises(CDDAnnotationError, match="did not return a search ID"):
        submit_cdd_batch([("protein_1", "MSTNPKPQR")])


def test_submit_cdd_batch_request_failure_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """Network failures during submission should raise CDDAnnotationError."""
    monkeypatch.setattr(
        "annotation.cdd.requests.post",
        Mock(side_effect=requests.RequestException("network down")),
    )

    with pytest.raises(CDDAnnotationError, match="submission failed"):
        submit_cdd_batch([("protein_1", "MSTNPKPQR")])


def test_poll_cdd_batch_waits_through_running_then_succeeds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Status 3 (still running) must be polled again; status 0 must stop polling.

    This is the exact bug being fixed: the old code treated the first
    ("still running") response as the final result.
    """
    running = Mock()
    running.raise_for_status.return_value = None
    running.text = "#status\t3\tmsg\tJob is still running\n"
    done = Mock()
    done.raise_for_status.return_value = None
    done.text = "#status\t0\n#status\tsuccess\n"
    post_mock = Mock(side_effect=[running, done])
    monkeypatch.setattr("annotation.cdd.requests.post", post_mock)
    sleep_calls: list[float] = []

    poll_cdd_batch("QM3-qcdsearch-ABCDEF", sleep_fn=sleep_calls.append)

    assert post_mock.call_count == 2
    assert sleep_calls == [CDD_POLL_INTERVAL_SECONDS, CDD_POLL_INTERVAL_SECONDS]


@pytest.mark.parametrize(
    ("status", "expected_message"),
    [
        ("1", "Invalid search ID"),
        ("2", "No effective input"),
        ("4", "Queue manager"),
        ("5", "Data is corrupted"),
    ],
)
def test_poll_cdd_batch_raises_on_terminal_error_status(
    monkeypatch: pytest.MonkeyPatch, status: str, expected_message: str
) -> None:
    """Every documented non-running, non-success status code must raise."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.text = f"#status\t{status}\n"
    monkeypatch.setattr("annotation.cdd.requests.post", Mock(return_value=response))

    with pytest.raises(CDDAnnotationError, match=expected_message):
        poll_cdd_batch("QM3-qcdsearch-ABCDEF", sleep_fn=lambda _seconds: None)


def test_poll_cdd_batch_times_out_while_still_running(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A job stuck at status 3 forever must eventually raise, not hang forever."""
    running = Mock()
    running.raise_for_status.return_value = None
    running.text = "#status\t3\tmsg\tJob is still running\n"
    monkeypatch.setattr("annotation.cdd.requests.post", Mock(return_value=running))

    with pytest.raises(CDDAnnotationError, match="did not complete within"):
        poll_cdd_batch(
            "QM3-qcdsearch-ABCDEF",
            poll_interval=1.0,
            max_wait=2.0,
            sleep_fn=lambda _seconds: None,
        )


def test_poll_cdd_batch_request_failure_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """Network failures during a status check should raise CDDAnnotationError."""
    monkeypatch.setattr(
        "annotation.cdd.requests.post",
        Mock(side_effect=requests.RequestException("network down")),
    )

    with pytest.raises(CDDAnnotationError, match="status check failed"):
        poll_cdd_batch("QM3-qcdsearch-ABCDEF", sleep_fn=lambda _seconds: None)


def test_poll_cdd_batch_malformed_response_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """A response without a recognizable #status line should raise, not loop forever."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.text = "not a batch cd-search response at all"
    monkeypatch.setattr("annotation.cdd.requests.post", Mock(return_value=response))

    with pytest.raises(CDDAnnotationError, match="unexpected response"):
        poll_cdd_batch("QM3-qcdsearch-ABCDEF", sleep_fn=lambda _seconds: None)


def test_fetch_cdd_batch_results_returns_response_text(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The raw completed-job response text should be returned as-is."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.text = "#status\t0\n\nQuery\tHit type\t...\n"
    monkeypatch.setattr("annotation.cdd.requests.post", Mock(return_value=response))

    text = fetch_cdd_batch_results("QM3-qcdsearch-ABCDEF")

    assert text == response.text


def test_fetch_cdd_batch_results_request_failure_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Network failures during result retrieval should raise CDDAnnotationError."""
    monkeypatch.setattr(
        "annotation.cdd.requests.post",
        Mock(side_effect=requests.RequestException("network down")),
    )

    with pytest.raises(CDDAnnotationError, match="result retrieval failed"):
        fetch_cdd_batch_results("QM3-qcdsearch-ABCDEF")


def test_parse_cdd_batch_response_assigns_hits_to_correct_query() -> None:
    """Multi-query Concise Results must route each hit to its own protein_id.

    This exact text (aside from re-wrapping) was captured from a live,
    completed 2-sequence Batch CD-Search job during the CDD investigation --
    it is real NCBI output, not a guessed format.
    """
    text = (
        "#Batch CD-search tool\tNIH/NLM/NCBI\n"
        "#cdsid\tQM3-qcdsearch-384A85A1D6118CB1-934BD050CEA111D\n"
        "#datatype\thitsConcise Results\n"
        "#status\t0\n"
        "#Start time\t2026-09-01T14:59:40\tRun time\t0:00:00:04\n"
        "#status\tsuccess\n"
        "\n"
        "Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\t"
        "Short name\tIncomplete\tSuperfamily\n"
        "Q#1 - >WP_011024006.1\tspecific\t441720\t2\t198\t8.87211e-118\t332.168\t"
        "COG2117\tCOG2117\t - \tcl42515\n"
        "Q#2 - >WP_011023824.1\tspecific\t441771\t9\t78\t1.01218e-14\t62.9264\t"
        "COG2168\tTusB\tC\tcl46469\n"
    )

    hits_by_query = parse_cdd_batch_response(text)

    assert set(hits_by_query) == {"WP_011024006.1", "WP_011023824.1"}
    hit1 = hits_by_query["WP_011024006.1"][0]
    assert hit1.accession == "COG2117"
    assert hit1.name == "COG2117"
    assert hit1.evalue == pytest.approx(8.87211e-118)
    assert hit1.bitscore == pytest.approx(332.168)
    assert hit1.start == 2
    assert hit1.end == 198
    hit2 = hits_by_query["WP_011023824.1"][0]
    assert hit2.accession == "COG2168"
    assert hit2.name == "TusB"
    assert hit2.evalue == pytest.approx(1.01218e-14)
    assert hit2.start == 9
    assert hit2.end == 78


def test_parse_cdd_batch_response_no_hits_returns_empty_dict() -> None:
    """A completed job with no domain hits should return an empty mapping, not raise."""
    text = (
        "#status\t0\n"
        "Query\tHit type\tPSSM-ID\tFrom\tTo\tE-Value\tBitscore\tAccession\t"
        "Short name\tIncomplete\tSuperfamily\n"
    )

    assert parse_cdd_batch_response(text) == {}


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
