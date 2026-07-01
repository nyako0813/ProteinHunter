"""Tests for UniProt annotation helpers."""

from __future__ import annotations

from unittest.mock import Mock

import pytest
import requests

from annotation.uniprot import (
    extract_uniprot_accession,
    search_uniprot_by_protein_id,
)
from core.cache import JsonCache
from core.exceptions import UniProtAnnotationError


def test_uniprot_successful_response_parsing(monkeypatch: pytest.MonkeyPatch) -> None:
    """A UniProt result should be parsed into compact metadata."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.json.return_value = {
        "results": [
            {
                "primaryAccession": "P12345",
                "uniProtkbId": "TEST_PROTEIN",
                "entryType": "UniProtKB reviewed (Swiss-Prot)",
                "proteinDescription": {
                    "recommendedName": {"fullName": {"value": "Test protein"}}
                },
                "organism": {"scientificName": "Test organism"},
            }
        ]
    }
    get_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.uniprot.requests.get", get_mock)

    metadata = search_uniprot_by_protein_id("protein_1", timeout=5)

    assert metadata == {
        "query": "protein_1",
        "accession": "P12345",
        "id": "TEST_PROTEIN",
        "protein_name": "Test protein",
        "organism": "Test organism",
        "reviewed": True,
    }
    assert extract_uniprot_accession(metadata) == "P12345"
    get_mock.assert_called_once()
    assert get_mock.call_args.kwargs["timeout"] == 5


def test_uniprot_no_results_returns_accession_none(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No UniProt result should return metadata with accession None."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.json.return_value = {"results": []}
    monkeypatch.setattr("annotation.uniprot.requests.get", Mock(return_value=response))

    metadata = search_uniprot_by_protein_id("missing")

    assert metadata["query"] == "missing"
    assert metadata["accession"] is None
    assert extract_uniprot_accession(metadata) is None


def test_uniprot_cache_hit_avoids_request(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Cached UniProt metadata should be returned without a request."""
    cache = JsonCache(tmp_path)
    cache.set(
        "uniprot",
        "protein_1",
        {
            "query": "protein_1",
            "accession": "Q99999",
            "id": "CACHED",
            "protein_name": "Cached protein",
            "organism": "Cached organism",
            "reviewed": False,
        },
    )
    get_mock = Mock()
    monkeypatch.setattr("annotation.uniprot.requests.get", get_mock)

    metadata = search_uniprot_by_protein_id("protein_1", cache=cache)

    assert metadata["accession"] == "Q99999"
    assert metadata["id"] == "CACHED"
    get_mock.assert_not_called()


def test_uniprot_request_error_raises_annotation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Request failures should be wrapped in UniProtAnnotationError."""
    monkeypatch.setattr(
        "annotation.uniprot.requests.get",
        Mock(side_effect=requests.RequestException("network down")),
    )

    with pytest.raises(UniProtAnnotationError, match="UniProt search failed"):
        search_uniprot_by_protein_id("protein_1")


def test_uniprot_invalid_response_raises_annotation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unexpected UniProt JSON shapes should raise UniProtAnnotationError."""
    response = Mock()
    response.raise_for_status.return_value = None
    response.json.return_value = {"results": "not a list"}
    monkeypatch.setattr("annotation.uniprot.requests.get", Mock(return_value=response))

    with pytest.raises(UniProtAnnotationError, match="unexpected response"):
        search_uniprot_by_protein_id("protein_1")
