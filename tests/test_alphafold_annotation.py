"""Tests for AlphaFold annotation helpers."""

from __future__ import annotations

from unittest.mock import Mock

import pytest
import requests

from annotation.alphafold import (
    build_alphafold_url,
    check_alphafold_exists,
    get_alphafold_url_if_exists,
)
from core.cache import JsonCache
from core.exceptions import AlphaFoldAnnotationError


def test_build_alphafold_url_with_accession_and_none() -> None:
    """AlphaFold URLs should be built only for non-empty accessions."""
    assert build_alphafold_url("P12345") == "https://alphafold.ebi.ac.uk/entry/P12345"
    assert build_alphafold_url(None) is None
    assert build_alphafold_url("") is None


def test_alphafold_exists_true(monkeypatch: pytest.MonkeyPatch) -> None:
    """HTTP success should mean an AlphaFold entry exists."""
    response = Mock(status_code=200)
    head_mock = Mock(return_value=response)
    monkeypatch.setattr("annotation.alphafold.requests.head", head_mock)

    assert check_alphafold_exists("P12345", timeout=5) is True
    assert get_alphafold_url_if_exists("P12345") == (
        "https://alphafold.ebi.ac.uk/entry/P12345"
    )
    assert head_mock.call_args_list[0].kwargs["timeout"] == 5


def test_alphafold_404_returns_false(monkeypatch: pytest.MonkeyPatch) -> None:
    """A 404 response should mean no AlphaFold prediction is available."""
    monkeypatch.setattr(
        "annotation.alphafold.requests.head",
        Mock(return_value=Mock(status_code=404)),
    )

    assert check_alphafold_exists("P00000") is False
    assert get_alphafold_url_if_exists("P00000") is None


def test_alphafold_cache_hit_avoids_request(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Cached AlphaFold availability should avoid the network request."""
    cache = JsonCache(tmp_path)
    cache.set("alphafold", "P12345", True)
    head_mock = Mock()
    monkeypatch.setattr("annotation.alphafold.requests.head", head_mock)

    assert check_alphafold_exists("P12345", cache=cache) is True
    assert get_alphafold_url_if_exists("P12345", cache=cache) == (
        "https://alphafold.ebi.ac.uk/entry/P12345"
    )
    head_mock.assert_not_called()


def test_alphafold_request_error_raises_annotation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Request failures should be wrapped in AlphaFoldAnnotationError."""
    monkeypatch.setattr(
        "annotation.alphafold.requests.head",
        Mock(side_effect=requests.RequestException("network down")),
    )

    with pytest.raises(AlphaFoldAnnotationError, match="AlphaFold DB check failed"):
        check_alphafold_exists("P12345")


def test_alphafold_unexpected_status_raises_annotation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unexpected HTTP status codes should raise AlphaFoldAnnotationError."""
    monkeypatch.setattr(
        "annotation.alphafold.requests.head",
        Mock(return_value=Mock(status_code=500)),
    )

    with pytest.raises(AlphaFoldAnnotationError, match="status 500"):
        check_alphafold_exists("P12345")
