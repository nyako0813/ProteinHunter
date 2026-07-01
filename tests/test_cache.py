"""Tests for the JSON cache helper."""

from __future__ import annotations

from pathlib import Path

import pytest

from core.cache import JsonCache
from core.exceptions import CacheError


def test_set_get_has_delete_and_count(tmp_path: Path) -> None:
    """Basic cache operations should work for one namespace."""
    cache = JsonCache(tmp_path)

    assert cache.get("pfam", "protein_1") is None
    assert cache.has("pfam", "protein_1") is False
    assert cache.count("pfam") == 0

    cache.set("pfam", "protein_1", {"domain": "PF00001", "score": 12.5})

    assert cache.has("pfam", "protein_1") is True
    assert cache.get("pfam", "protein_1") == {"domain": "PF00001", "score": 12.5}
    assert cache.count("pfam") == 1
    assert cache.delete("pfam", "protein_1") is True
    assert cache.delete("pfam", "protein_1") is False
    assert cache.get("pfam", "protein_1") is None
    assert cache.count("pfam") == 0


def test_namespaces_are_separated(tmp_path: Path) -> None:
    """Entries with the same key should stay separate by namespace."""
    cache = JsonCache(tmp_path)

    cache.set("cdd", "protein_1", {"source": "cdd"})
    cache.set("pfam", "protein_1", {"source": "pfam"})

    assert cache.get("cdd", "protein_1") == {"source": "cdd"}
    assert cache.get("pfam", "protein_1") == {"source": "pfam"}
    assert (tmp_path / "cdd.json").exists()
    assert (tmp_path / "pfam.json").exists()


def test_clear_namespace(tmp_path: Path) -> None:
    """Clearing one namespace should not affect another namespace."""
    cache = JsonCache(tmp_path)

    cache.set("uniprot", "protein_1", "P12345")
    cache.set("alphafold", "protein_1", "https://example.test/model")

    cache.clear_namespace("uniprot")

    assert cache.count("uniprot") == 0
    assert cache.get("uniprot", "protein_1") is None
    assert cache.get("alphafold", "protein_1") == "https://example.test/model"


def test_clear_all(tmp_path: Path) -> None:
    """Clearing all should empty every JSON cache file."""
    cache = JsonCache(tmp_path)

    cache.set("cdd", "protein_1", {"ok": True})
    cache.set("pfam", "protein_2", ["PF00001"])

    cache.clear_all()

    assert cache.count("cdd") == 0
    assert cache.count("pfam") == 0
    assert cache.get("cdd", "protein_1") is None
    assert cache.get("pfam", "protein_2") is None


def test_corrupted_json_raises_cache_error(tmp_path: Path) -> None:
    """Invalid JSON should raise CacheError instead of being ignored."""
    cache = JsonCache(tmp_path)
    (tmp_path / "cdd.json").write_text("{not valid json", encoding="utf-8")

    with pytest.raises(CacheError, match="not valid JSON"):
        cache.get("cdd", "protein_1")


def test_json_file_must_contain_object(tmp_path: Path) -> None:
    """A namespace cache file should contain a JSON object."""
    cache = JsonCache(tmp_path)
    (tmp_path / "pfam.json").write_text("[1, 2, 3]", encoding="utf-8")

    with pytest.raises(CacheError, match="must contain a JSON object"):
        cache.count("pfam")
