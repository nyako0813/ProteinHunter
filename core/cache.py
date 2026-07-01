"""JSON file cache helpers for ProteinHunter.

The cache stores each annotation namespace in its own JSON file. This keeps the
storage format simple, readable, and easy to remove when a cached result needs
to be rebuilt.
"""

from __future__ import annotations

import json
from json import JSONDecodeError
from pathlib import Path
from typing import Any

from core.exceptions import CacheError


JsonValue = dict[str, Any] | list[Any] | str | int | float | bool | None


class JsonCache:
    """A small namespace-based cache backed by UTF-8 JSON files."""

    def __init__(self, cache_dir: str | Path) -> None:
        """Create a cache that stores JSON files inside ``cache_dir``."""
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)

    def get(self, namespace: str, key: str) -> JsonValue:
        """Return a cached value, or None when the namespace or key is absent."""
        cache_data = self._read_namespace(namespace)
        return cache_data.get(key)

    def set(self, namespace: str, key: str, value: Any) -> None:
        """Save a value under ``namespace`` and ``key``."""
        cache_data = self._read_namespace(namespace)
        cache_data[key] = value
        self._write_namespace(namespace, cache_data)

    def has(self, namespace: str, key: str) -> bool:
        """Return True when ``key`` exists in ``namespace``."""
        cache_data = self._read_namespace(namespace)
        return key in cache_data

    def delete(self, namespace: str, key: str) -> bool:
        """Delete a cached key and return True when something was removed."""
        cache_data = self._read_namespace(namespace)

        if key not in cache_data:
            return False

        del cache_data[key]
        self._write_namespace(namespace, cache_data)
        return True

    def clear_namespace(self, namespace: str) -> None:
        """Remove all cached entries from one namespace."""
        self._write_namespace(namespace, {})

    def clear_all(self) -> None:
        """Remove all entries from every JSON cache file in the cache directory."""
        for cache_file in self.cache_dir.glob("*.json"):
            cache_file.write_text("{}\n", encoding="utf-8")

    def count(self, namespace: str) -> int:
        """Return the number of cached entries in one namespace."""
        return len(self._read_namespace(namespace))

    def _namespace_path(self, namespace: str) -> Path:
        """Return the JSON file path for a namespace."""
        return self.cache_dir / f"{namespace}.json"

    def _read_namespace(self, namespace: str) -> dict[str, JsonValue]:
        """Read a namespace JSON file and return its dictionary contents."""
        cache_file = self._namespace_path(namespace)

        if not cache_file.exists():
            return {}

        try:
            data = json.loads(cache_file.read_text(encoding="utf-8"))
        except JSONDecodeError as exc:
            message = (
                f"The cache file '{cache_file}' is not valid JSON. "
                "Please delete it or rebuild the cache."
            )
            raise CacheError(message) from exc

        if not isinstance(data, dict):
            message = (
                f"The cache file '{cache_file}' must contain a JSON object. "
                "Please delete it or rebuild the cache."
            )
            raise CacheError(message)

        return data

    def _write_namespace(self, namespace: str, data: dict[str, JsonValue]) -> None:
        """Write one namespace dictionary to disk as formatted JSON."""
        cache_file = self._namespace_path(namespace)
        cache_file.write_text(
            json.dumps(data, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )


__all__: tuple[str, ...] = (
    "JsonCache",
    "JsonValue",
)
