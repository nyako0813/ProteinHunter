"""Tests for analysis/domain_family_map.py."""

from __future__ import annotations

from pathlib import Path

import pytest

from annotation.uniprot_bulk import UniProtDomainInfo
from core.exceptions import ConfigError
from analysis.domain_family_map import (
    DomainFamilyCategory,
    DomainFamilyMap,
    load_domain_family_map,
)

SAMPLE_YAML = """
version: test
categories:
  activator:
    pfam: [PF24167]
    interpro: [IPR055834]
    supfam: []
  carrier:
    pfam: [PF12724]
    interpro: []
    supfam: [SSF52218]
category_rules:
  - left: activator
    right: carrier
    note: "activator x carrier"
"""


def test_default_domain_family_map_loads(tmp_path: Path) -> None:
    path = tmp_path / "domain_family_map.yaml"
    path.write_text(SAMPLE_YAML, encoding="utf-8")

    loaded = load_domain_family_map(path)

    assert loaded.version == "test"
    assert loaded.categories["activator"] == DomainFamilyCategory(
        pfam=frozenset({"PF24167"}), interpro=frozenset({"IPR055834"}), supfam=frozenset()
    )
    assert loaded.category_rules == (("activator", "carrier", "activator x carrier"),)


def test_categories_for_matches_on_pfam_interpro_or_supfam(tmp_path: Path) -> None:
    path = tmp_path / "domain_family_map.yaml"
    path.write_text(SAMPLE_YAML, encoding="utf-8")
    loaded = load_domain_family_map(path)

    activator_info = UniProtDomainInfo(pfam=frozenset({"PF24167"}))
    carrier_info = UniProtDomainInfo(supfam=frozenset({"SSF52218"}))
    unrelated_info = UniProtDomainInfo(pfam=frozenset({"PF00001"}))

    assert loaded.categories_for(activator_info) == {"activator"}
    assert loaded.categories_for(carrier_info) == {"carrier"}
    assert loaded.categories_for(unrelated_info) == set()
    assert loaded.categories_for(None) == set()


def test_find_category_match_detects_either_orientation(tmp_path: Path) -> None:
    path = tmp_path / "domain_family_map.yaml"
    path.write_text(SAMPLE_YAML, encoding="utf-8")
    loaded = load_domain_family_map(path)

    forward = loaded.find_category_match({"activator"}, {"carrier"})
    backward = loaded.find_category_match({"carrier"}, {"activator"})

    assert forward == ("activator", "carrier", "activator x carrier")
    assert backward == ("activator", "carrier", "activator x carrier")


def test_find_category_match_returns_none_when_no_rule_matches(tmp_path: Path) -> None:
    path = tmp_path / "domain_family_map.yaml"
    path.write_text(SAMPLE_YAML, encoding="utf-8")
    loaded = load_domain_family_map(path)

    assert loaded.find_category_match({"activator"}, {"activator"}) is None
    assert loaded.find_category_match(set(), set()) is None


def test_missing_file_raises_config_error(tmp_path: Path) -> None:
    with pytest.raises(ConfigError):
        load_domain_family_map(tmp_path / "missing.yaml")


def test_invalid_yaml_raises_config_error(tmp_path: Path) -> None:
    path = tmp_path / "broken.yaml"
    path.write_text("categories: [this is not, a mapping", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_domain_family_map(path)


def test_category_rule_referencing_unknown_category_rejected(tmp_path: Path) -> None:
    path = tmp_path / "domain_family_map.yaml"
    path.write_text(
        """
version: test
categories:
  activator:
    pfam: [PF24167]
category_rules:
  - left: activator
    right: does_not_exist
""",
        encoding="utf-8",
    )
    with pytest.raises(ConfigError):
        load_domain_family_map(path)


def test_shipped_v1_domain_family_map_loads() -> None:
    path = Path(__file__).resolve().parent.parent / "config" / "domain_family_map.v1.yaml"
    loaded = load_domain_family_map(path)

    assert loaded.version == "v1"
    assert loaded.categories["atp_dependent_activator"].pfam == frozenset({"PF24167"})
    assert loaded.categories["electron_carrier"].pfam == frozenset({"PF12724"})
    match = loaded.find_category_match({"atp_dependent_activator"}, {"electron_carrier"})
    assert match is not None
    assert match[0] == "atp_dependent_activator"
    assert match[1] == "electron_carrier"
