"""Tests for the Rockhopper operon evidence bridge's runtime lookup half
(analysis/rockhopper_operon_bridge.py::RockhopperOperonBundle /
load_rockhopper_operon_bundle). The generation half (FASTQ acquisition,
Rockhopper invocation, GFF-based old_locus_tag resolution) is exercised
manually, not by these tests -- it requires network access, Java, and
Rockhopper itself, matching how the STRING/coexpression bridges' own
live-download paths are left untested here too.
"""

from __future__ import annotations

import json
from pathlib import Path

from analysis.rockhopper_operon_bridge import (
    RockhopperOperonHit,
    load_rockhopper_operon_bundle,
)


def _write_cache(path: Path, entries: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(entries))


def test_disabled_returns_empty_bundle_without_reading_cache(tmp_path: Path) -> None:
    """enabled=False must short-circuit before ever touching the cache path."""
    cache_path = tmp_path / "rockhopper_operons.json"
    # Deliberately do not create cache_path -- if load_rockhopper_operon_bundle
    # tried to read it while disabled, this would raise instead of returning
    # an empty bundle.
    bundle = load_rockhopper_operon_bundle(False, cache_path)

    assert bundle.lookup("MA_0001", "MA_0002") is None
    assert bundle.n_samples == 0
    assert bundle.warnings == ()


def test_missing_cache_file_degrades_to_empty_bundle_with_warning(tmp_path: Path) -> None:
    cache_path = tmp_path / "rockhopper_operons.json"

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    assert bundle.lookup("MA_0001", "MA_0002") is None
    assert bundle.n_samples == 0
    assert len(bundle.warnings) == 1
    assert "not found" in bundle.warnings[0]


def test_malformed_cache_file_degrades_to_empty_bundle_with_warning(tmp_path: Path) -> None:
    cache_path = tmp_path / "rockhopper_operons.json"
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text("{not valid json")

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    assert bundle.lookup("MA_0001", "MA_0002") is None
    assert bundle.n_samples == 0
    assert len(bundle.warnings) == 1


def test_pair_grouped_in_one_sample_is_a_hit(tmp_path: Path) -> None:
    cache_path = tmp_path / "rockhopper_operons.json"
    _write_cache(
        cache_path,
        [
            {
                "sample": "LK57",
                "condition": "acetate",
                "operon_groups": [["MA_3896", "MA_3897", "MA_3898", "MA_3899"]],
            }
        ],
    )

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    assert bundle.n_samples == 1
    assert bundle.lookup("MA_3896", "MA_3899") == RockhopperOperonHit(sample="LK57", condition="acetate")
    # Order-independent.
    assert bundle.lookup("MA_3899", "MA_3896") == RockhopperOperonHit(sample="LK57", condition="acetate")


def test_pair_known_but_never_co_grouped_is_missing_not_zero(tmp_path: Path) -> None:
    """The core asymmetric-design test: unlike STRING's "evaluated, zero"
    case for a known-but-absent pair, a Rockhopper non-merge must resolve
    to None (MISSING) here even though both genes clearly appear in the
    cache (just never together) -- see RockhopperOperonBundle.lookup's
    docstring and claude/phase6e_rockhopper_lk57_validation.md's Mtp
    false-negative finding.
    """
    cache_path = tmp_path / "rockhopper_operons.json"
    _write_cache(
        cache_path,
        [
            {
                "sample": "LK57",
                "condition": "acetate",
                "operon_groups": [["MA_4164", "MA_9999"], ["MA_4165", "MA_8888"]],
            }
        ],
    )

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    # Both MA_4164 and MA_4165 appear in the cache, just never together.
    assert bundle.lookup("MA_4164", "MA_4165") is None


def test_protein_absent_from_every_sample_is_missing(tmp_path: Path) -> None:
    cache_path = tmp_path / "rockhopper_operons.json"
    _write_cache(
        cache_path,
        [{"sample": "LK57", "condition": "acetate", "operon_groups": [["MA_0001", "MA_0002"]]}],
    )

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    assert bundle.lookup("MA_9999", "MA_0001") is None


def test_or_aggregation_across_samples(tmp_path: Path) -> None:
    """A pair grouped in a later sample but not the first is still a hit --
    OR-aggregation across samples/conditions (Phase 6f spec section 3).
    """
    cache_path = tmp_path / "rockhopper_operons.json"
    _write_cache(
        cache_path,
        [
            {"sample": "LK57", "condition": "acetate", "operon_groups": [["MA_0001", "MA_0002"]]},
            {"sample": "LK21", "condition": "methanol", "operon_groups": [["MA_0003", "MA_0004"]]},
        ],
    )

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    assert bundle.n_samples == 2
    assert bundle.lookup("MA_0003", "MA_0004") == RockhopperOperonHit(sample="LK21", condition="methanol")
    # Not grouped in either sample.
    assert bundle.lookup("MA_0001", "MA_0004") is None


def test_multi_gene_group_yields_all_pairwise_combinations(tmp_path: Path) -> None:
    cache_path = tmp_path / "rockhopper_operons.json"
    _write_cache(
        cache_path,
        [
            {
                "sample": "LK57",
                "condition": "acetate",
                "operon_groups": [["MA_4546", "MA_4547", "MA_4548", "MA_4549", "MA_4550"]],
            }
        ],
    )

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    assert bundle.lookup("MA_4546", "MA_4550") is not None
    assert bundle.lookup("MA_4547", "MA_4549") is not None


def test_self_pair_and_empty_tags_never_hit(tmp_path: Path) -> None:
    cache_path = tmp_path / "rockhopper_operons.json"
    _write_cache(
        cache_path,
        [{"sample": "LK57", "condition": "acetate", "operon_groups": [["MA_0001", "MA_0002"]]}],
    )

    bundle = load_rockhopper_operon_bundle(True, cache_path)

    assert bundle.lookup("MA_0001", "MA_0001") is None
    assert bundle.lookup("", "MA_0002") is None
    assert bundle.lookup("MA_0001", "") is None
