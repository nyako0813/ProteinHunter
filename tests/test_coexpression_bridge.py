"""Tests for the GEO coexpression evidence bridge (analysis/coexpression_bridge.py)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import requests

from analysis.coexpression_bridge import (
    GSE77738_STEADY_STATE_RPKM_COLUMNS,
    load_gse64349_coexpression_bundle,
    load_gse77738_coexpression_bundle,
)
from core.cache import JsonCache

# A handful of real steady-state column names, enough samples to get a
# well-defined (non-degenerate) correlation without needing all 13.
_SAMPLE_COLUMNS = sorted(GSE77738_STEADY_STATE_RPKM_COLUMNS)[:6]


def _seed_gse77738(cache_dir: Path) -> Path:
    """Write a small synthetic GSE77738_ReadCounts.xls-equivalent workbook."""
    path = cache_dir / "coexpression" / "GSE77738_ReadCounts.xls"
    path.parent.mkdir(parents=True, exist_ok=True)

    gene_loci = ["MA0001", "MA0002", "MA0003", "MA0004", "MAt4684"]
    gene_names = ["cdc6", "-", "-", "-", "-"]
    # MA0001/MA0002 move together (perfectly correlated); MA0003 moves
    # opposite; MA0004 is flat (zero variance); MAt4684 is a tRNA feature
    # that must never resolve to an old_locus_tag.
    base = [10, 20, 30, 40, 50, 60]
    values = {
        "MA0001": base,
        "MA0002": [v * 2 for v in base],
        "MA0003": [max(base) - v for v in base],
        "MA0004": [25, 25, 25, 25, 25, 25],
        "MAt4684": [1, 2, 3, 4, 5, 6],
    }
    data = {"Gene Locus": gene_loci, "Gene Name": gene_names}
    for col_index, col in enumerate(_SAMPLE_COLUMNS):
        data[col] = [values[gene][col_index] for gene in gene_loci]

    df = pd.DataFrame(data)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="Raw Read Counts", index=False)
        df.to_excel(writer, sheet_name="RPKM Normalized Read Counts", index=False)
    return path


def _seed_gse64349(cache_dir: Path) -> tuple[Path, Path]:
    """Write small synthetic TableS1 (wild-type) and TableS2 (parental+mutant) workbooks."""
    coexpr_dir = cache_dir / "coexpression"
    coexpr_dir.mkdir(parents=True, exist_ok=True)
    table1_path = coexpr_dir / "GSE64349_TableS1_GEO.xlsx"
    table2_path = coexpr_dir / "GSE64349_TableS2_GEO.xlsx"

    # TableS1: Feature ID mixes a gene symbol ('cdc6_1', resolved via
    # GSE77738's own Gene Name column) with bare locus ids.
    table1 = pd.DataFrame(
        {
            "Feature ID": ["cdc6_1", "MA0002", "MA0003"],
            "DMS - S1_R1 (single) (GE) - RPKM": [12, 24, 40],
            "DMS - S2_R1 (single) (GE) - RPKM": [14, 28, 38],
            "MeOH - S3_R1 (single) (GE) - RPKM": [10, 20, 42],
        }
    )
    table1.to_excel(table1_path, index=False)

    # TableS2: WWM82 (parental, kept) + delta-msrH (mutant, must be excluded).
    table2 = pd.DataFrame(
        {
            "Feature ID": ["cdc6_1", "MA0002", "MA0003"],
            "WWM82 (parental strain) - S4_R1 (single) (GE) - RPKM": [11, 22, 41],
            "delta-msrH - S5_R1 (single) (GE) - RPKM": [999, 999, 999],
        }
    )
    table2.to_excel(table2_path, index=False)
    return table1_path, table2_path


def test_disabled_returns_empty_bundle(tmp_path: Path) -> None:
    """No coexpression evidence should be produced when the bridge is disabled."""
    cache = JsonCache(tmp_path / "jsoncache")

    bundle77738 = load_gse77738_coexpression_bundle(False, ["MA_0001"], cache, tmp_path)
    bundle64349 = load_gse64349_coexpression_bundle(False, ["MA_0001"], cache, tmp_path)

    assert bundle77738.lookup("MA_0001", "MA_0002") is None
    assert bundle64349.lookup("MA_0001", "MA_0002") is None
    assert bundle77738.warnings == ()
    assert bundle64349.warnings == ()


def test_gse77738_correlated_genes_rank_above_anticorrelated(tmp_path: Path) -> None:
    """MA_0002 (moves in lockstep with MA_0001) should score far above MA_0003 (moves opposite)."""
    _seed_gse77738(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_gse77738_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)

    positive = bundle.lookup("MA_0001", "MA_0002")
    negative = bundle.lookup("MA_0001", "MA_0003")
    assert positive is not None and negative is not None
    assert positive.correlation == pytest.approx(1.0, abs=1e-3)
    assert negative.correlation < 0
    assert positive.percentile > negative.percentile
    assert 0.0 <= positive.percentile <= 1.0
    assert 0.0 <= negative.percentile <= 1.0


def test_gse77738_trna_feature_and_unknown_gene_are_missing(tmp_path: Path) -> None:
    """A tRNA feature (MAt####) never resolves to an old_locus_tag; an absent gene is MISSING."""
    _seed_gse77738(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_gse77738_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_4684") is None  # MAt4684 never becomes MA_4684
    assert bundle.lookup("MA_0001", "MA_9999") is None
    assert bundle.lookup("MA_9999", "MA_0001") is None


def test_gse77738_zero_variance_gene_is_unavailable_not_an_error(tmp_path: Path) -> None:
    """A flat-expression gene (MA0004) cannot be correlated; must not crash, and must not
    produce a false-positive lookup for that specific pair."""
    _seed_gse77738(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_gse77738_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_0004") is None
    assert any("zero-variance" in w for w in bundle.warnings)


def test_gse64349_excludes_mutant_includes_parental_strain(tmp_path: Path) -> None:
    """Delta-msrH samples must never influence the result; WWM82 (parental) samples must be pooled in."""
    _seed_gse77738(tmp_path)  # needed to build the symbol->locus table
    _seed_gse64349(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_gse64349_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)

    # 3 TableS1 wild-type samples + 1 TableS2 WWM82 (parental) sample = 4;
    # the delta-msrH column must not be counted.
    assert bundle.n_samples == 4
    result = bundle.lookup("MA_0001", "MA_0002")
    assert result is not None


def test_gse64349_resolves_gene_symbol_via_gse77738_lookup(tmp_path: Path) -> None:
    """'cdc6_1' in GSE64349's Feature ID column must resolve to MA_0001 via GSE77738's Gene Name."""
    _seed_gse77738(tmp_path)
    _seed_gse64349(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")

    bundle = load_gse64349_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)

    assert "MA_0001" in bundle.known_tags
    assert bundle.lookup("MA_0001", "MA_0002") is not None


def test_result_is_cached_and_reused_without_reparsing(tmp_path: Path) -> None:
    """A second call for the same query must not need the (possibly large) source file again."""
    _seed_gse77738(tmp_path)
    cache = JsonCache(tmp_path / "jsoncache")
    first = load_gse77738_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)
    first_value = first.lookup("MA_0001", "MA_0002")
    assert first_value is not None

    (tmp_path / "coexpression" / "GSE77738_ReadCounts.xls").unlink()

    second = load_gse77738_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)
    assert second.lookup("MA_0001", "MA_0002") == first_value
    assert second.warnings == ()


def test_download_failure_degrades_to_empty_bundle_without_raising(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A fully offline environment must degrade gracefully, never raise."""
    cache = JsonCache(tmp_path / "jsoncache")
    monkeypatch.setattr(
        "analysis.coexpression_bridge.requests.get",
        lambda *a, **k: (_ for _ in ()).throw(requests.RequestException("simulated network failure")),
    )

    bundle = load_gse77738_coexpression_bundle(True, ["MA_0001"], cache, tmp_path)

    assert bundle.lookup("MA_0001", "MA_0002") is None
    assert bundle.warnings  # at least one warning explaining the degradation
