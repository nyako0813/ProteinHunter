"""Tests for tools/scoring_recalibration.py (tier thresholds and leave-one-out comparison)."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from analysis.scoring_engine_config import TierThresholds
from tools import calibration_report as cr
from tools import scoring_recalibration as sr


def sample(pos_scores, neg_scores, groups=None, *, categories: int = 3) -> sr.Sample:
    """A Sample with one row per score; each positive is its own fold unless ``groups`` says otherwise."""
    groups = groups or [f"g{i}" for i in range(len(pos_scores))]
    positives = pd.DataFrame(
        {
            "group": groups,
            "query": [f"Q{i}" for i in range(len(pos_scores))],
            "candidate": [f"P{i}" for i in range(len(pos_scores))],
            "candidate_source": "Candidates",
            "candidate_rank": 1.0,
            "final_score": [float(v) for v in pos_scores],
            "interaction_score": [float(v) for v in pos_scores],
            "evidence_category_count": float(categories),
        }
    )
    negatives = pd.DataFrame(
        {
            "query": "Q0",
            "candidate": [f"N{i}" for i in range(len(neg_scores))],
            "candidate_source": "Candidates",
            "candidate_rank": 1.0,
            "final_score": [float(v) for v in neg_scores],
            "interaction_score": [float(v) for v in neg_scores],
            "evidence_category_count": float(categories),
        }
    )
    return sr.Sample(positives, negatives)


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------


def test_auc_counts_ties_half_and_is_nan_for_an_empty_side() -> None:
    assert sr.auc([5, 6, 7], [1, 2, 6]) == pytest.approx((2 + 2.5 + 3) / 9)  # 5 beats 1,2; 6 beats 1,2 and ties 6; 7 beats all three
    assert math.isnan(sr.auc([], [1.0]))


def test_percentile_among_negatives() -> None:
    assert sr.percentile_among(5, [1, 2, 6]) == pytest.approx(2 / 3)
    assert sr.percentile_among(2, [1, 2, 6]) == pytest.approx((1 + 0.5) / 3)


def test_youden_finds_the_best_cut_and_reports_ties() -> None:
    # t=5: TPR 1, FPR 1/3 -> J 2/3 (best); t=6: 2/3 - 1/3; t=7: 1/3 - 0.
    j, thresholds = sr.youden([5, 6, 7], [1, 2, 6])
    assert j == pytest.approx(2 / 3) and thresholds == [5.0]
    j, cut, _ = sr.youden_cut([5, 6, 7], [1, 2, 6])
    assert cut == pytest.approx((2 + 5) / 2)  # midway to the next lower observed score

    _, tied = sr.youden([2, 4], [1, 3])  # t=2: 1-1... perfect ties produce more than one maximiser
    assert len(tied) >= 1


def test_jackknife_standard_error_matches_the_textbook_formula() -> None:
    # replicates 1,2,3: (k-1)/k * sum((x-mean)^2) = 2/3 * 2
    assert sr.jackknife_se([1, 2, 3]) == pytest.approx(math.sqrt(4 / 3))
    assert math.isnan(sr.jackknife_se([1.0]))


# ---------------------------------------------------------------------------
# Leave-one-out comparison
# ---------------------------------------------------------------------------


def test_identical_variants_are_never_a_clear_improvement() -> None:
    base = sample([10, 20, 30, 40], [15, 25, 35])
    twin = sample([10, 20, 30, 40], [15, 25, 35])

    analysis = sr.loo_analysis({"base": base, "twin": twin}, "base")

    twin_result = next(v for v in analysis.variants if v.name == "twin")
    assert twin_result.delta_auc == 0 and twin_result.z == 0 and not twin_result.clear_improvement
    assert (twin_result.folds_better, twin_result.folds_worse, twin_result.folds_equal) == (0, 0, 4)


def test_hand_computed_leave_one_out_values() -> None:
    base = sample([10, 20, 30, 40], [15, 25, 35])
    better = sample([110, 120, 130, 140], [15, 25, 35])

    analysis = sr.loo_analysis({"base": base, "better": better}, "base")

    b, v = analysis.variants
    assert b.full_auc == pytest.approx(0.5) and v.full_auc == pytest.approx(1.0)
    # baseline leave-one-out AUCs: 6/9, 5/9, 4/9, 3/9 -> mean 0.5; held-out percentiles 0, 1/3, 2/3, 1 -> mean 0.5 (= full AUC)
    assert b.loo_mean_auc == pytest.approx(0.5) and b.held_out_percentile == pytest.approx(0.5)
    assert v.loo_mean_auc == pytest.approx(1.0) and v.held_out_percentile == pytest.approx(1.0)
    assert v.delta_auc == pytest.approx(0.5)
    # held-out percentiles: baseline 0, 1/3, 2/3, 1 vs candidate 1, 1, 1, 1 -> better in 3 folds, equal in the last (both perfect)
    assert (v.folds_better, v.folds_equal, v.folds_worse) == (3, 1, 0)
    # two-sample jackknife: positives 3/4*0.06173 + negatives 2/3*0.03125 = 0.06713 -> se 0.2591, z 1.93
    assert v.delta_se == pytest.approx(math.sqrt(0.75 * (2 * (1 / 6) ** 2 + 2 * (1 / 18) ** 2) + (2 / 3) * (2 * 0.125**2)))
    assert v.z == pytest.approx(0.5 / v.delta_se) and v.z < sr.CLEAR_Z
    assert not v.clear_improvement  # a big gain on four pairs is still not enough evidence
    # every fold prefers the higher-AUC variant when choosing on the remaining pairs
    assert analysis.chosen_counts == {"base": 0, "better": 4}
    assert analysis.nested_selection_percentile == pytest.approx(1.0) and analysis.baseline_percentile == pytest.approx(0.5)


def test_a_large_consistent_gain_is_flagged_clear() -> None:
    rng = np.random.default_rng(0)
    pos = rng.normal(40, 8, 12)
    neg = rng.normal(35, 8, 15)
    base = sample(pos, neg)
    better = sample(pos + 60, neg)  # every positive moves far above every negative

    analysis = sr.loo_analysis({"base": base, "better": better}, "base")

    v = analysis.variants[1]
    assert v.delta_auc > 0.2 and v.z >= sr.CLEAR_Z and v.folds_better >= 10 and v.clear_improvement


def test_a_worse_variant_is_never_clear_and_nested_selection_prefers_the_baseline() -> None:
    base = sample([50, 60, 70, 80, 90], [10, 20, 30, 40, 55])
    worse = sample([5, 15, 25, 35, 45], [10, 20, 30, 40, 55])

    analysis = sr.loo_analysis({"base": base, "worse": worse}, "base")

    w = analysis.variants[1]
    assert w.delta_auc < 0 and not w.clear_improvement and w.folds_worse > 0
    assert analysis.chosen_counts["base"] == 5


def test_a_fold_leaves_out_both_directions_of_a_pair() -> None:
    """NifD->NifK and NifK->NifD share a group, so they leave together (no leakage between the two directions)."""
    base = sample([10, 20, 30], [15], groups=["A|B", "A|B", "C|D"])

    analysis = sr.loo_analysis({"base": base}, "base")

    assert analysis.n_groups == 2 and analysis.n_positive_rows == 3
    assert list(analysis.fold_table["held_out_pair"]) == ["A|B", "C|D"]
    # leaving A|B out leaves only the 30 (beats 15): AUC 1; leaving C|D out leaves 10, 20: AUC 0.5
    assert list(analysis.fold_table["base_loo_auc"]) == pytest.approx([1.0, 0.5])


def test_misaligned_samples_are_rejected_and_align_samples_fixes_them() -> None:
    a = sample([10, 20, 30], [5, 6])
    b = sample([10, 20, 30], [5, 6, 7])

    with pytest.raises(ValueError, match="not aligned"):
        sr.loo_analysis({"a": a, "b": b}, "a")

    aligned, notes = sr.align_samples({"a": a, "b": b}, "final_score")
    assert len(aligned["b"].negatives) == 2 and any("b:" in note for note in notes)
    sr.loo_analysis(aligned, "a")


def test_unknown_baseline_is_rejected() -> None:
    with pytest.raises(ValueError, match="baseline"):
        sr.loo_analysis({"a": sample([1, 2], [0])}, "zzz")


# ---------------------------------------------------------------------------
# Tier thresholds
# ---------------------------------------------------------------------------


def test_classify_uses_the_engines_own_rule() -> None:
    thresholds = TierThresholds(tier3_min_score=25.0)  # pinned: independent of the shipped defaults

    assert sr.classify(72, 3, thresholds) == "Tier1_VeryStrong"
    assert sr.classify(72, 2, thresholds) == "Tier2_Strong"  # too few categories for Tier1
    assert sr.classify(55, 2, thresholds) == "Tier2_Strong"
    assert sr.classify(55, 1, thresholds) == "Tier3_Moderate"
    assert sr.classify(30, 1, thresholds) == "Tier3_Moderate"
    assert sr.classify(10, 5, thresholds) == "Tier4_Weak"
    assert sr.classify(float("nan"), 3, thresholds) == "Unclassified"


def test_tier_breakdown_counts_both_groups_under_given_thresholds() -> None:
    data = sample([60, 45, 20], [30, 12, 55])

    current = sr.tier_breakdown(data, TierThresholds(tier3_min_score=25.0))
    lowered = sr.tier_breakdown(data, TierThresholds(tier2_min_score=40.0, tier3_min_score=25.0))

    assert (current.loc["Tier2_Strong", "Tier A"], current.loc["Tier2_Strong", "Negatives"]) == (1, 1)
    assert (lowered.loc["Tier2_Strong", "Tier A"], lowered.loc["Tier2_Strong", "Negatives"]) == (2, 1)
    assert current.loc["Tier3_Moderate", "Negatives"] == 1 and lowered.loc["Tier4_Weak", "Tier A"] == 1


def test_safety_net_effect_counts_rows_shown_beyond_the_top_n() -> None:
    scores = [90, 80, 70, 60, 55, 45, 40, 35, 30, 20]
    table = pd.DataFrame(
        {
            "query_old_locus_tag": "Q1",
            "candidate_old_locus_tag": [f"C{i}" for i in range(10)],
            "candidate_rank": list(range(1, 11)),
            "final_score": scores,
            "evidence_category_count": 3,
        }
    )

    current = sr.safety_net_effect(table, TierThresholds(), max_per_query=3)
    lowered = sr.safety_net_effect(table, TierThresholds(tier2_min_score=35.0), max_per_query=3)

    # Tier1: 90,80,70 (>=70). Tier2 (>=50): 60,55. So 5 Tier1/2 rows, 3 in the top 3 -> 2 more are added.
    assert current.iloc[0].to_dict() == {"query": "Q1", "tier1_2_rows": 5, "shown": 5, "shown_beyond_top": 2}
    # Tier2 from 35: adds 45, 40, 35 -> 8 rows, 5 beyond the top 3.
    assert lowered.iloc[0].to_dict() == {"query": "Q1", "tier1_2_rows": 8, "shown": 8, "shown_beyond_top": 5}


# ---------------------------------------------------------------------------
# Loading, matching and the CLI
# ---------------------------------------------------------------------------


def make_table(shift: float = 0.0) -> pd.DataFrame:
    rows = [
        ("Q1", "P1", 60), ("Q1", "P2", 55), ("Q2", "Q1", 50),  # Tier A pairs (Q1-P1, Q1-P2, Q1-Q2 via Q2 -> Q1)
        ("Q1", "N1", 20), ("Q1", "N2", 30), ("Q1", "N3", 25),  # negatives for query Q1
    ]
    return pd.DataFrame(
        {
            "query_old_locus_tag": [r[0] for r in rows],
            "candidate_old_locus_tag": [r[1] for r in rows],
            "query_id": [r[0] for r in rows],
            "candidate_source": "Candidates",
            "candidate_rank": [1, 2, 3, 4, 5, 6],
            "final_score": [r[2] + shift for r in rows],
            "final_score_tier": "Tier3_Moderate",
            "evidence_category_count": 3,
            "interaction_score": [r[2] + shift for r in rows],
        }
    )


def curated_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "protein_a_old_locus_tag": ["Q1", "Q1", "Q1", "Q1"],
            "protein_b_old_locus_tag": ["P1", "P2", "Q2", "X9"],
            "tier": ["A", "A", "A", "B"],
            "protein_a_label": "", "protein_b_label": "", "class": "", "confidence": "", "source": "",
        }
    )


def test_build_sample_matches_tier_a_pairs_in_both_directions_and_the_negatives() -> None:
    negatives = pd.DataFrame({"old_locus_tag": ["N1", "N2", "N3", "N4"], "query_old_locus_tag": "Q1", "source": "", "af3_classification": "", "af3_ipTM": ""})

    built = sr.build_sample(make_table(), curated_frame(), negatives)

    assert sorted(zip(built.positives["query"], built.positives["candidate"])) == [("Q1", "P1"), ("Q1", "P2"), ("Q2", "Q1")]
    assert built.positives.set_index("candidate")["final_score"].to_dict() == {"P1": 60.0, "P2": 55.0, "Q1": 50.0}
    assert sorted(built.negatives["candidate"]) == ["N1", "N2", "N3"]  # N4 is not in the results
    assert set(built.positives["group"]) == {"P1|Q1", "P2|Q1", "Q1|Q2"}
    assert built.positives["evidence_category_count"].eq(3).all()


def test_load_table_rejects_a_csv_without_the_required_columns(tmp_path: Path) -> None:
    path = tmp_path / "bad.csv"
    path.write_text("a,b\n1,2\n")

    with pytest.raises(cr.CalibrationInputError, match="missing column"):
        sr.load_table(path)


def test_cli_writes_the_summary_and_fold_table(tmp_path: Path) -> None:
    make_table().to_csv(tmp_path / "base.csv", index=False)
    make_table(shift=5.0).to_csv(tmp_path / "up.csv", index=False)
    curated_frame().to_csv(tmp_path / "curated.csv", index=False)
    pd.DataFrame({"old_locus_tag": ["N1", "N2", "N3"], "source": "x"}).to_csv(tmp_path / "neg.csv", index=False)

    code = sr.main(
        [
            "--curated", str(tmp_path / "curated.csv"), "--negatives", str(tmp_path / "neg.csv"),
            "--baseline", f"base={tmp_path / 'base.csv'}", "--candidate", f"up={tmp_path / 'up.csv'}",
            "--negatives-query", "Q1", "--tier2", "40", "--out", str(tmp_path / "out"),
        ]
    )

    assert code == 0
    summary = (tmp_path / "out" / "recalibration_summary.md").read_text(encoding="utf-8")
    assert "## Phase 1" in summary and "## Phase 2" in summary and "### final_score" in summary and "### interaction_score" in summary
    assert "proposed" in summary and "Youden's J" in summary
    folds = pd.read_excel(tmp_path / "out" / "loo_folds.xlsx", sheet_name="final_score")
    assert list(folds["held_out_pair"]) == ["P1|Q1", "P2|Q1", "Q1|Q2"]
    assert sr.main(["--curated", str(tmp_path / "curated.csv"), "--negatives", str(tmp_path / "neg.csv"), "--baseline", f"base={tmp_path / 'base.csv'}", "--negatives-query", "Q1", "--out", str(tmp_path / "out")]) == 2  # exists, no --force
