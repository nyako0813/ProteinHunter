"""Recalibration analysis for the v2 scoring model: tier thresholds and cap/weight candidates.

Analysis only -- nothing here changes scoring. It reads per-variant tables of
scored query/candidate rows (one variant = one pipeline run with different caps
or weights), matches them to the curated Tier A pairs and the known negatives
exactly as tools/calibration_report.py does, and reports:

Phase 1  -- ``final_score`` of Tier A vs negatives, the tier each lands in under
            the current thresholds, Youden's J thresholds, and how a threshold
            change alters the always-shown Tier1/Tier2 set of the reports.
Phase 2  -- leave-one-out cross-validation over a *small, pre-specified* list of
            candidate settings against the baseline. No grid search: with about
            ten positive pairs an optimiser would mostly fit noise. A candidate
            is only flagged as a clear improvement if its AUC gain over the
            baseline is large relative to its jackknife standard error AND the
            gain holds in most folds; otherwise the honest result is "keep the
            current values".

Input tables (CSV) need the columns
``query_old_locus_tag, candidate_old_locus_tag, candidate_source, query_id,
candidate_rank, final_score, final_score_tier, evidence_category_count,
interaction_score`` (a ``04_Score_Breakdown``-like export with locus tags), or
a ProteinHunter results workbook (read through tools/calibration_report.py).

Usage
  python tools/scoring_recalibration.py --curated claude/experimental_interactions_curated.csv \
      --negatives claude/calibration/af3_negatives_MA_4115.csv \
      --baseline baseline=rows_baseline.csv \
      --candidate A=rows_A.csv --candidate B=rows_B.csv \
      --out claude/calibration/2026-09-21_recalibration
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from tools import calibration_report as cr

#: A candidate must beat the baseline by at least this many jackknife standard errors ...
CLEAR_Z = 2.4  # about a two-sided 5% level after allowing for three candidate comparisons
#: ... and by at least this AUC margin ...
CLEAR_MIN_DELTA = 0.02
#: ... and be better in at least this share of the folds where the two differ.
CLEAR_FOLD_SHARE = 0.75

TABLE_COLUMNS = (
    "query_old_locus_tag",
    "candidate_old_locus_tag",
    "candidate_source",
    "candidate_rank",
    "final_score",
    "final_score_tier",
    "evidence_category_count",
    "interaction_score",
)


# ---------------------------------------------------------------------------
# Loading and matching
# ---------------------------------------------------------------------------


def load_table(path: Path) -> pd.DataFrame:
    """A variant's scored rows from a CSV export or a results workbook."""
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xlsm"}:
        pairs = cr.read_results(path).pairs
        return pairs.rename(columns={"query_tag": "query_old_locus_tag", "candidate_tag": "candidate_old_locus_tag"})
    frame = pd.read_csv(path)
    missing = [column for column in TABLE_COLUMNS if column not in frame.columns]
    if missing:
        raise cr.CalibrationInputError(f"{path}: missing column(s) {missing}")
    return frame


def _as_results(table: pd.DataFrame) -> cr.Results:
    pairs = table.rename(columns={"query_old_locus_tag": "query_tag", "candidate_old_locus_tag": "candidate_tag"}).copy()
    for column in ("final_score", "interaction_score", "interaction_priority_score", "candidate_rank"):
        pairs[column] = pd.to_numeric(pairs[column], errors="coerce") if column in pairs.columns else np.nan
    return cr.Results(pairs=pairs, tag_by_protein={})


@dataclass(frozen=True)
class Sample:
    """The labelled pairs of one variant: positives carry a group (their unordered pair) so folds can leave a pair out."""

    positives: pd.DataFrame  # columns: group, query, candidate, <scores>
    negatives: pd.DataFrame  # columns: query, candidate, <scores>


def build_sample(table: pd.DataFrame, curated: pd.DataFrame, negatives: pd.DataFrame) -> Sample:
    """Tier A rows and negatives found in ``table`` (same matching rules as tools/calibration_report.py).

    ``negatives`` is the frame from ``cr.read_negatives`` (it already carries each negative's query).
    """
    results = _as_results(table)
    queries_present = set(results.pairs["query_tag"].dropna())
    fields = (
        "candidate_source", "candidate_rank", "final_score", "interaction_score", "interaction_priority_score",
        "evidence_category_count", "final_score_tier", "evidence_tier", "interaction_evidence_tier",
    )

    def record(row: pd.Series) -> dict[str, Any]:
        return {name: row.get(name, np.nan) for name in fields}

    positive_rows: list[dict[str, Any]] = []
    for _, item in curated[curated["tier"] == "A"].iterrows():
        a, b = item["protein_a_old_locus_tag"].strip(), item["protein_b_old_locus_tag"].strip()
        for query, candidate in ((a, b), (b, a)):
            if "/" in candidate or query not in queries_present:
                continue
            row = cr._lookup(results, query, candidate)
            if row is not None:
                positive_rows.append({"group": "|".join(sorted((a, b))), "query": query, "candidate": candidate, **record(row)})

    negative_rows: list[dict[str, Any]] = []
    for _, item in negatives.iterrows():
        query, candidate = item["query_old_locus_tag"].strip(), item["old_locus_tag"].strip()
        row = cr._lookup(results, query, candidate)
        if row is not None:
            negative_rows.append({"query": query, "candidate": candidate, **record(row)})

    def frame(rows: list[dict[str, Any]], columns: list[str]) -> pd.DataFrame:
        out = pd.DataFrame(rows, columns=columns)
        for column in ("final_score", "interaction_score", "interaction_priority_score", "evidence_category_count", "candidate_rank"):
            out[column] = pd.to_numeric(out[column], errors="coerce")
        return out

    common = ["query", "candidate", *fields]
    return Sample(frame(positive_rows, ["group", *common]), frame(negative_rows, common))


def align_samples(samples: dict[str, Sample], metric: str) -> tuple[dict[str, Sample], list[str]]:
    """Restrict every variant to the positives/negatives that are scored under all of them."""
    notes: list[str] = []
    pos_keys = None
    neg_keys = None
    for name, sample in samples.items():
        p = {(g, q, c) for g, q, c, v in zip(sample.positives["group"], sample.positives["query"], sample.positives["candidate"], sample.positives[metric]) if pd.notna(v)}
        n = {(q, c) for q, c, v in zip(sample.negatives["query"], sample.negatives["candidate"], sample.negatives[metric]) if pd.notna(v)}
        pos_keys = p if pos_keys is None else pos_keys & p
        neg_keys = n if neg_keys is None else neg_keys & n
    aligned: dict[str, Sample] = {}
    for name, sample in samples.items():
        pos = sample.positives[[(g, q, c) in pos_keys for g, q, c in zip(sample.positives["group"], sample.positives["query"], sample.positives["candidate"])]]
        neg = sample.negatives[[(q, c) in neg_keys for q, c in zip(sample.negatives["query"], sample.negatives["candidate"])]]
        pos = pos.sort_values(["group", "query", "candidate"]).reset_index(drop=True)
        neg = neg.sort_values(["query", "candidate"]).reset_index(drop=True)
        aligned[name] = Sample(pos, neg)
        if len(pos) != len(sample.positives) or len(neg) != len(sample.negatives):
            notes.append(f"{name}: kept {len(pos)}/{len(sample.positives)} positives and {len(neg)}/{len(sample.negatives)} negatives scored under every variant")
    return aligned, notes


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def auc(positive: Sequence[float], negative: Sequence[float]) -> float:
    """P(positive > negative) with ties counted half; NaN when either side is empty."""
    p = np.asarray(positive, dtype=float)
    n = np.asarray(negative, dtype=float)
    if len(p) == 0 or len(n) == 0:
        return float("nan")
    greater = (p[:, None] > n[None, :]).sum()
    ties = (p[:, None] == n[None, :]).sum()
    return float((greater + 0.5 * ties) / (len(p) * len(n)))


def percentile_among(value: float, negative: Sequence[float]) -> float:
    """Share of the negatives that ``value`` exceeds (ties half): where one held-out positive ranks."""
    n = np.asarray(negative, dtype=float)
    return float(((n < value).sum() + 0.5 * (n == value).sum()) / len(n))


def youden(positive: Sequence[float], negative: Sequence[float]) -> tuple[float, list[float]]:
    """Best Youden's J = TPR - FPR for the rule "score >= t", and every observed t reaching it."""
    p = np.asarray(positive, dtype=float)
    n = np.asarray(negative, dtype=float)
    best = -math.inf
    thresholds: list[float] = []
    for t in sorted(set(p.tolist()) | set(n.tolist())):
        j = float((p >= t).mean() - (n >= t).mean())
        if j > best + 1e-12:
            best, thresholds = j, [t]
        elif abs(j - best) <= 1e-12:
            thresholds.append(t)
    return best, thresholds


def youden_cut(positive: Sequence[float], negative: Sequence[float]) -> tuple[float, float, list[float]]:
    """(J, cut, thresholds): the cut sits midway between the best threshold and the next lower observed score."""
    j, thresholds = youden(positive, negative)
    best = thresholds[0]
    values = sorted(set(np.asarray(positive, dtype=float).tolist()) | set(np.asarray(negative, dtype=float).tolist()))
    below = [v for v in values if v < best]
    cut = (below[-1] + best) / 2 if below else best
    return j, cut, thresholds


def jackknife_se(replicates: Sequence[float]) -> float:
    """Jackknife standard error of a statistic from its leave-one-out replicates."""
    values = np.asarray(replicates, dtype=float)
    k = len(values)
    if k < 2:
        return float("nan")
    return float(math.sqrt((k - 1) / k * ((values - values.mean()) ** 2).sum()))


@dataclass(frozen=True)
class VariantResult:
    name: str
    full_auc: float
    loo_mean_auc: float  # mean over the leave-one-pair-out folds of the AUC on the remaining pairs
    loo_se: float  # jackknife standard error of that AUC
    held_out_percentile: float  # mean over folds of where the held-out pair ranks among the negatives
    delta_auc: float  # full AUC minus the baseline's
    delta_se: float  # two-sample jackknife standard error of that difference
    z: float
    folds_better: int
    folds_worse: int
    folds_equal: int
    clear_improvement: bool


@dataclass(frozen=True)
class LooAnalysis:
    metric: str
    n_positive_rows: int
    n_groups: int
    n_negatives: int
    variants: list[VariantResult]
    fold_table: pd.DataFrame  # one row per fold (held-out pair): fold AUC and held-out percentile per variant
    nested_selection_percentile: float  # mean held-out percentile when each fold picks its variant on the remaining pairs
    baseline_percentile: float
    chosen_counts: dict[str, int]


def _fold_stats(sample: Sample, metric: str) -> tuple[list[str], np.ndarray, dict[str, np.ndarray]]:
    groups = sample.positives["group"].tolist()
    return sorted(set(groups)), sample.positives[metric].to_numpy(float), {}


def loo_analysis(samples: dict[str, Sample], baseline: str, metric: str = "final_score") -> LooAnalysis:
    """Leave-one-pair-out comparison of every variant with ``baseline`` (all samples must be aligned)."""
    names = list(samples)
    if baseline not in samples:
        raise ValueError(f"baseline {baseline!r} is not among the variants {names}")
    reference = samples[baseline]
    groups_of_rows = reference.positives["group"].tolist()
    group_names = sorted(set(groups_of_rows))
    for name in names:
        if samples[name].positives["group"].tolist() != groups_of_rows or len(samples[name].negatives) != len(reference.negatives):
            raise ValueError("samples are not aligned (use align_samples)")

    pos = {name: samples[name].positives[metric].to_numpy(float) for name in names}
    neg = {name: samples[name].negatives[metric].to_numpy(float) for name in names}
    group_array = np.asarray(groups_of_rows)

    fold_auc = {name: [] for name in names}
    fold_percentile = {name: [] for name in names}
    for group in group_names:
        keep = group_array != group
        held = ~keep
        for name in names:
            fold_auc[name].append(auc(pos[name][keep], neg[name]))
            fold_percentile[name].append(float(np.mean([percentile_among(v, neg[name]) for v in pos[name][held]])))

    n_negatives = len(neg[baseline])
    results: list[VariantResult] = []
    base_full = auc(pos[baseline], neg[baseline])
    for name in names:
        full = auc(pos[name], neg[name])
        # Two-sample jackknife of the AUC difference: leave out each positive pair, then each negative.
        delta_pos = [fold_auc[name][i] - fold_auc[baseline][i] for i in range(len(group_names))]
        delta_neg = []
        for j in range(n_negatives):
            mask = np.arange(n_negatives) != j
            delta_neg.append(auc(pos[name], neg[name][mask]) - auc(pos[baseline], neg[baseline][mask]))
        k, m = len(delta_pos), n_negatives
        variance = 0.0
        if k > 1:
            variance += (k - 1) / k * float(((np.asarray(delta_pos) - np.mean(delta_pos)) ** 2).sum())
        if m > 1:
            variance += (m - 1) / m * float(((np.asarray(delta_neg) - np.mean(delta_neg)) ** 2).sum())
        se = math.sqrt(variance)
        delta = full - base_full
        z = delta / se if se > 0 else (math.inf if delta > 0 else (-math.inf if delta < 0 else 0.0))
        diffs = np.asarray(fold_percentile[name]) - np.asarray(fold_percentile[baseline])
        better, worse = int((diffs > 1e-12).sum()), int((diffs < -1e-12).sum())
        equal = len(diffs) - better - worse
        differing = better + worse
        clear = (
            name != baseline
            and delta >= CLEAR_MIN_DELTA
            and z >= CLEAR_Z
            and differing > 0
            and better / differing >= CLEAR_FOLD_SHARE
        )
        results.append(
            VariantResult(
                name=name,
                full_auc=full,
                loo_mean_auc=float(np.mean(fold_auc[name])),
                loo_se=jackknife_se(fold_auc[name]),
                held_out_percentile=float(np.mean(fold_percentile[name])),
                delta_auc=delta,
                delta_se=se,
                z=z,
                folds_better=better,
                folds_worse=worse,
                folds_equal=equal,
                clear_improvement=clear,
            )
        )

    # Nested selection: on each fold pick the variant with the best AUC on the remaining pairs
    # (ties go to the baseline), then score the held-out pair with that choice.
    chosen: dict[str, int] = {name: 0 for name in names}
    nested: list[float] = []
    for i in range(len(group_names)):
        best = max(fold_auc[name][i] for name in names)
        pick = baseline if fold_auc[baseline][i] >= best - 1e-12 else next(n for n in names if fold_auc[n][i] >= best - 1e-12)
        chosen[pick] += 1
        nested.append(fold_percentile[pick][i])

    fold_table = pd.DataFrame({"held_out_pair": group_names})
    for name in names:
        fold_table[f"{name}_loo_auc"] = fold_auc[name]
        fold_table[f"{name}_held_out_percentile"] = fold_percentile[name]
    return LooAnalysis(
        metric=metric,
        n_positive_rows=len(groups_of_rows),
        n_groups=len(group_names),
        n_negatives=n_negatives,
        variants=results,
        fold_table=fold_table,
        nested_selection_percentile=float(np.mean(nested)),
        baseline_percentile=float(np.mean(fold_percentile[baseline])),
        chosen_counts=chosen,
    )


# ---------------------------------------------------------------------------
# Phase 1: tier thresholds
# ---------------------------------------------------------------------------


def classify(score: float | None, categories: float | int | None, thresholds: Any) -> str:
    """The engine's own tier rule (analysis.scoring_engine._classify_tier) for a score and category count."""
    from analysis.scoring_engine import _classify_tier

    if score is None or (isinstance(score, float) and math.isnan(score)):
        return "Unclassified"
    count = 0 if categories is None or (isinstance(categories, float) and math.isnan(categories)) else int(categories)
    return _classify_tier(float(score), count, thresholds)


def tier_breakdown(sample: Sample, thresholds: Any) -> pd.DataFrame:
    """Rows = tiers, columns = how many Tier A pairs and negatives fall in each under ``thresholds``."""
    labels = ["Tier1_VeryStrong", "Tier2_Strong", "Tier3_Moderate", "Tier4_Weak", "Unclassified"]
    counts = {"Tier A": dict.fromkeys(labels, 0), "Negatives": dict.fromkeys(labels, 0)}
    for label, frame in (("Tier A", sample.positives), ("Negatives", sample.negatives)):
        for score, cats in zip(frame["final_score"], frame["evidence_category_count"]):
            counts[label][classify(score, cats, thresholds)] += 1
    return pd.DataFrame(counts, index=labels)


def safety_net_effect(table: pd.DataFrame, thresholds: Any, max_per_query: int = 15) -> pd.DataFrame:
    """Per query: rows in Tier1/Tier2 and how many of them the reports show *beyond* the top ``max_per_query``.

    Uses output.report_v2.select_top_candidates_per_query, the function the Word/Notion reports
    and the Excel word_report_link column use, so the count is what a report would actually add.
    """
    from output.report_v2 import select_top_candidates_per_query

    rows = table.to_dict("records")
    for row in rows:
        row["final_score_tier"] = classify(row.get("final_score"), row.get("evidence_category_count"), thresholds)
        row["query_id"] = row.get("query_old_locus_tag") or row.get("query_id")
    selected = select_top_candidates_per_query(rows, max_per_query)
    by_query: dict[str, dict[str, int]] = {}
    for row in rows:
        query = str(row["query_id"])
        entry = by_query.setdefault(query, {"tier1_2_rows": 0, "shown": 0, "shown_beyond_top": 0})
        if row["final_score_tier"] in {"Tier1_VeryStrong", "Tier2_Strong"}:
            entry["tier1_2_rows"] += 1
    for row in selected:
        entry = by_query[str(row["query_id"])]
        entry["shown"] += 1
        if (row.get("candidate_rank") or 0) > max_per_query:
            entry["shown_beyond_top"] += 1
    frame = pd.DataFrame(by_query).T.reset_index().rename(columns={"index": "query"})
    return frame.sort_values("query").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def _f(value: float, digits: int = 3) -> str:
    return "n/a" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:.{digits}f}"


def render_loo(analysis: LooAnalysis, baseline: str) -> list[str]:
    lines = [
        f"Metric: `{analysis.metric}`. {analysis.n_positive_rows} Tier A rows in {analysis.n_groups} folds (one fold = one unordered pair; "
        f"both directions of a pair leave together) vs {analysis.n_negatives} negatives.",
        "",
        "| variant | full AUC | LOO AUC (mean ± jackknife SE) | held-out percentile | ΔAUC vs baseline (SE) | z | folds better / equal / worse | clear improvement |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in analysis.variants:
        is_base = r.name == baseline
        delta = "baseline" if is_base else f"{r.delta_auc:+.3f} ({_f(r.delta_se)})"
        z = "" if is_base else _f(r.z, 2)
        folds = "" if is_base else f"{r.folds_better} / {r.folds_equal} / {r.folds_worse}"
        verdict = "" if is_base else ("**yes**" if r.clear_improvement else "no")
        lines.append(f"| {r.name} | {_f(r.full_auc)} | {_f(r.loo_mean_auc)} ± {_f(r.loo_se)} | {_f(r.held_out_percentile)} | {delta} | {z} | {folds} | {verdict} |")
    chosen = ", ".join(f"{k}: {v}" for k, v in analysis.chosen_counts.items() if v)
    lines += [
        "",
        f"Nested selection (each fold picks its best variant on the other pairs, ties to the baseline): held-out percentile "
        f"{_f(analysis.nested_selection_percentile)} vs baseline {_f(analysis.baseline_percentile)}; picks per fold: {chosen}.",
        "",
        f"A variant is called a clear improvement only if ΔAUC ≥ {CLEAR_MIN_DELTA}, z ≥ {CLEAR_Z} and it is better in at least "
        f"{int(CLEAR_FOLD_SHARE * 100)}% of the folds where it differs. Nothing is fitted: the variants were fixed in advance.",
        "",
    ]
    return lines


def render_phase1(baseline_sample: Sample, baseline_table: pd.DataFrame, current: Any, proposed: Any | None, max_per_query: int) -> list[str]:
    pos = baseline_sample.positives["final_score"].dropna().to_numpy(float)
    neg = baseline_sample.negatives["final_score"].dropna().to_numpy(float)
    j, cut, thresholds = youden_cut(pos, neg)
    lines = [
        f"`final_score` of Tier A (n={len(pos)}): min {pos.min():.1f}, median {np.median(pos):.1f}, max {pos.max():.1f}; "
        f"negatives (n={len(neg)}): min {neg.min():.1f}, median {np.median(neg):.1f}, max {neg.max():.1f}. AUC {_f(auc(pos, neg))}.",
        "",
        f"Youden's J is maximal ({j:.3f}) for the rule `final_score >= {thresholds[0]:.1f}` (cut midway to the next lower observed score: {cut:.1f}).",
        "",
        "Sensitivity (share of Tier A at or above t) and false-positive rate (share of negatives at or above t):",
        "",
        "| t | Tier A ≥ t | negatives ≥ t | J |",
        "|---:|---:|---:|---:|",
    ]
    for t in range(20, 71, 5):
        tpr, fpr = float((pos >= t).mean()), float((neg >= t).mean())
        lines.append(f"| {t} | {tpr:.2f} | {fpr:.2f} | {tpr - fpr:+.2f} |")
    lines += ["", "Tier of every Tier A pair and negative under the **current** thresholds:", ""]
    lines += ["| tier | Tier A | negatives |", "|---|---:|---:|"]
    for label, row in tier_breakdown(baseline_sample, current).iterrows():
        lines.append(f"| {label} | {int(row['Tier A'])} | {int(row['Negatives'])} |")
    if proposed is not None:
        lines += ["", "…and under the **proposed** thresholds:", "", "| tier | Tier A | negatives |", "|---|---:|---:|"]
        for label, row in tier_breakdown(baseline_sample, proposed).iterrows():
            lines.append(f"| {label} | {int(row['Tier A'])} | {int(row['Negatives'])} |")
    for title, thresholds_obj in (("current", current), ("proposed", proposed)):
        if thresholds_obj is None:
            continue
        effect = safety_net_effect(baseline_table, thresholds_obj, max_per_query)
        lines += ["", f"Always-shown safety net (Tier1/Tier2 rows shown regardless of rank, top {max_per_query} per query) under the **{title}** thresholds:", ""]
        lines += ["| query | Tier1/2 rows | shown | shown beyond top " + str(max_per_query) + " |", "|---|---:|---:|---:|"]
        for _, r in effect.iterrows():
            lines.append(f"| {r['query']} | {int(r['tier1_2_rows'])} | {int(r['shown'])} | {int(r['shown_beyond_top'])} |")
    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_named(values: Iterable[str]) -> dict[str, Path]:
    named: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise cr.CalibrationInputError(f"expected NAME=PATH, got {value!r}")
        name, path = value.split("=", 1)
        named[name.strip()] = Path(path.strip())
    return named


def run(
    curated: Path,
    negatives: Path,
    baseline: tuple[str, Path],
    candidates: dict[str, Path],
    out: Path,
    *,
    proposed_thresholds: Any | None = None,
    max_per_query: int = 15,
    negatives_query: str = "MA_4115",
    force: bool = False,
) -> dict[str, Path]:
    from analysis.scoring_engine_config import TierThresholds

    curated_frame = cr.read_curated(curated)
    negatives_frame = cr.read_negatives(negatives, negatives_query)
    tables = {baseline[0]: load_table(baseline[1]), **{name: load_table(path) for name, path in candidates.items()}}
    samples = {name: build_sample(table, curated_frame, negatives_frame) for name, table in tables.items()}

    targets = {"summary": out / "recalibration_summary.md", "folds": out / "loo_folds.xlsx"}
    existing = [p.name for p in targets.values() if p.exists()]
    if existing and not force:
        raise cr.CalibrationInputError(f"{out} already contains {', '.join(existing)}; use a new --out directory or --force")

    lines = ["# Scoring recalibration analysis", "", "Analysis only; no scoring value was changed by this run.", ""]
    lines += ["## Phase 1: final_score tier thresholds", ""]
    lines += render_phase1(samples[baseline[0]], tables[baseline[0]], TierThresholds(), proposed_thresholds, max_per_query)

    aligned_notes: list[str] = []
    fold_frames: dict[str, pd.DataFrame] = {}
    lines += ["## Phase 2: leave-one-out comparison of candidate settings", ""]
    #: interaction_priority_score is also compared when every table carries it: it is where the coexpression
    #: weights act (they cancel out of interaction_score, see the module notes in the report).
    metrics = ["final_score", "interaction_score"]
    if all(sample.positives["interaction_priority_score"].notna().any() and sample.negatives["interaction_priority_score"].notna().any() for sample in samples.values()):
        metrics.append("interaction_priority_score")
    for metric in metrics:
        aligned, notes = align_samples(samples, metric)
        aligned_notes += notes
        analysis = loo_analysis(aligned, baseline[0], metric)
        lines += [f"### {metric}", ""] + render_loo(analysis, baseline[0])
        fold_frames[metric] = analysis.fold_table
    if aligned_notes:
        lines += ["Alignment notes: " + "; ".join(sorted(set(aligned_notes))), ""]

    out.mkdir(parents=True, exist_ok=True)
    targets["summary"].write_text("\n".join(lines), encoding="utf-8")
    with pd.ExcelWriter(targets["folds"], engine="openpyxl") as writer:
        for metric, frame in fold_frames.items():
            frame.to_excel(writer, sheet_name=metric[:31], index=False)
    return targets


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Tier-threshold and cap/weight recalibration analysis (analysis only).")
    parser.add_argument("--curated", type=Path, required=True)
    parser.add_argument("--negatives", type=Path, required=True)
    parser.add_argument("--baseline", required=True, help="NAME=PATH of the baseline variant's scored rows")
    parser.add_argument("--candidate", action="append", default=[], help="NAME=PATH of a candidate variant (repeatable)")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--negatives-query", default="MA_4115")
    parser.add_argument("--max-per-query", type=int, default=15)
    parser.add_argument("--tier1", type=float, help="proposed tier1_min_score")
    parser.add_argument("--tier2", type=float, help="proposed tier2_min_score")
    parser.add_argument("--tier3", type=float, help="proposed tier3_min_score")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    try:
        from analysis.scoring_engine_config import TierThresholds

        base = _parse_named([args.baseline])
        (baseline_name, baseline_path), = base.items()
        proposed = None
        if any(v is not None for v in (args.tier1, args.tier2, args.tier3)):
            defaults = TierThresholds()
            proposed = TierThresholds(
                tier1_min_score=args.tier1 if args.tier1 is not None else defaults.tier1_min_score,
                tier2_min_score=args.tier2 if args.tier2 is not None else defaults.tier2_min_score,
                tier3_min_score=args.tier3 if args.tier3 is not None else defaults.tier3_min_score,
            )
        targets = run(
            args.curated, args.negatives, (baseline_name, baseline_path), _parse_named(args.candidate), args.out,
            proposed_thresholds=proposed, max_per_query=args.max_per_query, negatives_query=args.negatives_query, force=args.force,
        )
    except cr.CalibrationInputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    for path in targets.values():
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
