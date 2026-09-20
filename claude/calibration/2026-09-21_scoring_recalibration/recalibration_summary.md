# Scoring recalibration analysis

Analysis only; no scoring value was changed by this run.

## Phase 1: final_score tier thresholds

`final_score` of Tier A (n=11): min 19.7, median 44.3, max 62.5; negatives (n=28): min 9.7, median 23.8, max 49.0. AUC 0.896.

Youden's J is maximal (0.802) for the rule `final_score >= 39.7` (cut midway to the next lower observed score: 37.9).

Sensitivity (share of Tier A at or above t) and false-positive rate (share of negatives at or above t):

| t | Tier A ≥ t | negatives ≥ t | J |
|---:|---:|---:|---:|
| 20 | 0.91 | 0.57 | +0.34 |
| 25 | 0.91 | 0.43 | +0.48 |
| 30 | 0.91 | 0.21 | +0.69 |
| 35 | 0.91 | 0.14 | +0.77 |
| 40 | 0.82 | 0.11 | +0.71 |
| 45 | 0.36 | 0.04 | +0.33 |
| 50 | 0.27 | 0.00 | +0.27 |
| 55 | 0.18 | 0.00 | +0.18 |
| 60 | 0.18 | 0.00 | +0.18 |
| 65 | 0.00 | 0.00 | +0.00 |
| 70 | 0.00 | 0.00 | +0.00 |

Tier of every Tier A pair and negative under the **current** thresholds:

| tier | Tier A | negatives |
|---|---:|---:|
| Tier1_VeryStrong | 0 | 0 |
| Tier2_Strong | 3 | 0 |
| Tier3_Moderate | 7 | 12 |
| Tier4_Weak | 1 | 16 |
| Unclassified | 0 | 0 |

…and under the **proposed** thresholds:

| tier | Tier A | negatives |
|---|---:|---:|
| Tier1_VeryStrong | 0 | 0 |
| Tier2_Strong | 3 | 0 |
| Tier3_Moderate | 7 | 4 |
| Tier4_Weak | 1 | 24 |
| Unclassified | 0 | 0 |

Always-shown safety net (Tier1/Tier2 rows shown regardless of rank, top 15 per query) under the **current** thresholds:

| query | Tier1/2 rows | shown | shown beyond top 15 |
|---|---:|---:|---:|
| MA_0688 | 1 | 15 | 0 |
| MA_3898 | 3 | 15 | 0 |
| MA_3899 | 3 | 15 | 0 |
| MA_3998 | 6 | 15 | 0 |
| MA_4115 | 0 | 15 | 0 |
| MA_4165 | 1 | 15 | 0 |

Always-shown safety net (Tier1/Tier2 rows shown regardless of rank, top 15 per query) under the **proposed** thresholds:

| query | Tier1/2 rows | shown | shown beyond top 15 |
|---|---:|---:|---:|
| MA_0688 | 1 | 15 | 0 |
| MA_3898 | 3 | 15 | 0 |
| MA_3899 | 3 | 15 | 0 |
| MA_3998 | 6 | 15 | 0 |
| MA_4115 | 0 | 15 | 0 |
| MA_4165 | 1 | 15 | 0 |

## Phase 2: leave-one-out comparison of candidate settings

### final_score

Metric: `final_score`. 11 Tier A rows in 10 folds (one fold = one unordered pair; both directions of a pair leave together) vs 28 negatives.

| variant | full AUC | LOO AUC (mean ± jackknife SE) | held-out percentile | ΔAUC vs baseline (SE) | z | folds better / equal / worse | clear improvement |
|---|---:|---:|---:|---:|---:|---:|---|
| baseline | 0.896 | 0.896 ± 0.052 | 0.886 | baseline |  |  |  |
| A | 0.896 | 0.896 ± 0.052 | 0.886 | +0.000 (0.000) | 0.00 | 0 / 10 / 0 | no |
| B | 0.893 | 0.893 ± 0.051 | 0.882 | -0.003 (0.005) | -0.70 | 0 / 9 / 1 | no |
| C | 0.896 | 0.896 ± 0.071 | 0.886 | +0.000 (0.030) | 0.00 | 5 / 3 / 2 | no |
| D | 0.890 | 0.889 ± 0.051 | 0.879 | -0.006 (0.008) | -0.82 | 0 / 9 / 1 | no |

Nested selection (each fold picks its best variant on the other pairs, ties to the baseline): held-out percentile 0.861 vs baseline 0.886; picks per fold: baseline: 8, C: 2.

A variant is called a clear improvement only if ΔAUC ≥ 0.02, z ≥ 2.4 and it is better in at least 75% of the folds where it differs. Nothing is fitted: the variants were fixed in advance.

### interaction_score

Metric: `interaction_score`. 11 Tier A rows in 10 folds (one fold = one unordered pair; both directions of a pair leave together) vs 28 negatives.

| variant | full AUC | LOO AUC (mean ± jackknife SE) | held-out percentile | ΔAUC vs baseline (SE) | z | folds better / equal / worse | clear improvement |
|---|---:|---:|---:|---:|---:|---:|---|
| baseline | 0.906 | 0.906 ± 0.036 | 0.896 | baseline |  |  |  |
| A | 0.906 | 0.906 ± 0.036 | 0.896 | +0.000 (0.000) | 0.00 | 0 / 10 / 0 | no |
| B | 0.906 | 0.906 ± 0.036 | 0.896 | +0.000 (0.000) | 0.00 | 0 / 10 / 0 | no |
| C | 0.896 | 0.896 ± 0.064 | 0.886 | -0.010 (0.038) | -0.26 | 5 / 3 / 2 | no |
| D | 0.899 | 0.899 ± 0.036 | 0.889 | -0.006 (0.006) | -1.02 | 0 / 8 / 2 | no |

Nested selection (each fold picks its best variant on the other pairs, ties to the baseline): held-out percentile 0.864 vs baseline 0.896; picks per fold: baseline: 9, C: 1.

A variant is called a clear improvement only if ΔAUC ≥ 0.02, z ≥ 2.4 and it is better in at least 75% of the folds where it differs. Nothing is fitted: the variants were fixed in advance.

### interaction_priority_score

Metric: `interaction_priority_score`. 11 Tier A rows in 10 folds (one fold = one unordered pair; both directions of a pair leave together) vs 28 negatives.

| variant | full AUC | LOO AUC (mean ± jackknife SE) | held-out percentile | ΔAUC vs baseline (SE) | z | folds better / equal / worse | clear improvement |
|---|---:|---:|---:|---:|---:|---:|---|
| baseline | 0.987 | 0.987 ± 0.013 | 0.986 | baseline |  |  |  |
| A | 0.987 | 0.987 ± 0.013 | 0.986 | +0.000 (0.000) | 0.00 | 0 / 10 / 0 | no |
| B | 0.987 | 0.987 ± 0.013 | 0.986 | +0.000 (0.000) | 0.00 | 0 / 10 / 0 | no |
| C | 0.987 | 0.987 ± 0.013 | 0.986 | +0.000 (0.000) | 0.00 | 0 / 10 / 0 | no |
| D | 0.987 | 0.987 ± 0.013 | 0.986 | +0.000 (0.000) | 0.00 | 0 / 10 / 0 | no |

Nested selection (each fold picks its best variant on the other pairs, ties to the baseline): held-out percentile 0.986 vs baseline 0.986; picks per fold: baseline: 10.

A variant is called a clear improvement only if ΔAUC ≥ 0.02, z ≥ 2.4 and it is better in at least 75% of the folds where it differs. Nothing is fitted: the variants were fixed in advance.
