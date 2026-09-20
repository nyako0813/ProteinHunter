# Calibration summary

Generated: 2026-09-21 00:00 by `tools/calibration_report.py`. Analysis only: nothing here changes scoring.

## Inputs

- curated pairs: `tests/fixtures/curated_pairs_calibration_report.csv`
- negatives: `claude/calibration/af3_negatives_MA_4115.csv`
- results workbook: `data/output/calibration_a2_mcr_2026-09-20.xlsx`
- run provenance (`calibration_a2_mcr_2026-09-20.run_provenance.yaml`):
  - config hash: `2881e05d`
  - code: 5.0 (git 397dacf+dirty)
  - reference genome check: passed
  - generated: 2026-09-20T23:57:09
  - settings: scoring_model=v2_evidence_based; ranking_metric=interaction_priority_score; max_candidates_per_query=2600; geo_coexpression=True; consider_cross_species_matches=True; candidate_sources=candidates,candidates_relaxed,no_hit

## Match coverage

| tier | status | curated rows/directions |
|---|---|---:|
| A | matched | 8 |
| B | candidate_not_in_results | 1 |
| B | matched | 5 |
| B | no_query_in_results | 14 |
| excluded | excluded | 2 |

Negatives found in the results: 21 of 28.

## Scores by group (mean / median)

| score | Tier A (n, mean / median) | Tier B (n, mean / median) | Negatives (n, mean / median) |
|---|---|---|---|
| `interaction_score` | 8, 35.15 / 33.32 | 5, 8.20 / 0.00 | 21, 9.70 / 6.14 |
| `interaction_priority_score` | 8, 56.00 / 58.40 | 5, 38.83 / 39.96 | 21, 27.95 / 28.52 |
| `string_neighborhood` | 8, 0.62 / 0.72 | 5, 0.02 / 0.00 | 19, 0.01 / 0.00 |
| `string_cooccurrence` | 8, 0.15 / 0.11 | 5, 0.19 / 0.00 | 19, 0.06 / 0.00 |
| `coexpression_gse77738` | 8, 0.52 / 0.39 | 5, 0.50 / 0.61 | 19, 0.26 / 0.19 |
| `coexpression_gse64349` | 3, 0.85 / 0.77 | 2, 0.56 / 0.56 | 19, 0.36 / 0.35 |
| `final_score` | 8, 46.89 / 46.66 | 5, 27.40 / 23.33 | 21, 26.15 / 25.14 |

## Separation from the negatives

AUC = probability that a random positive scores above a random negative (0.5 = no separation, ties count half); p is the two-sided Mann-Whitney U test.

| score | Tier A vs neg: AUC | p | Tier B vs neg: AUC | p |
|---|---:|---:|---:|---:|
| `interaction_score` | 0.917 | <0.001 (asymptotic) | 0.333 | 0.268 (asymptotic) |
| `interaction_priority_score` | 0.976 | <0.001 (exact) | 0.762 | 0.079 (exact) |
| `string_neighborhood` | 0.931 | <0.001 (asymptotic) | 0.626 | 0.208 (asymptotic) |
| `string_cooccurrence` | 0.658 | 0.137 (asymptotic) | 0.611 | 0.376 (asymptotic) |
| `coexpression_gse77738` | 0.796 | 0.018 (asymptotic) | 0.726 | 0.139 (exact) |
| `coexpression_gse64349` | 0.965 | 0.013 (asymptotic) | 0.632 | 0.610 (exact) |
| `final_score` | 0.881 | 0.002 (asymptotic) | 0.505 | 1.000 (asymptotic) |

Read with care: the samples are small, the negatives come from one query, and the positives from several, so a low p or a high AUC is a prompt to look, not a calibration fit.
