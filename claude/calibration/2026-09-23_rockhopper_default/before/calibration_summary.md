# Calibration summary

Generated: 2026-09-23 21:26 by `tools/calibration_report.py`. Analysis only: nothing here changes scoring.

## Inputs

- curated pairs: `claude/experimental_interactions_curated.csv`
- negatives: `claude/calibration/af3_negatives_MA_4115.csv`
- results workbook: `claude/calibration/2026-09-23_rockhopper_default/calib_results_rockhopper_false.xlsx`
- run provenance (`calib_results_rockhopper_false.run_provenance.yaml`):
  - config hash: `96a432d5`
  - code: 5.0 (git 6f6f988)
  - reference genome check: passed
  - generated: 2026-09-23T21:13:05
  - settings: scoring_model=v2_evidence_based; ranking_metric=interaction_priority_score; max_candidates_per_query=2600; geo_coexpression=True; consider_cross_species_matches=True; candidate_sources=candidates,candidates_relaxed,no_hit

## Match coverage

| tier | status | curated rows/directions |
|---|---|---:|
| A | matched | 20 |
| B | candidate_not_in_results | 2 |
| B | matched | 35 |
| excluded | excluded | 2 |

Negatives found in the results: 21 of 28.

## Scores by group (mean / median)

| score | Tier A (n, mean / median) | Tier B (n, mean / median) | Negatives (n, mean / median) |
|---|---|---|---|
| `interaction_score` | 20, 31.42 / 29.52 | 35, 28.02 / 22.13 | 21, 9.70 / 6.14 |
| `interaction_priority_score` | 20, 53.42 / 56.15 | 35, 55.20 / 54.94 | 21, 27.95 / 28.52 |
| `string_neighborhood` | 20, 0.44 / 0.62 | 35, 0.26 / 0.18 | 19, 0.01 / 0.00 |
| `string_cooccurrence` | 20, 0.31 / 0.24 | 35, 0.54 / 0.74 | 19, 0.06 / 0.00 |
| `coexpression_gse77738` | 20, 0.59 / 0.61 | 35, 0.72 / 0.73 | 19, 0.26 / 0.19 |
| `coexpression_gse64349` | 10, 0.83 / 0.80 | 24, 0.35 / 0.21 | 19, 0.36 / 0.35 |
| `final_score` | 20, 43.24 / 43.66 | 35, 41.76 / 38.82 | 21, 26.15 / 25.14 |

## Separation from the negatives

AUC = probability that a random positive scores above a random negative (0.5 = no separation, ties count half); p is the two-sided Mann-Whitney U test.

| score | Tier A vs neg: AUC | p | Tier B vs neg: AUC | p |
|---|---:|---:|---:|---:|
| `interaction_score` | 0.914 | <0.001 (asymptotic) | 0.737 | 0.003 (asymptotic) |
| `interaction_priority_score` | 0.962 | <0.001 (asymptotic) | 0.906 | <0.001 (asymptotic) |
| `string_neighborhood` | 0.913 | <0.001 (asymptotic) | 0.865 | <0.001 (asymptotic) |
| `string_cooccurrence` | 0.732 | 0.007 (asymptotic) | 0.868 | <0.001 (asymptotic) |
| `coexpression_gse77738` | 0.816 | <0.001 (asymptotic) | 0.889 | <0.001 (asymptotic) |
| `coexpression_gse64349` | 0.958 | <0.001 (asymptotic) | 0.373 | 0.159 (asymptotic) |
| `final_score` | 0.876 | <0.001 (asymptotic) | 0.778 | <0.001 (asymptotic) |

Read with care: the samples are small, the negatives come from one query, and the positives from several, so a low p or a high AUC is a prompt to look, not a calibration fit.
