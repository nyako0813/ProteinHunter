# Calibration summary

Generated: 2026-09-21 00:01 by `tools/calibration_report.py`. Analysis only: nothing here changes scoring.

## Inputs

- curated pairs: `claude/experimental_interactions_curated.csv`
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
| A | matched | 11 |
| B | candidate_not_in_results | 1 |
| B | matched | 4 |
| B | no_query_in_results | 14 |
| excluded | excluded | 2 |

Negatives found in the results: 21 of 28.

## Scores by group (mean / median)

| score | Tier A (n, mean / median) | Tier B (n, mean / median) | Negatives (n, mean / median) |
|---|---|---|---|
| `interaction_score` | 11, 33.65 / 29.94 | 4, 2.76 / 0.00 | 21, 9.70 / 6.14 |
| `interaction_priority_score` | 11, 56.20 / 57.51 | 4, 34.07 / 37.73 | 21, 27.95 / 28.52 |
| `string_neighborhood` | 11, 0.47 / 0.65 | 4, 0.02 / 0.00 | 19, 0.01 / 0.00 |
| `string_cooccurrence` | 11, 0.31 / 0.26 | 4, 0.06 / 0.00 | 19, 0.06 / 0.00 |
| `coexpression_gse77738` | 11, 0.57 / 0.47 | 4, 0.44 / 0.46 | 19, 0.26 / 0.19 |
| `coexpression_gse64349` | 6, 0.82 / 0.79 | 1, 0.31 / 0.31 | 19, 0.36 / 0.35 |
| `final_score` | 11, 46.13 / 44.29 | 4, 23.18 / 23.33 | 21, 26.15 / 25.14 |

## Separation from the negatives

AUC = probability that a random positive scores above a random negative (0.5 = no separation, ties count half); p is the two-sided Mann-Whitney U test.

| score | Tier A vs neg: AUC | p | Tier B vs neg: AUC | p |
|---|---:|---:|---:|---:|
| `interaction_score` | 0.922 | <0.001 (asymptotic) | 0.179 | 0.049 (asymptotic) |
| `interaction_priority_score` | 0.983 | <0.001 (asymptotic) | 0.702 | 0.231 (exact) |
| `string_neighborhood` | 0.921 | <0.001 (asymptotic) | 0.559 | 0.580 (asymptotic) |
| `string_cooccurrence` | 0.751 | 0.011 (asymptotic) | 0.513 | 0.958 (asymptotic) |
| `coexpression_gse77738` | 0.823 | 0.004 (asymptotic) | 0.684 | 0.286 (exact) |
| `coexpression_gse64349` | 0.956 | 0.001 (asymptotic) | 0.316 | 0.700 (exact) |
| `final_score` | 0.896 | <0.001 (asymptotic) | 0.393 | 0.529 (asymptotic) |

Read with care: the samples are small, the negatives come from one query, and the positives from several, so a low p or a high AUC is a prompt to look, not a calibration fit.
