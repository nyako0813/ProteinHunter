# Repository analysis & implementation plan: interaction_scoring v2

Written against `ProteinHunter_v5 × ProteinInteractionHunter 統合設計書 v1.0`
(the ChatGPT design specification), section 56 ("開発開始時の最初の要求").
Audited commit: `9a26dbe` (ProteinHunter_v5). This document reports what
exists, what changed, and what is deliberately left for a later pass.

## 1. Repository analysis (pre-change)

| Area | File(s) | Finding |
|---|---|---|
| Evidence extraction | `analysis/interaction_scoring.py` | Raw signals (BLAST source, GFF distance, source-set Jaccard, keyword match) computed and immediately converted to fixed points in one pass; no separate raw/normalized/status representation. |
| Score calculation | same, `_score_pair` | Plain sum of five components, weights configurable but caps/missing-handling absent. `0` used both for "evaluated, no support" and "could not evaluate". |
| Config | `config.py`, `InteractionScoringConfig` / `InteractionScoringWeightsConfig` | Hand-rolled dataclasses + manual YAML validation (~1000 lines total); no category concept, no evidence-status concept. |
| Output | `output/excel.py`, `INTERACTION_PAIR_COLUMNS` | Row dicts turned into a `pandas.DataFrame` via an explicit column list; missing dict keys become `NaN` automatically, which made additive (non-breaking) column growth possible. |
| Tests | `tests/test_interaction_scoring.py` (598 lines), `tests/test_scoring.py`, `tests/test_ortholog_filter.py` | 258 tests passing at audit time (baseline recorded before any change). |
| homolog removal | `analysis/ortholog_filter.py`, `ortholog_filter` config section | BLAST negative-hit strength classification (`strong/medium/weak/none`) + exclusion mode. Independent of interaction scoring; explicitly out of scope for this change per user instruction. |
| ProteinInteractionHunter | separate repository | MVP-1 evidence-based design (categories, caps, `EvidenceEvent.status`) used as a *reference model only*; no code imported (see the design specification's own stated boundary, and this document's companion `ProteinHunter_v5_scoring_v2_design.md`). |

## 2. What changed

New, additive-only files (nothing existing was deleted):

- `core/evidence.py` -- `EvidenceStatus`, `EvidenceComponent` (Phase 1).
- `analysis/scoring_engine.py` -- `score_candidate`, `rank_candidates`, `ScoreBreakdown`, `CategoryScore` (Phase 2).
- `analysis/scoring_engine_config.py` -- `ScoringEngineConfig` + YAML loader (Phase 2).
- `analysis/functional_complementarity_rules.py` + `config/functional_complementarity_rules.v1.yaml` -- externalized keyword-pair ruleset (Phase 3).
- `config/scoring_engine.example.yaml` -- example/default engine tuning file.
- `docs/scoring_engine_v2.md` -- user-facing explanation of the new model.
- `tests/test_evidence.py`, `tests/test_scoring_engine.py`, `tests/test_scoring_engine_config.py`, `tests/test_functional_complementarity_rules.py` -- new unit tests.

Modified files:

- `analysis/interaction_scoring.py` -- added a parallel `v2_evidence_based` code path (`_score_pair_v2`, `_rank_source_candidates_v2`, `_build_evidence_components_v2`, and status-aware counterparts of the existing raw-signal helpers). The original `legacy_additive` functions (`_score_pair`, `_rank_source_candidates`, `COMPLEMENTARY_TERM_PAIRS`, `MEANINGFUL_KEYWORDS`, `DESCRIPTION_STOPWORDS`) are untouched. `INTERACTION_PAIR_COLUMNS` gained six new trailing columns (`scoring_model`, `evidence_tier`, `formal_score_available`, `evidence_category_count`, `evidence_component_count`, `available_weight_total`) that stay empty for legacy rows.
- `config.py` -- `InteractionScoringConfig` gained three new fields with defaults (`scoring_model="legacy_additive"`, `scoring_engine_config=None`, `functional_complementarity_ruleset=None`), plus matching validation/loading. Every existing constructor call and every existing test that builds this dataclass keeps working unmodified because the new fields are optional.
- `config.yaml` -- documented (but not enabled) the new `scoring_model` key.
- `tests/test_interaction_scoring.py` -- `interaction_config()` test helper extended with the same three optional keyword arguments; existing calls unaffected. Added 7 new integration tests exercising `run_interaction_scoring(..., scoring_model="v2_evidence_based")` end to end.
- `CHANGELOG.md` -- entry for this change.

`analysis/ortholog_filter.py` and everything under `ortholog_filter:` in config: **not touched**.

## 3. Risk assessment

| Risk | Mitigation |
|---|---|
| Breaking the 258 existing tests | Every new field has a backward-compatible default; legacy code paths are untouched, not refactored in place. Full suite re-run after every phase (see below). |
| Config surface growing too large in `config.py` | New tunables (caps, penalties, tier thresholds) live in a separate `scoring_engine_config.py` + standalone YAML, not inside the already-1000-line `config.py` validator. |
| Silent double counting | Enforced structurally: `co_occurrence` and `domain_complementarity` share one category (`functional_annotation`) with one cap; `analysis/scoring_engine.py::_score_categories` cannot let a category exceed its configured cap regardless of how many components fire. |
| Silent zeroing of missing evidence | Enforced structurally: `EvidenceComponent.unavailable(...)` always has `effective_weight == 0` and `contribution == 0`, and is excluded from `available_weight`/`total_cap` in `score_candidate`. Unit-tested directly (`test_missing_evidence_is_excluded_from_denominator_not_zeroed`). |
| Uncalibrated weights presented as final | `docs/scoring_engine_v2.md` and the YAML files' own comments state explicitly that the numbers are a placeholder migration of the old point budget, not a calibrated model. |

## 4. Test results

```text
Baseline (before any change): 258 passed
After Phase 1 (evidence model):            + 17 tests, 0 regressions
After Phase 2 (scoring engine + config):    + 19 tests, 0 regressions
After Phase 3 (wiring + ruleset + v2 tests): + 16 tests, 0 regressions
Final: 308 passed, 0 failed
```

## 5. Recommended order (as executed)

1. Baseline test run.
2. Evidence model (`core/evidence.py`) -- pure, dependency-free, easiest to get right first.
3. Scoring engine (`analysis/scoring_engine.py` + `scoring_engine_config.py`) -- depends only on the evidence model, testable with synthetic components, no BLAST/GFF/Excel involved.
4. Wire v5's real evidence into the engine as an opt-in path, externalize the keyword ruleset.
5. Regression run, commit on a feature branch.

## 6. Explicitly out of scope for this change (tracked for a later pass)

Per the design specification's own phase breakdown (sections 47, Phase 4/6/7/8/9/10):

- Phase 4 -- importing ProteinInteractionHunter's own evidence types (orthology, phylogenetic profile, known-interaction databases). The two projects remain code-independent; see `ProteinHunter_v5_scoring_v2_design.md` section 1 for the reasoning.
- Phase 5 (partial) -- a true "Unified Score" combining a separate ProteinHunter score and a separate InteractionHunter score; only ProteinHunter's own evidence is modeled here.
- Phase 6/7/8 -- the redesigned multi-sheet Excel index/final-score-first layout and the single-file Word report with table of contents and Excel cross-links. The current change only adds columns to the existing `Interaction_*` sheets.
- Phase 9/10 -- full real-data end-to-end run and an old-vs-new score/rank comparison report.
- Weight/cap calibration against known positive/negative interaction pairs (design specification section 37) -- no such labeled data was available at the time of this change.
