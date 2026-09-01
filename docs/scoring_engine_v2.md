# Scoring model v2 (evidence-based interaction scoring)

Status: implemented, opt-in, disabled by default (`scoring_model: legacy_additive`).

## Why

The original interaction scorer (`analysis/interaction_scoring.py`, scoring
model `legacy_additive`) adds fixed points for each signal it finds. It has
two accuracy problems, described in the integration design specification:

1. It cannot tell "evaluated this and found nothing" apart from "could not
   evaluate this at all" -- both score `0`. A candidate with no GFF
   coordinates or no annotation text is silently treated the same as one
   that was checked and genuinely has no supporting evidence.
2. `co_occurrence_score` and `domain_complementarity_score` both draw on
   overlapping annotation signals and are added independently, which can
   double count the same underlying evidence.

## What changed

`scoring_model: v2_evidence_based` routes the same raw signals (BLAST
source classification, GFF gene distance, BLAST source overlap, domain and
description keyword complementarity) through:

- `core/evidence.py` -- `EvidenceComponent` / `EvidenceStatus`. Every signal
  carries a status (`AVAILABLE`, `MISSING`, `NOT_APPLICABLE`, ...). Only
  `AVAILABLE` components contribute to a score; everything else is
  excluded from the denominator instead of being counted as zero.
- `analysis/scoring_engine.py` -- groups components into categories
  (`source_classification`, `genomic_context`, `functional_annotation`),
  normalizes each category against only its available evidence, and caps
  each category's contribution so correlated components (`co_occurrence`
  and `domain_complementarity`, both in `functional_annotation`) cannot
  exceed one shared budget. Produces a `ScoreBreakdown` with a 0-100
  `final_score` (or `None` when evidence is insufficient), an `Evidence_Tier`
  (`Tier1_VeryStrong` .. `Tier4_Weak`, or `Unclassified`), and a full,
  auditable per-component trace.
- `analysis/functional_complementarity_rules.py` +
  `config/functional_complementarity_rules.v1.yaml` -- the keyword pair
  table that used to be a Python constant (`COMPLEMENTARY_TERM_PAIRS`) is
  now a versioned, editable YAML ruleset.

The original `legacy_additive` code path (`_score_pair`,
`COMPLEMENTARY_TERM_PAIRS`, etc.) is untouched. Nothing about
`ortholog_filter.py` or the BLAST positive/negative classification changed
in either mode.

## Enabling it

```yaml
interaction_scoring:
  scoring_model: v2_evidence_based
  # optional overrides:
  # scoring_engine_config: "config/scoring_engine.example.yaml"
  # functional_complementarity_ruleset: "config/functional_complementarity_rules.v1.yaml"
```

## New Excel/TSV columns (v2 rows only; blank for legacy rows)

`scoring_model`, `evidence_tier`, `formal_score_available`,
`evidence_category_count`, `evidence_component_count`,
`available_weight_total`.

`interaction_priority_score` is `None`/blank when `formal_score_available`
is `false` (insufficient evidence for a formal score) instead of a
misleading number.

## What is not done yet

- The category caps, per-component weights, and tier thresholds in
  `analysis/scoring_engine_config.py` / `config/scoring_engine.example.yaml`
  reproduce the old point budget as a safe default. They are not
  biologically calibrated. See the design specification, section 37, for
  the calibration plan once known positive/negative pairs are available.
- `alphafold_readiness_score` is still computed and shown for reference in
  both modes, but it is intentionally excluded from the v2 total score
  (pair length is not interaction evidence).
- Integrating ProteinInteractionHunter's own evidence types (orthology,
  phylogenetic profile, known-interaction databases) is out of scope for
  this change; see the integration design specification, section 1, for
  why the two projects stay code-independent.
