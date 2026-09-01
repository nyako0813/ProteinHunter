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
- `ortholog_filter.py`'s existing negative-hit classification
  (`record.negative_hit_strength`, unchanged) is read as a `negative_hit_strength`
  evidence component (category `source_reliability`, `is_negative=True`).
  It contributes nothing when there is no negative BLAST hit
  (`NOT_APPLICABLE`, not a phantom `0`), and otherwise subtracts a capped
  penalty from the final score instead of the legacy behavior of hard
  filtering. `ortholog_filter.py` itself is not modified; only its
  already-computed output is reused.
- `analysis/pih_evidence_bridge.py` (see "ProteinInteractionHunter bridge"
  below) optionally folds in evidence from a separate tool's own output
  file.

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
  # pih_evidence_bundle: "path/to/candidate_evidence_bundle.jsonl"
```

## ProteinInteractionHunter (PIH) bridge (optional)

ProteinInteractionHunter is a separate, independent project; its own design
specification states it must never import, or be imported by, this
project. `analysis/pih_evidence_bridge.py` respects that boundary: it never
imports PIH's package, and only reads PIH's own machine-readable output
file (`candidate_evidence_bundle.jsonl`, produced by
`protein-interaction-hunter generate-candidates`) as plain JSON lines --
the same way any other external tool's report could be read.

PIH's own "integrated scoring" groups evidence into five categories:
`genomic_context`, `functional_annotation`, `cellular_compatibility`,
`evolutionary`, and `direct_interaction`. This project already computes its
own `genomic_context` (GFF gene distance) and `functional_annotation`
(BLAST source overlap + domain/description keywords) independently, so
folding PIH's versions of those same two categories in would double count
the same kind of signal computed two different ways. The bridge therefore
only imports the three categories this project has no equivalent for:
`cellular_compatibility` (localization/topology), `evolutionary`
(orthology / phylogenetic profile), and `direct_interaction` (gene fusion /
known interaction databases). Each is added as its own evidence component,
prefixed `pih_` (`pih_cellular_compatibility`, `pih_evolutionary`,
`pih_direct_interaction`) so it is always visually and structurally
distinct from this project's own categories of a conceptually similar
type.

Set `interaction_scoring.pih_evidence_bundle` to the path of PIH's
`candidate_evidence_bundle.jsonl` to enable it. Leaving it unset runs
without any PIH evidence at all -- identical to before this bridge
existed. A configured but missing or malformed bundle file degrades
gracefully (each bad line, or the whole missing file, becomes a run
warning) rather than aborting the run, since this evidence source is
optional and best-effort. Because PIH and this project were not designed
to share an identifier convention, the bridge tries several ID spellings
on both sides (`protein_id`, `old_locus_tag`, and a version-stripped
`protein_id`) when matching a query/candidate pair to PIH's records.

Category caps for the bridged categories live in
`analysis/scoring_engine_config.py::DEFAULT_CATEGORY_CAPS` /
`config/scoring_engine.example.yaml` (`pih_cellular_compatibility: 5`,
`pih_evolutionary: 10`, `pih_direct_interaction: 20`) and must match
`analysis.pih_evidence_bridge.BRIDGED_PIH_CATEGORY_CAPS` --
`tests/test_scoring_engine_config.py::test_example_config_matches_defaults`
guards this so the example file cannot silently drift out of sync with the
code again (a category with no configured cap raises a `ConfigError` the
first time it actually fires for a real pair, not at load time, so a
forgotten cap would otherwise fail silently until real bridge data was
used).

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
- The PIH bridge's per-category weights (`PIH_CATEGORY_WEIGHTS` /
  `pih_*` caps) are provisional, matching the same "reproduce the old point
  budget" placeholder philosophy as the rest of v2 -- not yet calibrated
  against known interacting/non-interacting pairs.
- The PIH bridge is read-only, file-based, and one-directional (v5 reads
  PIH's output; PIH is never given v5's output, and neither project
  imports the other's code). Running PIH itself, keeping its output fresh,
  and resolving identifier mismatches beyond the handful of ID spellings
  this bridge already tries remain the user's / a future orchestration
  layer's responsibility.
- Phases 5 and later of the integration design specification (unified
  ranking refinements, Excel/Word report redesign, PIH-vs-v5 comparison
  reports, end-to-end runs on real data, weight calibration) are not yet
  started.
