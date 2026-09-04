# Phase 6/7/8 (design spec numbering) Stage 2: Word report generation — investigation and design proposal

Status: **investigation and design proposal only, nothing implemented**, per
the request that started this document. Continues directly from
`claude/phase678_excel_word_redesign_investigation.md` (Stage 1, merged as
PR #13) and follows the same process used for that document and for
`claude/phase6_external_evidence_design.md` / `claude/phase6b_coexpression_design.md`.

**Terminology / provenance note, carried over from Stage 1's document**:
this session still does not have the literal design-spec source text (see
Stage 1 doc, item 1 — confirmed absent from this repo and from every local
path searched at the time). What follows treats the section numbers and
paraphrased rules the request quotes directly (design spec §29 "one Word
file per run", §34 "Why this candidate ranks highly" / "Biological
Interpretation", §35 "do not assert this protein IS the target enzyme,
say it ranked highly given currently available evidence", §39 "do not
carelessly remove existing functionality / show limitations honestly",
§45 "reproducibility: same input + config + version ⇒ same result") as
authoritative input, the same way Stage 1 treated the user's own Stage 1
directive as authoritative once the literal spec text was confirmed
unavailable. Everything else below is grounded in what is actually in the
repository today (code, config, CHANGELOG, and the Stage 1 real-data
verification), cited by file:line.

Confirmed already decided (carried over from the request, not
re-litigated here): `python-docx` (MIT), the A+C+E cross-link combination
(Excel→Word external hyperlink + Word bookmarks/Excel links + a static
summary table embedded in Word), "8. Candidate Details" as a top-candidate
summary rather than a full transcription, ranking by `final_score`, and
support for both `scoring_model` values.

---

## 1. Multiple queries in a single Word file

**No implementation obstacle found.** The data this report would consume
is already query-partitioned and already sorted by query first:

- `FINAL_SCORE_COLUMNS` and `SCORE_BREAKDOWN_COLUMNS` both carry `query_id`
  / `query_protein_id` as leading columns (`output/excel.py:86-107`,
  `output/excel.py:129-158`).
- `report_v2.py::rerank_final_score_rows` already sorts consolidated rows
  with `query_id` as the primary sort key before `candidate_rank`
  (`output/report_v2.py:149-172`, key tuple starts
  `str(row["query_id"])`).

So a Word-generation module can consume the exact same consolidated rows
`build_workbook_sheets()` already produces (`output/report_v2.py:290`),
group them by `query_id` (a plain `itertools.groupby` over already-sorted
rows, no new sorting/joining logic needed), and emit one subsection per
query under both "7. Candidate Ranking" and "8. Candidate Details" —
matching the request's own reading of §29 ("one file per pipeline run",
not "one file per query"). Recommended heading scheme: `7.1`, `7.2`, ...
one per query in "7. Candidate Ranking", and `8.1`, `8.2`, ... mirroring
the same query order in "8. Candidate Details" (Word's native heading
levels drive its built-in Table of Contents field for free — see the M1
note below on why that field still needs one manual "update field" step
in Word itself).

Three things worth flagging before implementation, none of them blockers:

1. **Query order must be deterministic and config-driven, not
   incidentally stable.** §45 requires the same input/config/version to
   reproduce the same result — that includes section *order*, not just
   section *content*. `rerank_final_score_rows`'s sort key is
   `query_id` (string), so today's order is alphabetical by query ID, not
   necessarily the order queries were listed in
   `interaction_scoring.query_proteins` (`config.yaml:144`). That is
   already deterministic (good for §45) but may read oddly to a
   biologist opening the report (e.g. `MA_0688` before `MA_4115` even if
   the config listed `MA_4115` first). Worth a one-line decision at M4
   time: keep alphabetical-by-`query_id` (simplest, already deterministic,
   zero new code) or thread config order through
   (`build_workbook_sheets` would need to carry it) — not a design
   blocker either way, just needs a conscious pick rather than depending
   on whatever a dict happens to iterate in.
2. **Scale.** The Stage 1 real-data run produced 4,627 classified records
   / 5 queries / **23,130** consolidated `02_Final_Score` rows (Stage 1
   doc, "Real-data verification" section). A Word report that mirrored
   every row, for every query, would be enormous and unreadable as a
   single file — this is exactly why item 2 below (Top-N per query) is
   load-bearing for item 1, not a separate concern: "one file, sectioned
   by query" is only tractable if each query's section is bounded.
3. **Word's Table of Contents field is not self-populating from
   `python-docx`.** `python-docx` can insert the TOC field code
   (`w:fldSimple`/`TOC \o "1-3" \h \z \u`) but cannot compute page
   numbers itself (well-documented library limitation, not specific to
   this project) — Word shows "Right-click to update field" until the
   user does so once, or the document is saved with
   `updateFields`/`w:settings` set so Word auto-updates on open. Confirm
   which behavior is acceptable during the M1 smoke test (see milestones).

---

## 2. How many top candidates per query

**The exact Tier1–4 / candidate-count distribution from the PR #13
real-data run is not available to check.** Per the Stage 1 document, that
run's diagnostic outputs lived under `.cache/geo_investigation/`
(gitignored, never committed) and its comparison `git worktree` was
removed after use — both confirmed still true: `.cache/geo_investigation`
does not exist, `data/output/` and `data/local_backup/` are empty, and no
`.xlsx` artifact from that run survives anywhere on this machine. This
document does not fabricate a distribution it cannot show; recommendation
below is built from what *is* independently known and is real-data-derived,
not guessed, plus an explicit offer to close the remaining gap.

**What is known, and it is more informative than a raw tier histogram
would be:**

- Tier thresholds (`config/scoring_engine.example.yaml:62-68`, this is the
  *default* — nothing in the repo overrides it):

  | Tier | Score floor | Category-count floor |
  |---|---|---|
  | Tier1_VeryStrong | ≥ 70 | ≥ 3 |
  | Tier2_Strong | ≥ 50 | ≥ 2 |
  | Tier3_Moderate | ≥ 25 | ≥ 1 |
  | Tier4_Weak | (below Tier3) | — |
  | Unclassified | no formal score at all | — |

  The same file's own header comment: *"None of these numbers are
  biologically calibrated yet -- they reproduce the point budget of the
  old additive scorer as a safe starting point... recalibrate them (see
  the design specification, section 37) rather than guessing"*
  (`config/scoring_engine.example.yaml:9-12`). This is the codebase's own
  admission that tier boundaries are provisional, not a settled ground
  truth — directly relevant to whether a tier cutoff should gate Word
  inclusion (see conclusion below).

- The Stage 1 real-data verification's matched-config comparison (36
  pairs: 8 curated Tier A true-positive pairs + 28 AlphaFold3-confirmed
  true-negative pairs, same STRING+GEO-enabled config as the Final Score
  integration phase) measured actual `final_score` values, not just
  pass/fail:  **POS mean 42.884, NEG mean 25.680** (`CHANGELOG.md:73`,
  reproducing the Final Score integration phase's own +17.20 separation
  exactly).

**Reading these two facts together is the actual finding of this
section:**

- The known-true-**positive** mean (42.884) sits **below the Tier2 floor
  (50)**, let alone Tier1 (70) — under current, self-admittedly
  uncalibrated weights, a real interacting pair lands in Tier3/Tier4
  territory *on average*, not Tier1/Tier2.
- The known-true-**negative** mean (25.680) sits **almost exactly at the
  Tier3 floor (25)** — i.e. the tier boundary that would separate "shown
  in the report" from "not shown" under a naive Tier1–3 cutoff sits right
  on top of where confirmed non-interacting pairs already average.

Two consequences follow directly, and both argue against tier-based
inclusion as the primary mechanism:

1. **A "Tier1–2 only" cutoff (the request's first example) would likely
   under-fill or empty out sections for real queries** — the only
   real-data mean available for genuine positives doesn't clear Tier2.
2. **Any tier-based cutoff, including Tier1–3, does not yet cleanly
   separate real candidates from real non-candidates** — the ~17-point
   POS/NEG mean gap straddles the Tier3 boundary itself, not a clean
   margin above it.

**Recommendation: rank-based Top-N per query, config-driven, not a tier
cutoff** — the same style the codebase already uses for
`max_candidates_per_query` (`config.yaml:187`, an existing
config-not-hardcoded knob for a structurally similar problem). Show
`evidence_tier`/`final_score_tier` as a visible *label* on every included
candidate (this is what makes the tier's honesty value land — see item 4
— without using it as a silent gate that could empty a section or hide
the fact that inclusion and confidence are different axes). A concrete
starting point:

- Default **Top 15 per query**, exposed as a new config key (naming to
  match the existing `paths.output_excel` / `output_word` pattern
  discussed in item 5, e.g.
  `interaction_scoring.word_report.max_candidates_per_query: 15`) —
  large enough that genuine positives (which the real-data check shows
  can rank outside a very tight Top 5 given current score compression)
  are not routinely cut, small enough that 5 queries × 15 stays a
  legible single file (75 write-ups vs. 23,130 raw pairs).
- A safety-net rule on top of the count, **not instead of it**: always
  include any candidate that reaches Tier1_VeryStrong or Tier2_Strong for
  that query even if its rank falls past N — the count is a *floor* for
  "show enough to be useful," the tier is an *uncapped inclusion
  guarantee* for "never silently drop a candidate the pipeline itself
  is confident about," never the reverse (tier is not used to shrink the
  N candidates below the configured count).

This is offered as a reasoned starting default, not a final number — it
is explicitly an **open question** below, and this session can close the
data gap directly: the same `git worktree` + matched-config pattern Stage
1 used for its real-data verification could be re-run to produce an
actual Tier1–4 count histogram across all 23,130 rows (not just the 36
calibration pairs) before committing to a specific N, if that precision
is wanted before Stage 2 implementation starts. Not run automatically
here because it touches external APIs (STRING/GEO/CDD/Pfam/UniProt/
AlphaFold) and, per Stage 1's own account, is not fast — better to
confirm you want it than to spend that time unasked.

---

## 3. "Why this candidate ranks highly" / "Biological Interpretation" — deterministic template design

**Constraint taken as given (§45):** no LLM call, no network call, no
randomness — pure function of the row's already-computed fields, so the
same input always produces the same string, and the logic is unit
-testable the same way every other module in this codebase is (fixed
input → fixed expected string, no mocking of an external service needed
because there is no external service).

**Fields available to build from** (already present on every
`04_Score_Breakdown` row — see `SCORE_BREAKDOWN_COLUMNS`,
`output/excel.py:129-158`): `final_score`, `final_score_tier`,
`evidence_category_count`, `evidence_component_count`, `candidate_source`,
`negative_hit_strength`, and the per-category reference columns
`candidate_priority_score` (source_classification proxy),
`same_gene_neighborhood_score` (genomic_context proxy),
`functional_domain_score`, `evolutionary_score`,
`cellular_compatibility_score`, `interaction_evidence_score`. Each
reference column is `None`/blank when its category had zero available
evidence for this pair — *"blank means not evaluated, not scored zero"*
is already the established convention for every one of these columns
(`INTERACTION_SCORE_EXPLANATIONS`, `output/excel.py:294-462` — e.g. the
`evolutionary_score` entry: *"blank (not a scored zero) when no PIH
evolutionary evidence was evaluated for this pair"*). The template logic
below only has to branch on `is None`, reusing that same convention
rather than inventing a new one.

### 3a. "Why this candidate ranks highly"

Structure: one opening (tier) sentence, one evidence-enumeration
sentence built by looping only over non-`None` category columns, one
`candidate_source` sentence, and an optional `negative_hit_strength`
caveat sentence — four deterministic template pieces, concatenated.

**Opening sentence, branch on `final_score_tier`:**

```text
Tier1_VeryStrong:
  "This candidate reached the highest confidence tier for this query
  (Tier 1 — Very Strong), with a Final Score of {final_score:.1f}/100
  supported by evidence from {evidence_category_count} independent
  evidence categories."

Tier2_Strong:
  "This candidate reached the second-highest confidence tier for this
  query (Tier 2 — Strong), with a Final Score of {final_score:.1f}/100
  supported by evidence from {evidence_category_count} independent
  evidence categories."

Tier3_Moderate:
  "This candidate reached a moderate confidence tier for this query
  (Tier 3 — Moderate), with a Final Score of {final_score:.1f}/100.
  The supporting evidence is present but comes from fewer independent
  categories ({evidence_category_count}) than higher-tier candidates."

Tier4_Weak:
  "This candidate's Final Score of {final_score:.1f}/100 places it in
  the weakest confidence tier evaluated (Tier 4 — Weak). Its position in
  this list reflects a relative rank among the candidates scored for
  this query, not strong independent support."

Unclassified:
  "This candidate did not receive a formal Final Score for this query
  (insufficient evidence to meet the pipeline's minimum evidence
  requirement). It appears here only because it was among the
  higher-ranked candidates by whatever partial evidence was available;
  treat its inclusion as informational, not as a ranked, scored result."
```

**Evidence-enumeration sentence** (loop over available — i.e. non-`None`
— category columns; cap values read from the run's actual
`scoring_engine_config`, never hardcoded into the template string):

```text
"Contributing evidence categories: {", ".join(
    f"{label} ({value:.1f}/{cap:.0f})"
    for label, value, cap in available_categories
)}."

# label mapping (fixed, matches item 4's category names):
#   candidate_priority_score        -> "Sequence/Source Classification"
#   same_gene_neighborhood_score    -> "Genomic Context"
#   functional_domain_score         -> "Functional/Domain"
#   evolutionary_score              -> "Evolutionary"
#   cellular_compatibility_score    -> "Cellular Compatibility"
#   interaction_evidence_score      -> "Interaction"
```

**`candidate_source` sentence** (fixed text per value, six branches —
wording adapted from the existing `candidate_source` glossary entry,
`output/excel.py:429-435`, so the Word report and the Excel workbook
describe the same bucket the same way):

```text
Candidates:
  "This protein was classified as a strict positive candidate — a BLAST
  hit to a positive reference sequence with no negative-reference hit at
  all, the most stringent classification this pipeline produces."

Positive_all_sources:
  "This protein hit every configured positive reference source, with no
  negative-reference hit."

Candidates_relaxed:
  "This protein was classified as a positive candidate under relaxed
  criteria — it has a positive-reference hit, and tolerates a
  non-strong (medium or weak) negative-reference hit."

No_hit:
  "This protein had no BLAST hit to either the positive or negative
  reference sets. It is a lineage-specific or novel candidate not
  detected by sequence homology alone — its ranking here depends
  entirely on non-sequence evidence (genomic context, domain
  annotation, interaction evidence), not on a BLAST match."

Negative_unmatched:
  "This protein has no negative-reference hit, though it did not meet
  the stricter positive-match criteria above."

Negative_hit:
  "This protein has at least one negative-reference BLAST hit
  (strength: {negative_hit_strength}). It is included in this ranking
  despite that hit because of the other evidence enumerated above;
  interpret its position with additional caution."
```

**`negative_hit_strength` caveat** (only appended when the value is not
`"none"`, regardless of `candidate_source`, since a weak/medium hit can
occur outside the `Negative_hit` bucket too):

```text
if negative_hit_strength != "none" and candidate_source != "Negative_hit":
  "Note: this candidate also has a {negative_hit_strength} BLAST hit to
  a negative-reference sequence, which does not change its
  candidate_source classification but is a relevant caveat."
```

### 3b. "Biological Interpretation" — the §35-critical section

§35, as quoted in the request, is the binding constraint here: never
assert the candidate *is* the target enzyme/interaction partner; always
frame the statement as "ranked highly given currently available
evidence." The template enforces this by making the hedge the mandatory
opening clause (never optional, never rephrased away), then only adding
category-specific color for categories that actually fired, then closing
with an explicit, honest statement of what was **not** evaluated —
directly implementing §39's "show limitations honestly" for the one
place in the report a reader is most likely to over-read a claim.

```text
Opening (always present, every tier, no branch skips this):
  "Based on the evidence currently available to this pipeline,
  {candidate_id} ranked {rank_ordinal} of {n_candidates_for_query}
  evaluated candidates for query {query_id}. This reflects the evidence
  categories described below as currently available to the pipeline —
  it is not a confirmed identification of {candidate_id} as the target
  enzyme or interaction partner, and should be read as 'best-supported
  by currently available evidence,' not as a settled conclusion."

Category color clauses (only for non-None categories, one sentence
each, deliberately hedged — "consistent with", never "confirms" /
"proves" / "shows that X is Y"):

  same_gene_neighborhood_score (if present and > 0):
    "Its genomic proximity to the query gene is consistent with — but
    does not by itself establish — a shared operon or functional
    module."

  functional_domain_score (if present and > 0):
    "Shared or complementary functional domain annotations were found
    between this candidate and the query, consistent with a related or
    complementary biochemical role; domain similarity alone does not
    establish direct interaction."

  interaction_evidence_score (if present and > 0):
    "External interaction evidence (protein-protein interaction
    database and/or coexpression data) supports a functional
    association between this candidate and the query, independent of
    sequence similarity or genomic position."

  evolutionary_score (if present and > 0):
    "Evolutionary/phylogenetic profile evidence from the optional PIH
    bridge supports consistency between this candidate and the query
    across the genomes compared."

  cellular_compatibility_score (if present and > 0):
    "Cellular-compatibility evidence from the optional PIH bridge
    supports this candidate and the query being compatible with the
    same subcellular context."

Honesty closer (always present — this is the §39 tie-in):
  "Categories not listed above were not evaluated for this candidate in
  this run — a blank category means no evidence was available to
  assess, not that the evidence was checked and found absent. In
  particular: Evolutionary and Cellular Compatibility evidence require
  an optional external data bundle that {"was" if pih_bundle_configured
  else "was not"} supplied for this run, and Negative (biological
  contradiction) evidence is a reserved category with no implemented
  signal in this pipeline version (see Evidence Architecture, section
  5.7) — its absence here does not mean this candidate was checked for
  contradicting evidence and cleared."
```

The closer's PIH-bundle branch and the permanent Negative-evidence
sentence are not candidate-specific — they are the same two sentences on
every candidate in a given run (first depends only on whether
`pih_evidence_bundle` was configured for the run at all; second is
always true today) — so in practice this closer can be generated once
per run and reused, rather than recomputed per candidate; noted here so
the eventual implementation doesn't do 75× redundant string-building for
no reason.

### 3c. Where this logic should live

Recommend a new, pure module — e.g. `output/word_narrative.py` — that
takes a row `dict` (already shaped like a `04_Score_Breakdown` row) and a
small run-level context (which `scoring_engine_config` cap values were
used, whether a PIH bundle was configured) and returns the two strings.
No `python-docx` import in this module at all, mirroring how
`output/report_v2.py` already separates "shape the data" from
`output/excel.py`'s "write the workbook" — the narrative logic becomes
independently unit-testable (fixed row dict in, fixed string out, no
`.docx` file needs to exist to test it) the same way
`tests/test_scoring_engine.py` tests scoring without touching Excel.

---

## 4. Evidence Architecture (5.1–5.7)

Maps 1:1 to the seven per-candidate evidence categories in
`config/scoring_engine.example.yaml:33-43` (excluding `protein_hunter`
and `interaction`, which are Final-Score-level blends of the other
categories, not evidence categories themselves) and to the
`CATEGORY_EVIDENCE_SHEETS` grouping already used to build the Excel
`05`–`10` sheets (`output/excel.py:170-182`) — note `07_Evolutionary_Evidence`
already carries **two** categories (`pih_evolutionary` +
`pih_cellular_compatibility`) in the existing 12-sheet budget, per the
Stage 1 directive's own mapping; 5.1–5.7 keeps them as **two** subsections
even though they share one Excel sheet, since they are conceptually and
numerically distinct (caps 10 vs. 5).

Recommended fixed text per subsection (short, factual, reusing existing
glossary phrasing from `INTERACTION_SCORE_EXPLANATIONS` /
`config/scoring_engine.example.yaml` comments where it already exists,
rather than inventing new wording that could drift from the Excel
workbook's own explanations):

> **5.1 Sequence Evidence** (cap 30) — BLAST-based positive/negative
> classification and best-hit identity/coverage/E-value strength (see
> `candidate_source` and `candidate_priority_score`). Populated for
> essentially every candidate that has any BLAST hit at all — this is
> the pipeline's most consistently available evidence category.
>
> **5.2 Functional/Domain Evidence** (cap 20) — shared or complementary
> functional annotation (CDD/Pfam domains, description terms) between
> query and candidate. Populated whenever domain annotation is enabled
> and available for both proteins.
>
> **5.3 Genomic Context** (cap 25) — genomic proximity between query and
> candidate genes, from GFF coordinates; a genomic-distance signal used
> as positive evidence only (a distant candidate is never penalized for
> being far away, only not credited for being close). Populated whenever
> GFF neighborhood data is available for both genes.
>
> **5.4 Interaction Evidence** (caps 15 + 12 + 20 across three sources) —
> external protein-protein interaction evidence: STRING PPI
> (`external_ppi_evidence`), GEO transcript coexpression
> (`coexpression_evidence`), and the optional PIH direct-interaction
> bridge (`pih_direct_interaction`). Populated whenever the corresponding
> optional data source (STRING taxon ID, GEO coexpression, or a PIH
> evidence bundle) is configured for the run; each source is independent
> and any subset may be available.
>
> **5.5 Evolutionary Evidence** (cap 10) — phylogenetic/evolutionary
> profile consistency between candidate and query, sourced entirely from
> the optional ProteinInteractionHunter (PIH) evidence bridge. **This
> category requires a `pih_evidence_bundle` file to be supplied via
> configuration; when one is not supplied, every candidate's Evolutionary
> Evidence is reported as not evaluated, not as a scored zero, and it
> does not affect any candidate's ranking either way.** {run-level:
> state plainly whether a PIH bundle was configured for this specific
> run, and if not, that this category is blank throughout this report.}
>
> **5.6 Cellular Compatibility** (cap 5) — subcellular localization /
> compatibility evidence, also sourced from the optional PIH bridge, with
> the identical caveat as 5.5: not evaluated (not scored zero) whenever
> no PIH bundle is configured.
>
> **5.7 Negative Evidence** (reserved, no cap currently populated) — per
> the design specification's own definition (§7.7), this category is
> reserved for evidence that directly *contradicts* a candidate/query
> pairing — e.g. incompatible cellular localization, phylogenetic
> inconsistency, or functionally contradictory annotation. **No such
> signal is implemented in this pipeline version. This category is
> reported as not evaluated for every candidate, in every run, with no
> exceptions** — distinct from 5.5/5.6 above, which *can* become
> populated simply by supplying an optional data file; 5.7 cannot become
> populated by any configuration change available today. This is a
> deliberate design decision, not an oversight: an earlier implementation
> attempt used `negative_hit_strength` (shown elsewhere in this report,
> under `candidate_source`) as a stand-in for this category and found, on
> real-data verification, that it measures a different thing — how
> broadly a candidate protein is conserved across negative reference
> genomes (a phylogenetic-novelty signal) — and using it as a
> contradiction penalty incorrectly punished well-conserved true
> interaction partners (see `claude/final_score_integration_investigation.md`
> and `output/excel.py`'s `final_score` glossary entry). The two signals
> are kept visibly separate in this report for that reason, not merged
> back together for convenience.

This treatment follows §39's two instructions together rather than
picking one: existing functionality (the `negative_hit_strength`/
Negative-Evidence distinction, deliberately re-established after a prior
design confusion) is not quietly re-merged for the sake of a tidier Word
section, and the genuine limitation (5.5–5.7 sparse or empty in most
runs) is stated plainly with *why*, rather than omitted or euphemized.

---

## 5. `requirements.txt`

Current file (`requirements.txt:1-19`) uses one library per logical
group, `>=`-floor pins, grouped with a one-line comment where the
library's purpose isn't self-evident from its name. `python-docx` (PyPI
name `python-docx`, import name `docx`) is not present — confirmed
already in Stage 1's investigation (`import docx` fails with
`ModuleNotFoundError` in this project's `.venv`).

Current PyPI stable release: **1.2.0** (released 2025-06-16, requires
Python ≥ 3.9 — compatible with this project's Python 3.12 requirement per
`AGENTS.md`). Recommended addition, matching the file's existing style:

```text
# Word report generation
python-docx>=1.2
```

`>=1.2` rather than an exact pin (`==1.2.0`) to match every other entry in
this file (all `>=`, none exact-pinned) — re-verify the actual latest
version with `pip index versions python-docx` at M1 implementation time
in case a newer release has shipped by then, rather than trusting this
document's snapshot indefinitely.

---

## Design proposal: M1–M6 for Stage 2

Mirrors Stage 1's own M1–M3 shape (small, test-covered, sequential
milestones) and picks up exactly where Stage 1's original M1–M6 sketch
(in the Stage 1 document, written before Stage 1 itself was scoped down
to Excel-only) left off for the Word half.

- **M1 — `python-docx` dependency + capability smoke test.** Add the
  `requirements.txt` line above, confirm it installs cleanly in the
  existing `.venv` (no known heavy transitive dependencies). Write a
  throwaway smoke-test document that exercises exactly the primitives
  Stage 1's document identified as needed and non-trivial: a heading-based
  Table of Contents field, at least one `w:bookmarkStart`/`w:bookmarkEnd`
  pair, and an external hyperlink (`docx.oxml` escape hatch) — confirm
  all three render and navigate correctly when opened in this machine's
  actual Word/LibreOffice before committing to the cross-link design in
  M5. Resolve the TOC "update field" UX question from item 1 here.
- **M2 — deterministic narrative module, no `python-docx` involved.**
  New `output/word_narrative.py` (or similar) implementing section 3's
  two template functions as pure functions of a row `dict` + run-level
  context. Fully unit-testable in isolation (`tests/test_word_narrative.py`,
  fixed rows in, fixed strings out) — this milestone can be built and
  fully tested before M1's `python-docx` smoke test is even resolved,
  since it has no dependency on the library at all.
- **M3 — Top-N candidate selection, config-driven.** Implement the
  config key and selection logic from item 2 (rank-based Top-N per query
  + Tier1/2 safety-net) as a pure function over the same consolidated
  rows `build_workbook_sheets()` already produces — no new data source,
  just a filter/sort step reusable by both the Word generator and (if
  useful later) a "top candidates" Excel view. Exact N and the
  safety-net rule confirmed with you first (open questions below).
- **M4 — Word document assembly, content only, no links yet.** One
  `.docx` per run, built from M2 + M3's outputs: title/index page,
  "5. Evidence Architecture" (5.1–5.7, static text from item 4, with the
  one run-level branch on whether a PIH bundle was configured),
  "7. Candidate Ranking" and "8. Candidate Details" sectioned per query
  per item 1, each candidate's two narrative paragraphs from M2, a
  summary table per query (design decision E from Stage 1's document —
  self-contained, doesn't require Excel to be open). No hyperlinks or
  bookmarks yet — confirms the content layer is correct and stable
  before adding the more fragile cross-link layer on top of it.
- **M5 — Excel↔Word cross-links.** Bookmarks per query section (and
  optionally per top-candidate row) in Word; the already-reserved
  `word_report_link` column in `02_Final_Score` (`output/excel.py:107`,
  `output/excel.py:461` — currently always blank, explicitly reserved
  for exactly this by Stage 1) populated with
  `report.docx#<BookmarkName>` external hyperlinks (option A+C from
  Stage 1's document); each Word section prints the companion
  `report.xlsx` filename as plain text (option D, since reliable
  Word→Excel deep-linking isn't available per Stage 1's own findings on
  option B). New config key `paths.output_word` alongside the existing
  `paths.output_excel` (`config.yaml:25`), same directory by default so
  relative hyperlink targets resolve without extra configuration.
- **M6 — real-data run + review.** Re-run the same 5-query
  (MA_0688/MA_4165/MA_3898/MA_3899/MA_4115) config Stage 1 used for its
  own real-data verification, generate both the Excel workbook and the
  new Word report from the same run, and confirm by hand: the Word file
  opens correctly with all query sections present and in the expected
  order, every Excel→Word link and every Word→Excel filename reference
  resolves, the narrative text reads correctly for at least one
  candidate in each tier (including Unclassified, if any exist in real
  data) and for both a `Negative_hit`-bucket candidate and a strict
  `Candidates`-bucket candidate, and file size/section count stay
  reasonable at the chosen Top-N. Document results in CHANGELOG.md the
  same way Stage 1 did, following the same
  backup/override/restore-free diagnostic-config pattern (separate
  `--config` file, `config.yaml` itself untouched, worktree removed
  after use if one is needed).

---

## Open questions (need your decision before implementation)

1. **Top-N value and the Tier1/2 safety-net rule (item 2).** Top 15/query
   plus an uncapped Tier1/2 floor is this document's reasoned starting
   proposal, built from the real POS/NEG mean scores and the tier
   thresholds — not a rediscovery of a number you already had in mind.
   If a precise Tier1–4 count histogram across the full 23,130-row
   real-data run would change this, this session can re-run the same
   `git worktree` + matched-config real-data check Stage 1 used, but it
   touches several external APIs and was not fast last time — say so
   explicitly if you want it run before Stage 2 starts.
2. **Query section ordering (item 1).** Keep the existing
   alphabetical-by-`query_id` order (already deterministic, zero new
   code) or thread `interaction_scoring.query_proteins`' config-file
   order through into the Word/Excel row-building path (small change,
   but a change) — your call, both satisfy §45.
3. **Word TOC field behavior (item 1, M1).** Whether "opens with an
   unfilled TOC that needs one right-click-update in Word" is acceptable,
   or whether the document should force `updateFields` in its settings so
   Word recalculates automatically on open (behavior varies slightly by
   Word version/host app — to be confirmed empirically at M1, not
   guessed here).
4. **`word_report.max_candidates_per_query` config location** — proposed
   as a new key nested under `interaction_scoring` in item 2, matching
   where `max_candidates_per_query` (the existing, structurally similar
   knob) already lives; confirm this placement rather than, say, a new
   top-level `word_report:` config block, before M3.
5. **PIH bundle for the M6 real-data run.** Sections 5.5/5.6's narrative
   text branches on whether a PIH bundle was configured. Confirm whether
   the M6 real-data run should include one (if one exists/is ready) so
   the "populated" branch of that text gets exercised at least once
   before Stage 2 is called done, or whether it's acceptable for M6 to
   only ever exercise the "not configured" branch, same as every prior
   real-data run in this project so far.
