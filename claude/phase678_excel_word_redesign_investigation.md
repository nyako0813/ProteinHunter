# Phase 6/7/8 (design spec numbering): Excel 12-sheet redesign + Word report — investigation

Status: **investigation and design proposal only, nothing implemented**.
Mirrors the process used for Phase 6a/6b (`claude/phase6_external_evidence_design.md`,
`claude/phase6b_coexpression_design.md`) and for the earlier v2 scoring
engine work (`docs/implementation_plan_scoring_v2.md`). This document
covers investigation items 1-6 from the request, then a design proposal
with an M1... milestone split for a **future** implementation phase.

**Terminology note up front**: "Phase 6/7/8" here means the *design spec's
own* phase numbers (see below), not this project's own "Phase 6a"/"Phase
6b" (STRING PPI, GEO coexpression) from the last several sessions. Those
are a coincidental numbering collision — the design spec's Phase 6/7/8 is
about the Excel/Word report redesign specifically, unrelated to evidence
sources.

## 1. The original design spec is not in this repository

Searched the full repository, `.git` history is not needed since the file
was never committed, and common local folders (Downloads/Desktop/Documents)
for `ProteinHunter Integrated System 統合設計書 v1.0`,
`ProteinHunter_v5 × ProteinInteractionHunter 統合設計書 v1.0`, and every
filename variant found in citations. **It is not present anywhere on this
machine that this session can see.** Two files that *reference* it as a
companion document are also missing from the repo:
`docs/integrated_scoring_design.md` (cited by `core/evidence.py`) and
`ProteinHunter_v5_scoring_v2_design.md` (cited by
`docs/implementation_plan_scoring_v2.md`). `docs/implementation_plan_scoring_v2.md`
itself calls it "the ChatGPT design specification" -- this strongly
suggests it was an external document (a ChatGPT conversation/export) that
was read from and cited, but never checked into version control.

What *does* exist is a scattering of paraphrased citations across code
comments and two `docs/*.md` files. Every such citation found, quoted in
full with its source:

- `docs/implementation_plan_scoring_v2.md:3-4`: *"Written against
  `ProteinHunter_v5 × ProteinInteractionHunter 統合設計書 v1.0` (the
  ChatGPT design specification), section 56 ("開発開始時の最初の要求")."*
- `docs/implementation_plan_scoring_v2.md:72-78`, section "6. Explicitly
  out of scope for this change (tracked for a later pass)": *"Per the
  design specification's own phase breakdown (sections 47, Phase
  4/6/7/8/9/10):"*
  - *"Phase 4 -- importing ProteinInteractionHunter's own evidence types
    (orthology, phylogenetic profile, known-interaction databases). The
    two projects remain code-independent; see
    `ProteinHunter_v5_scoring_v2_design.md` section 1 for the reasoning."*
  - *"Phase 5 (partial) -- a true "Unified Score" combining a separate
    ProteinHunter score and a separate InteractionHunter score; only
    ProteinHunter's own evidence is modeled here."*
  - **"Phase 6/7/8 -- the redesigned multi-sheet Excel index/final-score-first
    layout and the single-file Word report with table of contents and
    Excel cross-links. The current change only adds columns to the
    existing `Interaction_*` sheets."** (this is the entire, complete
    description found anywhere of the task at hand -- one sentence, no
    further elaboration of "12 sheets" or a precise definition of
    "final-score-first" anywhere else in the repo)
  - *"Phase 9/10 -- full real-data end-to-end run and an old-vs-new
    score/rank comparison report."*
  - *"Weight/cap calibration against known positive/negative interaction
    pairs (design specification section 37) -- no such labeled data was
    available at the time of this change."*
- `docs/scoring_engine_v2.md:169-171`: *"Phases 5 and later of the
  integration design specification (unified ranking refinements,
  Excel/Word report redesign, PIH-vs-v5 comparison reports, end-to-end
  runs on real data, weight calibration) are not yet started."*
- `CHANGELOG.md:326`: *"`ProteinHunter_v5 × ProteinInteractionHunter 統合設計書
  v1.0` に基づく、相互作用スコアリングの改修。"*
- `CHANGELOG.md:284` / `docs/implementation_plan_sequence_evidence.md:9`:
  *"統合設計書3章・9章の指摘:「非常に強いBLAST hitと弱いBLAST hitが、同じ
  positive hitとして同程度に扱われる可能性がある」。"*
- Section numbers cited elsewhere in code comments (no further text
  quoted at each site, just the number): section 3.4 (externalized
  ruleset YAML, `analysis/functional_complementarity_rules.py:10`,
  `config/functional_complementarity_rules.v1.yaml:8`); sections 11-24
  (evidence-based scoring model overview, `analysis/scoring_engine.py:4`);
  section 21 (tie-break/dense-rank rule, `analysis/scoring_engine.py:232`,
  `tests/test_scoring_engine.py:181`); section 22 (protein_hunter_score /
  interaction_score split, `analysis/interaction_scoring.py:1084`,
  `analysis/scoring.py:25`, `config.py:207`, `config.yaml:244`); section 37
  (weight/cap calibration plan, `config/scoring_engine.example.yaml:12`,
  `docs/scoring_engine_v2.md:154`).

**Bottom line for item 1**: there is no fuller text to quote. "The
redesigned multi-sheet Excel index/final-score-first layout and the
single-file Word report with table of contents and Excel cross-links" is
the complete extent of what this repository knows about Phase 6/7/8. If
you have the original document (or the ChatGPT conversation it came from)
outside this repo, pasting the actual Phase 6/7/8 section text in would
materially change the design proposal below from "reasonable inference"
to "implementing what was specified."

## 2. Current `output/excel.py` structure

### Sheet inventory (as actually written by `write_classification_workbook`)

A full run can produce **up to 21 sheets** (not all always populated --
`Interaction_*` sheets only exist for enabled `candidate_sources`, and
`Interaction_Evidence_Detail`/`Interaction_Neighborhood` are opt-in):

| # | Sheet | Always present? | Row source | Role |
|---|---|---|---|---|
| 1 | `Index` | Yes | `INDEX_ROWS` (static, 10 rows) + `interaction_index_rows()` (dynamic, one row per actually-created `Interaction_*` sheet) | Navigation hub: sheet name, "condition", "meaning", "what to do" per sheet, plus a column-glossary block (see below) |
| 2 | `Candidates` | Yes | `positive_only_records` | Strict candidate set: positive hit, no negative hit |
| 3 | `Candidates_relaxed` | Yes | `candidates_relaxed_records` | Positive hit, no *strong* negative hit (tolerates medium/weak) |
| 4 | `Positive_all_sources` | Yes | `positive_all_sources_records` | Hits all configured positive sources, no negative hit |
| 5 | `Positive_source_summary` | Yes | derived from `positive_all_sources` | Compact view of positive-source breadth per target |
| 6 | `Negative_unmatched` | Yes | `negative_unmatched_records` | No negative hit at all (broader than Candidates) |
| 7 | `No_hit` | Yes | `no_hit_records` | Neither positive nor negative hit -- lineage-specific/novel candidates |
| 8 | `Negative_hit` | Yes | `negative_hit_records` | Any negative hit (superset of the three below) |
| 9 | `Negative_strong_hit` | Yes | `negative_strong_hit_records` | Strong negative hit (subset of `Negative_hit`) |
| 10 | `Negative_medium_hit` | Yes | `negative_medium_hit_records` | Medium negative hit, no strong (subset) |
| 11 | `Negative_weak_hit` | Yes | `negative_weak_hit_records` | Weak negative hit only (subset) |
| 12 | `Interaction_query` | If `interaction_scoring.enabled` | resolved query proteins | Confirms query resolution (protein_id/old_locus_tag/sequence -> resolved record) |
| 13-20 | `Interaction_Candidates`, `Interaction_Candidates_relaxed`, `Interaction_Positive_all`, `Interaction_Neg_unmatched`, `Interaction_No_hit`, `Interaction_Neg_hit`, `Interaction_Neg_strong`, `Interaction_Neg_medium`, `Interaction_Neg_weak` | One per enabled `candidate_sources` bucket (up to 9) | per-query, per-candidate scored pairs from that classification bucket | Ranked interaction candidates within that BLAST-classification bucket |
| 21 | `Interaction_Evidence_Detail` | Opt-in (`evidence_detail_sheet`) | one row per (query, candidate, category, component) for v2, or one wide row per pair for legacy | Component-level audit trail |
| 22 | `Interaction_Neighborhood` | Conditional (produced whenever GFF neighborhood data exists) | genomic-distance-ranked pairs | Distance-only view, independent of scoring model |

(Numbered 1-22 above for reference; `Index` + 10 classification sheets + up
to 9 `Interaction_*` bucket sheets + 2 optional detail sheets = up to 22 in
the most permissive configuration used for the recent calibration
diagnostic run; the *default* config produces `Index` + 10 classification
+ 3 enabled `Interaction_*` sheets (`candidates`/`candidates_relaxed`/`no_hit`)
+ `Interaction_query` = 15.)

### Column structure per sheet family

- **Classification sheets** (`Candidates` and the other 9 BLAST-bucket
  sheets): `EXCEL_COLUMNS`, 41 columns -- `protein_id`, `description`,
  `old_locus_tag`, `total_score`, `score_components`, `score_reasons`,
  domain annotation columns (CDD/Pfam sources, names, accessions,
  descriptions, counts), sequence length, positive/negative hit counts and
  best-hit details, `negative_hit_strength` and per-strength counts,
  motifs, UniProt/AlphaFold reference columns, `notes`, and
  positive-source breadth columns.
- **`Positive_source_summary`**: `POSITIVE_SOURCE_SUMMARY_COLUMNS`, 16
  columns -- a compact projection of the above (id/description/locus,
  `negative_hit` flag, source-breadth columns, best-hit summaries).
- **`Interaction_query`**: `query_id`, `input_protein_id`,
  `input_old_locus_tag`, `resolved_protein_id`, `resolved_old_locus_tag`,
  `sequence_length`, `resolution_status`, `description`, `notes`.
- **`Interaction_*` bucket sheets**: `INTERACTION_PAIR_COLUMNS`, 34 columns
  (38 if `include_sequences_in_excel`) -- query/candidate identity,
  genomic-neighborhood fields (contig/coordinates/strand/distance),
  `interaction_priority_score` and its sub-scores
  (`candidate_priority_score`, `same_gene_neighborhood_score`,
  `co_occurrence_score`, `domain_complementarity_score`,
  `alphafold_readiness_score`), `string_ppi_score` (legacy-only),
  v2-only columns (`scoring_model`, `evidence_tier`,
  `formal_score_available`, `evidence_category_count`,
  `evidence_component_count`, `available_weight_total`), and the Phase 5
  reference-only pair: `protein_hunter_score` (+ components/reasons) and
  `interaction_score` (+ `interaction_evidence_tier`).
- **`Interaction_Evidence_Detail`**: two shapes depending on
  `scoring_model`. v2 (`INTERACTION_EVIDENCE_DETAIL_V2_COLUMNS`, long
  format): one row per component --
  `query_id/query_protein_id/query_old_locus_tag/candidate_protein_id/
  candidate_old_locus_tag/candidate_source/candidate_rank/category/
  component_name/status/raw_value/normalized_value/weight/category_cap/
  is_negative/explanation`. Legacy (`INTERACTION_EVIDENCE_DETAIL_LEGACY_COLUMNS`,
  wide format): one row per pair with the five legacy sub-scores plus
  `string_ppi_score` and `interaction_score_reasons`.
- **`Interaction_Neighborhood`**: `INTERACTION_NEIGHBORHOOD_COLUMNS`, 18
  columns -- query/candidate genomic coordinates, `distance_bp`,
  `strand_relation`, `neighborhood_band` (a coarse near/mid/far bucket),
  independent of any scoring model.
- **`Index`**: 4 columns (sheet name / condition / meaning / what-to-do),
  with the sheet-name cell hyperlinked to `#'SheetName'!A1` (openpyxl
  internal-link pattern, see `_format_index_worksheet`/`_add_back_to_index_link`),
  followed by three appended glossary blocks in the same sheet:
  `NEGATIVE_EVIDENCE_EXPLANATIONS` (9 rows), `INTERACTION_SCORE_EXPLANATIONS`
  (13+ rows, growing with every phase -- STRING and coexpression columns
  were appended here in Phase 6a/6b), and `INTERACTION_SCORE_NOTES` (a
  single "Selection rule" cell joining several caveat/attribution
  sentences, including the CC BY 4.0 STRING credit and the GEO PMID
  credits added in Phase 6a/6b).

Every non-Index sheet gets a `"Back to Index"` hyperlink in cell A1
(`_add_back_to_index_link`), and the Index sheet's own sheet-name column
links forward into each sheet's A1 -- i.e., **bidirectional Excel
navigation between Index and every other sheet already exists today**,
via openpyxl's internal `#'SheetName'!A1` hyperlink syntax. This is the
same mechanism that would extend to Excel<->Word links (see item 6).

## 3. Word-generation library availability

**Not present.** `requirements.txt` currently lists: `biopython`,
`pandas`, `openpyxl`, `xlrd`, `numpy`, `requests`, `urllib3`, `PyYAML`,
`tqdm`, `xlsxwriter`, `typing_extensions`, `colorama`. No `.docx` library
of any kind. Confirmed `import docx` fails in the project's `.venv`
(`ModuleNotFoundError`).

**Recommendation: `python-docx`** (PyPI name `python-docx`, import name
`docx`). MIT-licensed, actively maintained, the de facto standard for
programmatic `.docx` generation in Python. Directly supports everything
this phase needs: paragraphs/headings (for a generated table of contents
field or a manually-built one), tables (for a compact per-query summary,
mirroring what a "final-score-first" sheet would show), and -- via its
underlying `docx.oxml` escape hatch -- internal bookmarks and external
hyperlinks, which `python-docx` doesn't wrap at a high level but which are
a well-documented ~15-line helper pattern (add a `w:hyperlink`/`w:bookmarkStart`
element directly). An alternative, `docxtpl` (Jinja2 templating layered on
`python-docx`), is worth a mention but is optimized for filling a
human-designed template with data, not for programmatically generating a
variable number of per-query sections with dynamically-numbered bookmarks
-- `python-docx` directly is the better fit here. No other candidate
(`docx2python` is read-only; `pandoc`/LibreOffice-headless would add an
external-binary dependency this project has avoided so far, e.g. it
already implements its own Excel formatting rather than shelling out).

## 4. What "12-sheet structure" concretely means

**Not specified anywhere found** (see item 1) beyond "redesigned
multi-sheet Excel index." No document or comment gives a sheet-by-sheet
list, so the number 12 itself cannot be traced to source text in this
repo -- it may come from a part of the original design spec that was
never paraphrased into a committed file, or it may be information you
have from the original document/conversation that hasn't been shared with
this session. Reporting this now rather than presenting a guess as fact.

Given that, section 5 of this document (design proposal) includes a
concrete 12-sheet proposal *derived from the current ~21-sheet reality
above*, built by consolidating the most obviously redundant/overlapping
sheets. It should be treated as a **starting proposal for you to
confirm or correct against the real spec**, not a rediscovery of it.

## 5. What "final-score-first" means

Also not specified beyond the one-sentence citation. Two readings are
both consistent with that phrase and with how this codebase already
thinks about scores:

- **(a) Row-order reading**: within each Interaction sheet (or a
  consolidated one), sort so the highest-scoring candidate appears first,
  rather than any other implicit order. This already happens today --
  `candidate_rank`/`interaction_priority_score`-descending sort is the
  existing default (`ranking_metric: interaction_priority_score`, with
  `interaction_score` as the query-specific-only alternative added in
  Phase 5 M5) -- so under this reading, "final-score-first" would mostly
  be documentation/naming (making the existing sort order an explicit,
  named design property) rather than new behavior.
- **(b) Structural reading**: the *workbook* itself leads with a summary
  of top-scoring candidates *across* sheets/queries, before the detailed
  per-bucket sheets -- e.g. a new first sheet (or a rebuilt `Index`) that
  shows, per query, the single best-ranked candidate (or top 3-5) with
  its score and a jump-link into the detailed sheet, so a reader sees
  "here are the results that matter" before wading into 10+ classification
  sheets. This is the more literal reading of "index... layout" combined
  with "final-score-first," and is also the natural shape for the "single
  page per query" content a Word report's table of contents would want to
  mirror.

These are not mutually exclusive -- (a) is almost free (mostly already
true) and (b) is the substantive design work. The proposal in section 6
assumes both: existing sort order kept/confirmed, plus a new
summary-first sheet.

## 6. Excel<->Word cross-link technical options

| Option | How | Pros | Cons |
|---|---|---|---|
| **A. Excel -> Word via external hyperlink** | `cell.hyperlink = "report.docx"` (openpyxl supports external targets, not just the internal `#'Sheet'!A1` form already used for Index navigation) | Trivial, same API already in use; opens the Word file in the reader's default app | Cannot deep-link to a specific *location* inside the docx this way -- lands on page 1 unless combined with option C |
| **B. Word -> Excel via external hyperlink** | `docx.oxml` external relationship + `w:hyperlink` pointing at `results.xlsx` (or `results.xlsx#'SheetName'!A1` -- Windows Explorer/Office *sometimes* honors the `#Sheet!Cell` fragment for `.xlsx` targets, but this is not a documented, cross-platform-guaranteed behavior, unlike Excel's own internal links) | Same idea, reverse direction | Deep-linking into a specific sheet/cell from Word is unreliable across OS/Office versions; best treated as "opens the workbook," not "jumps to the row" |
| **C. Word bookmarks + Excel hyperlink into them** | Insert a `w:bookmarkStart`/`w:bookmarkEnd` pair per query (or per top candidate) in the Word doc; from Excel, `cell.hyperlink = "report.docx#BookmarkName"` | This *is* a reliable, standard Office mechanism (unlike the Excel-side fragment above) -- Word bookmark navigation via `#Name` in a hyperlink target works consistently | Word-side only; still can't deep-link Word->Excel by fragment reliably (falls back to option B's "opens the file" behavior) |
| **D. Shared stable IDs, no live hyperlink** | Give each query (and optionally each top candidate pair) a short stable id (e.g. `query_id`, already present in every sheet) and print it in both documents ("See Excel row ID Q3-C7" / Excel column `pair_id`) | Zero fragility -- survives file renames, moves, different machines, doesn't depend on any Office-version-specific link behavior | Manual lookup instead of a click; weakest UX of the four |
| **E. Embed a static table snapshot in Word instead of linking** | Word report includes its own compact table (top N candidates per query) built directly from the same data, rather than only linking out | Self-contained, readable without Excel open at all -- matches "single-file Word report" framing from the spec citation | Word doc can drift from the Excel file if either is regenerated independently; not a "link," a duplication |

**Recommendation**: **A + C + E together**, not any single option alone --
they solve different halves of the problem and the spec citation asks for
both directions ("Excel cross-links" implies Excel documents need to
reference the Word report, and a "table of contents" inside Word implies
Word needs internal navigation, which naturally extends to Excel via
option C once bookmarks already exist for the TOC). D is a reasonable
fallback everywhere a live link isn't feasible (e.g., no good way to
express "jump to this specific `Interaction_Candidates` row" from Word,
since Excel has no per-row bookmark concept the way Word has
`w:bookmarkStart` -- row/cell references are the practical ceiling, not
a named-anchor jump). B is possible but its reliability caveat should be
documented rather than presented as equivalent to Excel's internal
linking.

## Design proposal (M1... milestones for a future implementation pass)

Given items 4-5 are genuinely underspecified in what's available, this
proposal is deliberately staged so the earliest, cheapest milestones don't
depend on resolving the open questions below -- those get resolved (by
you) before the sheet-count/layout work that actually needs them.

- **M1 -- python-docx dependency, no content changes.** Add `python-docx`
  to `requirements.txt`, confirm it installs cleanly alongside the
  existing stack (no known conflicts expected -- it has no heavy
  transitive dependencies). Write a throwaway smoke-test doc to confirm
  bookmarks + external hyperlinks work as described in item 6 on this
  machine/Word version before committing to the cross-link design.
- **M2 -- 12-sheet Excel consolidation proposal, confirmed against you
  before touching code.** Concretely: propose collapsing the 9 possible
  `Interaction_*` bucket sheets into fewer sheets by adding a
  `candidate_source_bucket` column (the information already exists as
  `candidate_source`) rather than one sheet per bucket, and merging
  `Negative_strong/medium/weak_hit` into `Negative_hit` with a strength
  column the same way (they're already subsets, see the `negative_hit`
  redundancy finding from the last PR) -- this alone could take ~21 sheets
  down to roughly 10-12 without losing any information, purely by
  un-duplicating rows that are currently spread across multiple sheets by
  bucket membership. Exact target list to be confirmed with you (open
  question below) before implementation.
- **M3 -- "final-score-first" summary sheet.** A new sheet (or a rebuilt
  `Index`) presenting, per query, the top-ranked candidate(s) with score
  and an internal hyperlink into the detail row (reusing the existing
  `#'Sheet'!A1`-style pattern, extended to target a specific row rather
  than always A1). Reading (a)+(b) from item 5.
- **M4 -- Word report generation, content only.** One `.docx` per run
  (matching "single-file"), built from the same data the M3 summary sheet
  uses -- table of contents (Word's built-in TOC field, populated from
  heading styles: one heading per query), one section per query with its
  top candidates in a table (option E from item 6). No links yet.
- **M5 -- Excel<->Word cross-links.** Bookmarks per query section in Word
  (option C), external hyperlink cells in Excel's summary sheet pointing
  at `report.docx#QueryBookmark` (option A+C combined), and a
  `report.xlsx` reference (option D, printed as text) in each Word
  section since reliable Word->Excel deep-linking isn't available
  (option B's caveat).
- **M6 -- real-data run + review.** Generate both outputs for a real
  config (e.g. the MA_4115/Tier-A queries already used for calibration)
  and confirm the links actually work when double-clicked, on this
  machine, before calling the phase done.

## Open questions (need your decision before implementation)

1. **Do you have the actual Phase 6/7/8 section text from the original
   design spec (or the source ChatGPT conversation)?** Everything in
   sections 4-5 above is inference from one sentence, not a rediscovery of
   a real specification. If the real text exists, it should override the
   M2/M3 proposal directly rather than this session guessing further.
2. **Exact 12-sheet target list.** The M2 proposal above (merge the 9
   `Interaction_*` bucket sheets and the 3 `Negative_*_hit` sub-buckets
   down using existing category columns) is one specific way to reach
   "roughly 12," not the only one -- e.g. `Interaction_Evidence_Detail`
   and `Interaction_Neighborhood` could also be folded into the main
   ranking sheet as extra columns instead of staying separate, which
   would free up sheet-count budget for something else (a dedicated
   summary sheet, a query-comparison sheet, etc.).
3. **Scope: `v2_evidence_based` only, or does `legacy_additive` also need
   the new layout/Word report?** Given the amount of `v2`-only
   infrastructure already in place (evidence tiers, category caps,
   `Interaction_Evidence_Detail`), a "final-score-first" summary makes the
   most sense keyed on `interaction_score`/`evidence_tier`, which
   `legacy_additive` only partially supports (no tiering concept). Worth
   deciding whether Phase 6/7/8 is v2-only, matching how STRING/coexpression
   integration went.
4. **Word report depth**: full per-candidate detail (every row from every
   Interaction sheet, mirrored into Word) vs. top-N summary per query with
   Excel as the place for full detail. Section 6's option E assumes the
   latter (summary + link out) -- a full mirror would make the "single
   file" Word report very large for queries with thousands of scored
   candidates (as seen in the recent calibration run, individual buckets
   reached 2000+ rows).
5. **Cross-link robustness**: should the Excel/Word files assume they
   always live in the same directory (relative hyperlink target, breaks if
   moved independently) or use some other addressing scheme? No
   existing precedent in this codebase to follow either way.
