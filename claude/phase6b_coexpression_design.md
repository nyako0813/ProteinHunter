# Phase 6b: public coexpression evidence — investigation (pre-design)

Status: **investigation only, nothing approved yet**. Written before any code
so decisions can be made from real data, the same process used for
`claude/phase6_external_evidence_design.md` (Phase 6a). Do not start
implementation from this document alone — it records findings and open
options, not approved decisions.

## Background

Phase 6a (STRING PPI evidence, PR #7) is merged. Phase 6a's investigation
found three real GEO datasets for *M. acetivorans* C2A (see that document's
"Investigation: public coexpression data" section):

| Accession | Samples (GEO-listed) | Conditions | Strain |
|---|---:|---|---|
| GSE77738 | 61 | acetate / methanol / TMA growth + RNA half-life | wild-type |
| GSE64349 | 15 | methylated sulfur compounds (DMS/MMPA/MeSH) vs. MeOH | wild-type (+ a Δ*msrH* mutant subset, see below) |
| GSE66445 | 6 | methane vs. methanol | metabolically engineered — **excluded per instruction** |

This document covers GSE77738 and GSE64349 only, downloaded and inspected
directly (not just their GEO metadata).

## Investigation: file inventory

Both series' processed supplementary files were downloaded from
`ftp.ncbi.nlm.nih.gov/geo/series/...` and decompressed:

- `GSE77738_ReadCounts.xls` (Excel 97-2003, 3 sheets: Raw Read Counts, RPKM
  Normalized Read Counts, Upper Quartile Normalized Read Counts) — the actual
  per-sample expression matrix.
- `GSE77738_DifferentiallyExpressedGenes.MultiFactor.xlsx` — DE gene lists
  (edgeR/PoissonSeq/DESeq2), not a per-sample expression matrix.
- `GSE77738_HalfLives.xlsx` — per-gene RNA half-life estimates, not
  expression levels.
- `GSE64349_TableS1_GEO.xlsx` — wild-type DMS/MMPA/MeSH/MeOH comparison,
  Partek Genomics Suite export with embedded per-sample RPKM columns.
- `GSE64349_TableS2_GEO.xlsx` — a **separate** Δ*msrH* mutant-vs-parental
  comparison, same export format.

Only `ReadCounts.xls` and `TableS1`/`TableS2` contain per-sample expression
values; the other files are derived (DE calls, half-lives) and not directly
usable for a correlation matrix.

## Investigation: gene ID system

- **GSE77738** (`Gene Locus` column): `MA0001`, `MA0002`, ... — the classic
  locus tag **without the underscore** STRING/`old_locus_tag` uses
  (`MA_0001`). 4550/4613 rows (98.6%) match `^MA\d{4}$` and convert to
  `old_locus_tag` by simple string insertion (`MA4115` → `MA_4115`) — no
  lookup table needed. The other 63 rows are `MAt####` (tRNA genes), a
  different feature class, correctly excluded.
- **GSE64349 TableS1/S2** (`Feature ID` column): a **mixed** column — genes
  with a common name use the name (`cdc6_1`, `repA`, `ssrA_1`); genes without
  one fall back to the same `MA####` (no underscore) form as GSE77738. Of
  4725 rows, 4659 (98.6%) were resolved to an `old_locus_tag`: `MA####` rows
  converted directly as above, name rows resolved via a symbol→locus lookup
  table built from GSE77738's own `Gene Locus`/`Gene Name` columns (matching
  case-insensitively and stripping a trailing `_<digit>` disambiguator, e.g.
  `cdc6_1` → `cdc6` → `MA_0001`). The 66 unmapped rows are almost all
  `MAt####` tRNAs plus a handful of plasmid genes (`repA`, `pC2A_p3`,
  `pC2A_p4`) absent from GSE77738's chromosome-only gene list.
- Net effect: **no real ID-mapping project is needed** — both series already
  use this organism's classic locus tags (matching `ProteinRecord.old_locus_tag`,
  same conclusion Phase 6a reached for STRING), and the one non-trivial case
  (gene-symbol rows in GSE64349) is resolved with a small lookup table built
  from data already downloaded for the other series, not an external
  resource.

## Investigation: expression value type

Both series ship **raw counts and RPKM** (not just fold-changes):

- GSE77738: separate sheets for raw counts, RPKM, and upper-quartile
  normalized counts.
- GSE64349: per-sample "Expression values" (raw), "Normalized expression
  values", and "RPKM" columns, plus per-condition "Means"/"Normalized means"
  summary columns.

RPKM is the natural common currency for a cross-series comparison (already
length-normalized; both series are single-end short-read RNA-seq on the same
genome, so RPKM is comparable in principle). log2(RPKM + 1) was used for all
correlation work below, standard practice for count-based expression data.

## Investigation: sample counts and conditions — important correction

**GSE77738 is not 61 independent biological replicates.** Its
`!Sample_characteristics_ch1` metadata (fetched from the series matrix file)
shows this is primarily an **RNA-decay time course**: cells were grown to
mid-exponential phase, transcription was halted with actinomycin D, and RNA
was sampled at 0/5/10/20/30/60/120/240 minutes afterward to compute half-lives.
Actual breakdown of the 61 columns:

| time point | n samples |
|---|---:|
| 0 min (baseline, before decay) | 9 |
| 5 / 10 min | 3 + 3 |
| 20 / 30 / 60 / 120 min | 9 each |
| 240 min | 6 |
| *(no time point recorded — 4 extra steady-state samples)* | 4 |

Only the **9 samples at "0 min"** (3 methanol, 3 acetate, 3 TMA replicates),
plus **4 more samples with no time-point annotation at all** (2 methanol, 2
TMA — likely simple steady-state harvests outside the decay-chase design),
represent independent growth-condition measurements. The other 48 samples
are repeated measurements of a shrinking pool of the *same* starting
material at different times after a transcription-blocking drug — using them
as independent samples in a correlation matrix would mix real coexpression
signal with **shared decay-kinetics correlation** (nearly all transcripts
decline together after actinomycin D, which inflates correlation between
otherwise-unrelated gene pairs; this is visible in the results below).

**GSE64349 does not consist of 15 uniform wild-type samples either.**
`TableS1` (the DMS/MMPA/MeSH-vs-MeOH comparison actually described in the
series summary) has only **9** per-sample RPKM columns (DMS×2, MMPA×3,
MeSH×1, MeOH×3). The remaining 6 of the series' 15 GEO-listed samples are in
`TableS2`, a **different, separate comparison**: 3 replicates of "WWM82
(parental strain)" vs. 3 replicates of "Δ*msrH*" (a deletion mutant). The
series' own summary text ("mRNA from wild-type... compared to that grown on
MeOH") describes only the `TableS1` design; `TableS2` is a mutant-vs-parental
knockout experiment bundled into the same GEO series. This mirrors the
GSE66445 situation the user already flagged for exclusion (metabolically
engineered strain) — Δ*msrH*'s 3 samples should almost certainly be excluded
by the same logic; "WWM82 (parental strain)" is arguably wild-type-equivalent
but is a distinct genetic background used as that experiment's control, not
the same strain as GSE77738/TableS1's C2A — worth a scope decision (see
Design Question 3 below) rather than an automatic include.

**Net usable, genuinely independent wild-type sample count is much smaller
than the "76-82 samples" figure from Phase 6a's high-level GEO metadata
scan**: 9 (or up to 13) from GSE77738, 9 (or up to 12) from GSE64349 — on
the order of **18-25 samples total**, not 76-82. This does not block Phase
6b (see feasibility below), but it changes what statistical care the
normalized_value mapping needs (see Design Question 2).

## Investigation: platform / normalization compatibility

Both series were sequenced on **Illumina HiSeq 2000**, mapped to the same
reference assembly (`NC_003552`/`NC_002097` etc., GCF_000007345.1-equivalent
coordinates), and both ship RPKM. GSE64349 additionally used strain-level
taxid **188937** in its own GEO platform record (`GPL19569`) — the exact
STRING strain taxid Phase 6a's investigation required, confirming both
series really are this project's target organism and not a different
Methanosarcina strain. Technology-wise the two series are compatible; they
are **not** the same experiment design (GSE77738 = growth-substrate
comparison + decay chase; GSE64349/TableS1 = growth-substrate comparison
only) so a naive pool without per-series standardization is not advisable —
tested below.

## Investigation: is a genome-wide gene×gene correlation matrix feasible?

Yes, trivially — this was a real concern to check but not a real
bottleneck:

| dataset | genes | samples | corr matrix compute time | matrix size (float64) |
|---|---:|---:|---:|---:|
| GSE77738, 0-min only | 4550 | 9 | 0.11 s | 165.6 MB |
| GSE77738, all 61 (incl. decay) | 4550 | 61 | 0.12 s | 165.6 MB |
| GSE64349 TableS1 | 4419 (nonzero-var) | 9 | 0.13 s | 156.2 MB |

(`numpy.corrcoef` on log2(RPKM+1), on ordinary hardware.) A full N×N matrix
is not actually needed for the pipeline's use case, though — like the STRING
bridge, only **per-query rows** (one query protein vs. the Candidates
bucket) are ever looked up, which is a single dot-product per candidate and
effectively free. The 150-165 MB full-matrix figures above are informational
only; the recommended implementation (see Design Question 5) never
materializes or caches a full matrix.

## Investigation: MA_4115 vs. neighborhood — is there a usable signal?

Computed Pearson (and Spearman) correlation of log2(RPKM+1) between MA_4115
and its immediate genomic neighbors (MA_4114/4116/4117 — the same genes that
were STRING's top `nscore`/`pscore` partners for MA_4115 in Phase 6a),
against a background of 500 random other genes, in four scenarios:

| scenario | n | MA_4114 r | MA_4116 r | MA_4117 r | background mean r (std) | background P95 |
|---|---:|---:|---:|---:|---:|---:|
| GSE77738, 0-min only | 9 | 0.785 | 0.875 | 0.854 | 0.209 (0.241) | 0.738 |
| GSE77738, all 61 (incl. decay) | 61 | 0.774 | 0.807 | 0.772 | 0.147 (0.236) | 0.667 |
| GSE64349 TableS1 | 9 | 0.720 | 0.944 | 0.934 | **0.760 (0.179)** | 0.909 |
| Combined, z-scored per series | 18 | 0.753 | 0.909 | 0.894 | 0.487 (0.185) | 0.852 |

Two findings worth flagging clearly:

1. **The signal is real and strong in absolute terms** — MA_4115's
   correlation with its known neighbors (0.72-0.94) is consistently far
   above the *median* random gene pair, in every scenario. This corroborates
   STRING's own finding (Phase 6a) that MA_4114-4117 are MA_4115's strongest
   partners, from a completely independent data source (measured expression,
   not curated/predicted PPI).
2. **GSE64349 alone has a badly inflated background** (mean random-pair
   r = 0.76!, P95 = 0.91). With only 9 samples split across 4 growth
   conditions, most gene pairs move together simply because expression
   shifts happen in a few large, shared blocks (one block per condition
   switch) rather than because of real regulatory coupling — a classic
   small-n/few-groups correlation artifact. Against that background,
   MA_4115-MA_4114's r=0.72 is barely above the P50, not obviously
   meaningful, even though it is a "large" correlation in absolute terms.
   GSE77738 (more samples, more independent variation) does **not** have
   this problem (background mean 0.15-0.21). This is the concrete reason a
   fixed r→normalized_value mapping is not recommended — see Design
   Question 2.
3. Including the actinomycin-D decay time points (GSE77738 "all 61") does
   **not** obviously inflate MA_4115's own neighbor correlations much
   (0.77-0.81 vs. 0.79-0.88 for 0-min only) but it does compress the
   background tighter around a lower mean — consistent with shared decay
   kinetics adding correlated noise across all genes rather than a large
   spurious boost specifically for real neighbors. Not disqualifying, but a
   reason to prefer the 0-min-only (or a decay-corrected) subset if simplicity
   is preferred over using every sample.

## Investigation: GEO data terms of use

Attempted to fetch NCBI's GEO disclaimer/FAQ pages directly
(`ncbi.nlm.nih.gov/geo/info/disclaimer.html`, `.../faq.html`); both requests
were intercepted by an NCBI-side reCAPTCHA bot-check from this environment's
network, so the exact current wording could not be verified live. Stating
what is reliably known instead: GEO is an NCBI/NIH public database, records
are U.S. government/publicly-funded work with no login or usage fee, and
NCBI's long-standing published position is that GEO data may be freely
downloaded and reused. This is a materially different situation from STRING
(explicit CC BY 4.0 license, attribution normatively required) — GEO reuse
is unrestricted, but citing the original submitters is standard scientific
courtesy, not a license term. **Recommend**: if this ships, credit both
source publications by PMID (27852217 for GSE77738, 25691524 for GSE64349)
next to the existing STRING attribution on the Index sheet, consistent with
how this project already credits data sources, but do not represent it as a
license requirement the way the STRING CC BY 4.0 line is. Suggest someone
manually re-check `ncbi.nlm.nih.gov/geo/info/disclaimer.html` in a normal
browser before shipping, since this could not be verified programmatically.

## Design questions (open — none of this is decided)

These mirror the six items the request asked about. Presented as options
with a recommendation where the investigation above points somewhere
clearly, not as decisions.

### 1. Category / component naming

Recommend a **new v2 category**, e.g. `coexpression_evidence`, separate from
Phase 6a's `external_ppi_evidence` — not a fusion into it. Reasoning: STRING
scores in `external_ppi_evidence` are, per Phase 6a's own investigation,
"most likely dominated by transferred evidence, not primary data for
*M. acetivorans* itself" for `ascore`. What this phase adds is the opposite:
directly measured, organism-specific transcript coexpression. Conflating the
two under one cap would blur a meaningful distinction the pipeline already
draws elsewhere (e.g., `genomic_context`'s own GFF-distance component vs.
STRING's `string_neighborhood` are kept as two components specifically
because they're "conceptually related but methodologically different").

Component name(s): depends on Design Question 3's answer (one combined
component, e.g. `geo_coexpression`, vs. two, e.g.
`geo_coexpression_gse77738` / `geo_coexpression_gse64349`).

### 2. Mapping correlation → normalized_value (0-1)

**Not** a fixed linear/clamped map of raw r (e.g. `max(r, 0)`, or STRING's
`score/1000` style) — Finding 2 above shows this would treat GSE64349's
inflated-background r=0.72 the same as GSE77738's well-separated r=0.72,
when they mean very different things statistically.

**Recommended**: compute, per query gene, its correlation against *all*
other genes in that dataset (already free — same per-query row computation
used for the lookup itself), and express `normalized_value` as that pair's
**percentile rank** within the query's own background distribution (e.g. a
candidate at the 95th percentile of MA_4115's correlation distribution →
normalized_value ≈ 0.95). This self-calibrates per dataset and per query
gene without needing a hand-picked threshold, and was directly validated
above (the background distributions were computed exactly this way for the
feasibility check). Cost: still O(n_genes) per query, same order of
magnitude as the raw lookup, no meaningful overhead.

Simpler fallback if percentile-ranking is judged too much complexity for
this phase: `max(r, 0)` with an explicitly documented caveat (like STRING's
provisional cap) that this is not corrected for GSE64349's small-n inflation
— acceptable as an MVP but should say so plainly in the component's
docstring, the way Phase 6a explicitly flagged `escore`/`dscore` sparsity
rather than silently treating them as fine.

### 3. Handling two series — separate components, or pooled?

Both are technically viable; genuinely open, unlike Question 2.

- **Separate components** (`geo_coexpression_gse77738` +
  `geo_coexpression_gse64349`, sharing one category cap, mirroring
  `string_cooccurrence`/`string_neighborhood`): simpler to reason about and
  cache independently; GSE64349's noisier background stays contained to its
  own component instead of diluting GSE77738's cleaner signal.
- **Pooled, per-series z-scored** (tested above, n=18 effective samples):
  more statistical power (moderately larger n), and MA_4115's neighbor
  correlations stayed strong (0.75-0.91) after pooling — but the pooled
  background mean (0.487) sits between GSE77738's clean 0.15-0.21 and
  GSE64349's inflated 0.76, i.e. pooling **imports** some of GSE64349's
  inflation rather than diluting it away. Also couples the two datasets'
  caching/invalidation together.

Recommendation if forced to pick one: **separate components**. GSE64349's
inflation problem is a property of that dataset specifically (few samples,
few conditions) and is better handled by scoping it to its own
percentile-normalized component (Question 2's fix already neutralizes most
of the practical harm) than by mixing it into a pooled matrix.

### 4. MISSING vs. "evaluated, weak signal"

No new status is needed — this is exactly the pattern Phase 6a's decision 5
already established for STRING, and the same logic applies directly:

- **MISSING**: the query's or candidate's `old_locus_tag` is not present in
  that GEO dataset's gene list at all (e.g., filtered out during that
  study's own QC, or one of the small number of unmapped tRNA/plasmid
  features found above).
- **AVAILABLE**, always, once both genes are present in the dataset — even
  when the resulting percentile/correlation is low. A pair with a
  near-median or negative correlation is "measured, no coexpression signal
  found," which is informative and distinct from "never measured," exactly
  the same distinction already coded for STRING's cooccurrence/neighborhood
  lookup (`_string_cooccurrence_status_and_value` /
  `StringPpiBundle.lookup`).

### 5. Caching

Confirmed the STRING pattern (Phase 6a decision 6) is the right template and
translates directly:

- Raw supplementary files (`GSE77738_ReadCounts.xls`,
  `GSE64349_TableS1_GEO.xlsx`) downloaded once and kept as **ordinary local
  files** (not through `JsonCache`, which is JSON-only) under something like
  `cache_dir/geo_coexpression/`, re-parsed locally on demand — same reasoning
  as STRING's bulk `.txt.gz` files.
- **No full N×N matrix is ever cached** — at ~4500 genes it would be ~150+
  MB of float64 per dataset, far too large for `JsonCache`'s pretty-printed
  JSON files, and unnecessary: only the query gene's own row is ever needed.
- A **per-query** result (`{candidate_old_locus_tag: normalized_value, ...}`)
  cached via the existing `JsonCache` pattern, namespace `geo_coexpression`
  (or per-component namespace if Question 3 goes with separate components),
  keyed by `f"{dataset_id}:{query_old_locus_tag}"` — small, same shape as
  `string_ppi`'s cache entries.
- Unlike STRING, there is no live-API fallback to design — GEO's
  supplementary files are the only source, so a missing/unreachable file
  degrades straight to "coexpression evidence unavailable" (NOT_RUN or a
  warning), same "optional, best-effort" contract as the rest of Phase 6.

### 6. `legacy_additive` applicability

Defer, same as Phase 6a decision (STRING's `legacy_additive` counterpart was
its own M4 follow-up, done after the v2 components were live and calibrated).
No reason to do this earlier than v2 for coexpression either.

## What this investigation does *not* cover yet

- The `WWM82`/`Δmsrh` inclusion question (Design Question 3's neighbor issue,
  scope decision needed: is `TableS2`'s parental-strain subset in-bounds at
  all, given it's a different construct from GSE77738's plain C2A?).
- Real integration against the current 148-protein Candidates bucket (Phase
  6a's STRING investigation did this cross-check; not repeated here since no
  implementation exists yet to test against).
- Exact GEO license wording (blocked by bot-check; flagged above for manual
  verification before shipping).

## Investigation artifacts

Raw downloaded files and the verification script used above are under
`.cache/geo_investigation/` (gitignored `.cache/` directory, not part of the
repo) — kept for reference during design discussion, not meant to be
committed.
