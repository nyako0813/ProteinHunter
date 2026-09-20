# Phase 6c candidate: KEGG pathway co-occurrence — investigation (pre-design)

Status: **investigation only, nothing approved or implemented**. Written
before any code, per `patches/claude_code_instructions_kegg_pathway_investigation.md`,
using the same process as Phase 6a (`claude/phase6_external_evidence_design.md`)
and Phase 6b (`claude/phase6b_coexpression_design.md`): verify real data
before deciding whether the signal justifies the implementation cost. This
supersedes the preliminary WebFetch-based estimates quoted in the
instruction document (59% KO / 15% pathway coverage) — those numbers turn
out to be wrong in both directions (see below).

> **注記 (2026-09-20)**: 本文の候補バケツ件数(Candidates 127 / Candidates_relaxed 1,261、
> `old_locus_tag` 解決可能 1,239)は、陰性参照ゲノムの取り違え
> (`negative/Sulfolobus_solfataricus` に陽性ゲノムのコピーが入っていた。
> `claude/negative_reference_genome_mixup_investigation.md`)がある状態の出力
> (`data/output/aa.xlsx`)に基づく。修正後の再分類では Candidates 127→709、
> Candidates_relaxed 1,261→1,668。ただし本調査の結論(全5クエリがPATHWAY割り当てを
> 持たず共起チェックが計算不能、実装は見送り)は KEGG の PATHWAY 被覆に基づくもので、
> バケツの大きさには依存しないため、そのまま有効。バケツ件数を再利用する場合は再取得すること。
> 本文は当時の記録として未変更。

## Method

Bulk, organism-wide REST calls only, one call each, 1s+ apart, per KEGG's
rate guidance (see "Terms of use" below):

- `GET https://rest.kegg.jp/list/mac` — full gene list for *M. acetivorans*
  C2A (KEGG organism code `mac`, matches STRING's taxid 188937 and this
  project's `old_locus_tag` join key; no new identifier logic needed, same
  conclusion as Phase 6a).
- `GET https://rest.kegg.jp/link/ko/mac` — gene→KO links, whole organism.
- `GET https://rest.kegg.jp/link/pathway/mac` — gene→PATHWAY links, whole
  organism.
- `GET https://rest.kegg.jp/get/mac:MA_XXXX` for the 5 query genes only (5
  calls total, to get human-readable KO/BRITE names for the report — well
  within the 3 req/sec limit even without throttling).

The Candidates/Candidates_relaxed bucket used for co-occurrence checking
was read directly from the most recent local pipeline output
(`data/output/aa.xlsx`, sheet `03_Candidate_Overview`), not regenerated:
`candidate_source` there is **query-independent** (consolidated per
candidate protein, not per query — see that sheet's own column
description), so it is the same bucket regardless of which query protein
is configured. 1,261 proteins are classified as `Candidates` (127) or
`Candidates_relaxed` (1,134); 1,239 of those carry a resolvable
`old_locus_tag` (22 lack one and were excluded from the pathway
cross-reference).

## 1. Coverage rate (exact numbers)

| | count | % of all `list/mac` entries (4,662) | % of CDS-only (4,540) |
|---|---:|---:|---:|
| genes with KO | 1,854 | 39.8% | 40.8% |
| genes with PATHWAY | 989 | 21.2% | 21.8% |

`list/mac` breaks down as 4,540 CDS + 112 `gene` (features without a KEGG
protein entry, likely pseudogenes) + 10 rRNA. The local UniProt bulk cache
(`data/databases/uniprot/uniprotkb_methanosarcina_acetivorans.json`) has
4,880 entries — same order of magnitude, no reconciliation needed for this
investigation.

**Correction to the preliminary estimate**: the instruction document's
WebFetch-based numbers (~59% KO, ~15% pathway) were wrong in both
directions — actual KO coverage is meaningfully *lower* (40% vs. an
estimated 59%) and actual pathway coverage is somewhat *higher* (21-22%
vs. an estimated 15%). Directionally this still supports the same overall
picture as the preliminary estimate (KO coverage well under genome-wide
majority, pathway coverage a clear minority), just with different
magnitudes — worth flagging since a summary-model WebFetch pass is not a
substitute for pulling the real bulk data, consistent with why this
investigation was requested in the first place.

## 2. Pathway size distribution (specificity check)

113 distinct pathway IDs across the 989 pathway-annotated genes. Top 10 by
gene count:

| pathway | genes | note |
|---|---:|---|
| path:mac01100 | 571 | **KEGG global/overview map** ("Metabolic pathways") |
| path:mac01110 | 239 | **overview map** ("Biosynthesis of secondary metabolites") |
| path:mac01120 | 239 | **overview map** ("Microbial metabolism in diverse environments") |
| path:mac01200 | 158 | **overview map** ("Carbon metabolism") |
| path:mac00680 | 154 | Methane metabolism (specific) |
| path:mac01240 | 140 | **overview map** ("Biosynthesis of cofactors") |
| path:mac01230 | 115 | **overview map** ("Biosynthesis of amino acids") |
| path:mac02024 | 80 | Quorum sensing (specific) |
| path:mac03010 | 70 | Ribosome (specific) |
| path:mac02010 | 64 | ABC transporters (specific) |

The top of the distribution is dominated by KEGG's own "global/overview"
maps (the `011xx`/`012xx` numbered pathways, which are deliberately broad
composites covering large fractions of core metabolism), not biologically
specific pathways. `path:mac01100` alone covers **12% of the entire
genome** — matching two proteins on nothing but "both appear somewhere on
the global metabolic map" would be close to meaningless, the same
fold-level-over-matching failure mode already fixed once in
`analysis/domain_family_map.py` (v2/v3, see recent commit history) and
exactly the concern the instruction document raised in advance.

Excluding the 12 recognized overview-map IDs (`01100/01110/01120/01130/
01200/01210/01212/01220/01230/01232/01240/01250`), 102 pathways remain,
median size 10 genes, top specific pathway `path:mac00680` (Methane
metabolism) at 154 genes. 984 of the 989 pathway-annotated genes have at
least one non-overview pathway in addition to (or instead of) an overview
one — so a real implementation would need to **exclude the overview-map
IDs from the "same pathway" co-occurrence check**, not just filter by
size, since some overview maps (571, 239 genes) are far larger than any
specific pathway's max (154 genes) and would otherwise dominate scoring
with the least specific signal available.

## 3. Signal for the 5 current query genes

None of the 5 currently configured/candidate query genes have **any**
KEGG PATHWAY assignment:

| query | KO | PATHWAY | KEGG classification |
|---|---|---|---|
| MA_4115 | K07585 (tRNA methyltransferase) | **none** | BRITE: Unclassified — genetic information processing / Translation |
| MA_0795 | **none** | **none** | no ORTHOLOGY entry at all |
| MA_0826 (`hat`) | K07739 (elongator complex protein 3 / tRNA carboxymethyluridine synthase) | **none** | BRITE: Transfer RNA biogenesis (Elongator complex) |
| MA_0074 (`radA`) | K04484 (DNA repair protein RadB) | **none** | BRITE: DNA repair and recombination proteins (HR) |
| MA_0266 | K14415 (tRNA-splicing ligase RtcB) | **none** | BRITE: Transfer RNA biogenesis (tRNA splicing ligase complex) |

Since every query has zero PATHWAY assignments, the co-occurrence check
against the 1,239-member Candidates/Candidates_relaxed bucket is **not
computable for any of the 5 — 0/5, not "low count," literally no query to
run**. This is a stronger negative result than the instruction document's
preliminary MA_4115-only finding suggested: it isn't one weak query among
five, all five current queries share the same pattern (KO present, no
PATHWAY), and it's a coherent pattern, not noise — all five are
information-processing/genetic-machinery genes (tRNA modification/
splicing, DNA repair), the class of gene KEGG's PATHWAY maps are known to
under-represent relative to classical metabolic pathways (matches
`path:mac00680` Methane metabolism topping the specific-pathway list
above — KEGG's PATHWAY coverage for this organism skews toward central
carbon/methanogenesis metabolism, not genetic-information-processing
machinery).

## 4. Terms of use / rate limits

- **Academic use of `rest.kegg.jp` is free** ("Academic users may freely
  use the KEGG website"); **non-academic/commercial use requires a paid
  license** (contact Pathway Solutions). If this tool is ever offered as a
  hosted service to others (vs. run locally for research), an "Academic
  Service Provider License" would additionally be required — not relevant
  to the current local-tool usage pattern, but worth remembering if that
  changes.
- **Rate limit: hard cap of 3 requests/second, enforced** — "Please limit
  your API calls up to 3 times per second. Otherwise, your access will be
  blocked." Stricter and more concretely enforced than STRING's soft
  "wait ~1s" courtesy guidance from Phase 6a.
- **No free bulk/FTP download**, unlike STRING. KEGG's FTP bulk-download
  tier (comprehensive `pathway`/`brite`/`module`/`genes`/`fasta`/etc. data,
  updated weekly) requires a paid **academic subscription contract**; only
  a limited KEGG MEDICUS subset is free via FTP, which doesn't cover
  pathway/genes data. This means the STRING-style "download the whole
  per-organism dump once" pattern from Phase 6a's decision 6 is **not
  available** for KEGG without a subscription.
- **Practical substitute**: `link/ko/<org>` and `link/pathway/<org>` are
  already whole-organism, one-call bulk endpoints (not per-gene calls) —
  this is the officially intended way to get organism-wide coverage
  through the free REST API, confirmed by the response sizes here (1,856
  and 3,053 lines respectively for ~4,500 genes, one call each). The
  investigation's own two calls stayed well under the 3 req/sec cap even
  including the 5 confirmatory `get/mac:MA_XXXX` calls. Caching the result
  of these two calls as local files (same `JsonCache`-adjacent pattern as
  Phase 6a decision 6, keyed by organism code) would still avoid all
  repeat network calls in normal operation — the FTP subscription gap only
  matters if this ever needs cross-organism bulk data beyond `mac`.
- No explicit citation/attribution requirement found on either page
  checked (`kegg/rest/`, `kegg/legal.html`); KEGG's copyright notice
  ("Copyright 1995-2026 Kanehisa Laboratories") appears but no CC-style
  license text like STRING's CC BY 4.0. Cite the standard KEGG references
  (Kanehisa & Goto, and the REST API paper) as scientific courtesy if
  results are published, same reasoning already applied to the GEO
  datasets in Phase 6a/6b.

## 5. Is this worth implementing?

**No, not for the current query set — weaker case than STRING (Phase
6a), which was itself already a marginal call.**

Phase 6a's STRING investigation found *some* real signal to build on: 96%
of the Candidates bucket mapped to STRING IDs at all, and even though
MA_4115 itself had zero default-confidence partners, the general
`cooccurrence`/`neighborhood` channels had substantial edge coverage
(87%/58%) across the candidate×candidate network — the STRING decision to
implement anyway rested on "no signal for this one query, but real
capability for future queries."

KEGG pathway co-occurrence doesn't clear even that bar for the queries in
front of it:

- Every one of the 5 current/candidate queries has **zero** PATHWAY
  assignments — not "weak," literally nothing to compute a co-occurrence
  count from, for all 5, not just one.
- This isn't bad luck — it's structural. All 5 are genetic-information-
  processing genes (tRNA modification, tRNA splicing, DNA repair), and
  KEGG's PATHWAY maps for this organism are visibly metabolism-centric
  (the top specific pathway is Methane metabolism at 154 genes; the
  BRITE-only classification these 5 genes get instead — "Transfer RNA
  biogenesis," "DNA repair and recombination proteins" — exists precisely
  *because* KEGG doesn't have a PATHWAY map for that functional class).
  Unless a future query happens to be a classical metabolic enzyme, the
  same "no signal" result should be expected again.
- Only 21.8% of the CDS-only genome has any PATHWAY at all, so even a
  future metabolic-enzyme query starts from worse odds than KO coverage
  (40.8%) would suggest, and the implementation would additionally need
  bespoke handling to exclude KEGG's own overview/global maps (up to 571
  genes, 12% of the genome, on `path:mac01100` alone) to avoid
  reintroducing the exact low-specificity over-matching problem already
  fixed once in `domain_family_map.py`.
- The access pattern is also strictly worse than STRING's: no free bulk
  dump file (subscription-gated), a harder 3 req/sec cap, and no CC-style
  reuse license text found — all manageable via the two organism-wide
  REST calls used here, but more caveats to document and enforce in code
  than STRING's single `.txt.gz` download ever needed.

**Recommendation: do not implement `analysis/kegg_pathway_bridge.py` (or
equivalent) now.** Revisit only if/when a future query protein is a
recognizable classical-metabolism enzyme (the kind of gene that shows up
in `path:mac00680`/`02010`/`00230`-style specific pathways) rather than
information-processing machinery — at that point coverage is still capped
around 22% genome-wide, but at least the query itself would have a
PATHWAY entry to compute co-occurrence from, which none of the 5 current
queries do.
