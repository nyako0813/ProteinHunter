# Phase 6d candidate: BioGRID/IntAct interolog inference — investigation (pre-design)

Status: **investigation only, nothing approved or implemented**. Same
process as Phase 6c (`claude/phase6c_kegg_pathway_investigation.md`) and
Phase 6a/6b: verify real data before deciding whether the signal justifies
the implementation cost. Written per
`patches/claude_code_instructions_interolog_investigation.md`, which
itself supersedes a Cowork-side WebFetch pre-survey that hit EBI rate
limits (HTTP 429) and couldn't confirm several claims live — those are
re-verified here from a local environment.

## Background: what interolog inference means here

*M. acetivorans* itself has almost no direct experimental PPI data (Phase
6a's STRING investigation: `escore`/`dscore`/`fscore` all ~0 for this
organism). Interolog inference transfers an experimentally-confirmed
interaction pair from a well-studied organism (E. coli, yeast, human) onto
*M. acetivorans* via orthology: if protein A interacts with protein B in
E. coli, and *M. acetivorans* has orthologs of both A and B, infer the
pair might interact too. This is explicitly an **inference from ortholog
behavior, not direct evidence for this organism** — the same category of
weakness already flagged for `domain_complementarity` ("plausible domain
combination" ≠ "confirmed interaction," see
`claude/alphafold3_candidates_MA_0826.md`'s MA_0826×MA_0361/MA_0363
discussion). Any implementation must carry that caveat in the scoring
design, not just this document.

## Method

- **BLAST reference genome inventory**: read directly from
  `data/databases/positive/` and `data/databases/negative/` (no new
  BLAST run).
- **IntAct**: live REST calls to `https://www.ebi.ac.uk/intact/ws/...`,
  ~1.5s apart, well under any documented limit.
- **BioGRID**: bulk organism-split MITAB zip downloaded directly (no API
  key needed for the bulk file, only for the live REST API); inspected
  with Python's `zipfile`, no extraction to disk needed for the listing.
- **Ortholog proxy for the 5 queries**: the project has no existing BLAST
  reference for E. coli or yeast (see finding 1), and adding one and
  tuning `ortholog_filter` thresholds is explicitly out of scope for this
  investigation. As a stand-in *for this investigation only*, KEGG
  Orthology (KO) group membership — already pulled and cached during
  Phase 6c — was used to check whether a query gene's KO group has any
  E. coli (`eco`) or yeast (`sce`) member, and, where an
  IntAct-interacting partner was found, whether *that* partner's KO group
  has an *M. acetivorans* (`mac`) member. This is a reasonable
  function-based orthology proxy for a feasibility check, but a real
  implementation would use BLAST + the existing `ortholog_filter`
  strong/medium/weak thresholds (`config.yaml`), not KO, once/if E. coli
  or yeast genomes are added to `positive_dir`.

## 1. Existing BLAST reference genome composition (most important finding)

```
data/databases/positive/
  Methanocaldococcus_jannaschii
  Methanococcus_maripaludis

data/databases/negative/
  Aeropyrum_pernix
  Sulfolobus_solfataricus
  thermoplasma_acidophilum
```

**Neither E. coli nor yeast is present.** All five reference genomes
currently in use are archaea (two positive, three negative), chosen for
this project's own conservation-pattern logic (distinguishing "conserved
within closely related methanogens" from "conserved too broadly across
distant archaea"), not as a general-purpose ortholog bridge to
well-studied model organisms. The instruction document's optimistic case
— "ortholog mapping might come free as a byproduct of existing BLAST" —
**does not hold**. Confirmed independently by finding 3 below: even if it
did, BioGRID has no data for either existing positive-reference organism
anyway. A real interolog implementation needs a **new** BLAST reference
genome (E. coli and/or yeast), which is a real cost the instruction
document flagged as the more expensive branch — that is the branch this
investigation landed in.

> **注記 (2026-09-20)**: 上のフォルダ構成は記載どおりだが、当時 `negative/Sulfolobus_solfataricus`
> の中身は実際には *Methanococcus maripaludis*(陽性参照と同一ゲノム)だった
> (`claude/negative_reference_genome_mixup_investigation.md`)。現在は *Saccharolobus
> solfataricus* P2 に差し替え済みで、`config/reference_genomes.v1.yaml` と起動時チェックで
> 構成を検証している。本節の結論(E. coli/酵母のゲノムが無く、interolog には新しい
> 参照ゲノムが必要)は、どのアーキアが陰性側に入っているかに依存しないためそのまま有効。
> 本文は当時の記録として未変更。

## 2. IntAct API connectivity (live, from this environment)

No API key needed — confirmed empirically across 8 successful calls,
zero 429s, ~1.5s spacing. Findings that matter for implementation, beyond
simple "it works":

- **The gene-symbol search endpoint (`findInteractions/{name}`) is
  unreliable for anything but a quick smoke test.** `findInteractions/iscS
  ?taxonId=511145` returned 1,090 results spanning multiple E. coli
  strains (`83333` K12 *and* `83334` O157:H7) despite the `taxonId`
  filter. `findInteractions/thiS?taxonId=511145` returned **1,787,151**
  results dominated by human/Salmonella/unrelated proteins — `thiS` as
  free text apparently matches broadly across the whole database (likely
  tokenized substring matching, not an exact gene-symbol filter), and
  **the `taxonId` query parameter had no observable filtering effect in
  either test** (same result count and species mix with and without it).
- **Searching by UniProt accession instead of gene symbol is precise and
  reliable.** Re-run with accessions:
  - `findInteractions/P0A6B7` (E. coli IscS) → 101 total elements,
    correctly scoped to E. coli-only pairs including the expected
    `iscS↔iscU`, `iscS↔cyaY` (frataxin), `iscS↔fdx` (ferredoxin).
  - `findInteractions/O32583` (E. coli ThiS, the correct accession — an
    initial guess at `P0A891` returned 0 and had to be corrected via a
    live UniProt lookup) → exactly **1** result, `thiS↔thiF`, matching
    the PDB 1ZUD structural pair already catalogued in
    `claude/domain_family_map_v2_sulfur_relay_expansion.md`.
- **Response format**: flat JSON array under `content`, one row per
  interaction record, with `moleculeA`/`moleculeB` (gene symbol),
  `idA`/`idB` (UniProt AC), `taxIdA`/`taxIdB`, `speciesA`/`speciesB`,
  plus PSI-MI method/role metadata. No pagination trouble for
  small/precise queries (`size` defaults to 20/page but `totalElements`
  gives the true count).
- **License/attribution**: no CC-BY-style license text found on the
  IntAct-specific docs (JS-rendered gitbook page, not fetchable by
  WebFetch — same limitation Cowork's pre-survey hit). EBI's general
  terms of use (`ebi.ac.uk/about/terms-of-use`) state attribution is
  expected ("we expect attribution... in publications, services,
  products") and warn that automated use beyond the provided service
  interfaces, or usage heavy enough to degrade service for others, risks
  being blocked — no specific numeric rate limit stated, softer than
  KEGG's hard 3 req/sec cap.

**Conclusion for IntAct specifically**: technically usable without any
registration, but a real implementation must query by UniProt accession
(requires resolving each ortholog to its UniProt AC first — already
available via `annotation/uniprot.py`/`uniprot_bulk.py` for the *M.
acetivorans* side, and via a live UniProt lookup for the E. coli/yeast
side), never by gene-symbol text search.

## 3. BioGRID: API key and bulk file coverage

- **API key**: registration is free and self-service
  (`https://webservice.thebiogrid.org`) but requires a real first
  name/last name/email/project name submitted through a web form —
  correctly identified by the instruction document as something Claude
  Code should not do on the user's behalf. Deferred, as instructed;
  **the user needs to register this personally** if the live REST API
  (rather than the bulk file) is ever needed.
- **Bulk download works without any key.** The URL in the instruction
  document (`.../Release-Archive/BIOGRID-5.0.261/BIOGRID-ORGANISM-5.0.261
  .mitab.zip`) 404s/redirects to an HTML page now (confirms the
  JS-rendering problem Cowork's pre-survey hit — the downloads page is a
  JS app, its real file links aren't in static HTML). The working URL is
  `https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/
  BIOGRID-ORGANISM-LATEST.mitab.zip` — 191 MB, resolves to the same
  release (5.0.261) internally.
- **Contents**: 98 per-organism MITAB files.
  - **E. coli: present**, three variants —
    `Escherichia_coli_K12` (943 B, essentially empty),
    `Escherichia_coli_K12_MG1655` (1.8 MB, the useful one),
    `Escherichia_coli_K12_W3110` (109 MB, huge — likely one giant
    screen's worth of near-duplicate records, would need de-duplication
    before use).
  - **Yeast: present**, `Saccharomyces_cerevisiae_S288c` — 673 MB, by far
    the largest file in the archive (yeast is BioGRID's best-curated
    organism by a wide margin).
  - **Archaea: absent. Zero archaeal species anywhere in the 98-file
    list.** Checked specifically for every organism the instruction
    document asked about — *Haloferax volcanii*, *Sulfolobus* (any
    species), *Pyrococcus*, *Methanocaldococcus jannaschii*,
    *Methanosarcina* (any species) — none present. The only
    non-eukaryotic entries at all are a handful of bacteria (*E. coli*
    variants, *Bacillus subtilis*, *Mycobacterium tuberculosis*,
    *Streptococcus pneumoniae*). This directly answers finding 1's open
    question in the negative: even if `Methanocaldococcus_jannaschii` (already
    a `positive_dir` reference genome) had been usable as a BioGRID
    bridge species, **BioGRID has no data for it at all** — there was
    never a free path through the existing archaea references.
  - License: **MIT**, confirmed via
    `wiki.thebiogrid.org/doku.php/terms_and_conditions` — free
    commercial/automated-pipeline use, only obligation is retaining the
    copyright notice and license text. Matches the instruction document's
    claim and is indeed the least restrictive of STRING (CC BY 4.0),
    KEGG (academic-only, rate-capped), and BioGRID (MIT) — but this
    advantage only matters if BioGRID actually has data to use, and for
    this organism it doesn't.

**Conclusion for BioGRID specifically**: any interolog work would have to
go through E. coli or yeast (both present, well-curated) — never through
an archaeal bridge species, because BioGRID simply has none.

## 4. Real-data estimate for the 5 queries

Using the KO-orthology proxy described in Method:

| query | KO | eco/sce ortholog? | interacting partners (IntAct) | partners with *mac* ortholog | interolog candidates |
|---|---|---|---|---|---|
| MA_4115 | K07585 | none in eco or sce | — | — | not computable |
| MA_0795 | *(none)* | no KO at all | — | — | not computable |
| MA_0826 (`hat`, Elp3-type) | K07739 | **sce:YPL086C (ELP3)** | ELP2, ELP4, ELP6, KTI12 (clean Elongator-complex-only partner list, 129 curated records) | **0 of 4** — none of ELP2/ELP4/ELP6/KTI12 have any KO assignment in `mac` | **0** |
| MA_0074 (`radA`) | K04484 | none in eco or sce | — | — | not computable |
| MA_0266 (RtcB) | K14415 | **eco:b3421 (rtcB)** | 20 raw records, thematically incoherent (ahpC, appY, bluF, cusB, gyrA, groEL, rplC, rpsD, sucB, ...) — reads like generic large-scale screen background, not RtcB-specific biology | **4 of ~19 distinct partners** have a `mac` KO ortholog: GroEL→6 *mac* paralogs (MA_0086/0631/0857/1682/4386/4413), ribosomal L3→MA_1072, ribosomal S4→MA_1109, DNA gyrase A→MA_1583 | **9** raw candidate genes (3 of which — MA_0631, MA_0857, MA_4386 — already sit in the existing `Candidates_relaxed` bucket for unrelated reasons) |

Two clear, opposite failure modes, not one uniform "low signal" result:

- **MA_0826/ELP3**: the ortholog step and the "real curated interaction"
  step both succeed cleanly — ELP3's yeast partners are exactly the
  other Elongator complex subunits, a textbook-clean interaction set.
  But it produces **zero** candidates because archaea (at least this one)
  only carry the catalytic Elp3 subunit and lack the rest of the
  Elongator complex entirely — there is nothing on the *M. acetivorans*
  side to transfer the edge onto. Biologically correct, but useless for
  this pipeline.
- **MA_0266/RtcB**: produces a nonzero count (9), but the E. coli
  interaction list itself looks like large-scale-screen noise rather than
  curated, biologically specific interactions (no thematic coherence
  around RNA repair/ligase biology, dominated by chaperonin/ribosomal
  proteins that show up as promiscuous binders in almost any bacterial
  interactome screen). Naively surfacing these 9 genes as "interolog
  evidence for MA_0266" would be actively misleading — it's the same
  "promiscuous binder as false signal" problem STRING's `textmining`
  channel was excluded for in Phase 6a, just arrived at through a
  different database.
- The other 3 of 5 queries (MA_4115, MA_0795, MA_0074) don't even clear
  step 1 — no E. coli/yeast KO-ortholog found at all, so the method isn't
  computable for them, the same "0/5 have a signal to compute" pattern
  Phase 6c found for KEGG pathways (consistent with all 5 current queries
  being information-processing/genetic-machinery genes rather than
  classical enzymes with deep cross-domain conservation).

## 5. Is this worth implementing?

**No — weakest case of the three external-evidence investigations so
far (STRING implemented, KEGG rejected, this one also rejected, and for
a more fundamental reason than KEGG's).**

KEGG (Phase 6c) was rejected because the *query genes themselves* lacked
the needed annotation (no PATHWAY assignment), a problem that could
plausibly resolve itself for a different future query. Interolog
inference has a **structural** problem on top of that: even where a
query's ortholog and a real curated interaction both exist, the biology
itself may make the inferred edge untransferable (MA_0826/ELP3 case:
architecturally correct, biologically real, and still zero usable
output). That is not a data-coverage problem implementation could fix
later — it's the actual answer.

Cost comparison:

- **STRING (implemented)**: bulk per-organism dump exists and was
  directly usable with the existing `old_locus_tag` join key, no new
  BLAST reference genome needed.
- **KEGG (rejected)**: same join key, no new reference genome needed,
  free bulk-equivalent REST endpoints — rejected purely on signal
  strength (0/5 queries had any PATHWAY).
- **Interolog (this investigation, rejected)**: requires (a) a **new**
  BLAST reference genome (E. coli and/or yeast — not a config tweak,
  an actual new dataset to source, BLAST-index, and validate against
  `ortholog_filter` thresholds), (b) a UniProt-accession resolution step
  for the ortholog before IntAct can be queried reliably, (c) either a
  BioGRID API key (needs the user's personal registration) or careful
  de-duplication of the 109–673 MB per-organism MITAB files, and (d) — the
  real deal-breaker — scoring-design work to keep this from surfacing
  noise like the RtcB/GroEL case as if it were meaningful signal. That is
  substantially more implementation cost than either STRING or KEGG for a
  demonstrated 0-for-3-computable, 1-noisy, 1-clean-but-empty result
  across the current query set.

**Recommendation: do not implement `analysis/interolog_bridge.py` now.**
If revisited, it should wait for either (a) a future query protein that
is a conserved, complex-forming classical enzyme likely to have *both* a
close E. coli/yeast ortholog *and* archaea-conserved partners (unlike
Elongator, which is eukaryote/archaea-asymmetric) — the sulfur-relay
proteins already catalogued in
`claude/domain_family_map_v2_sulfur_relay_expansion.md` (IscS/IscU,
ThiS/ThiF, TusA/BCD/E, SufS/SufE) are plausible future candidates for
exactly this reason, since Phase 6a/6c/6d have now all separately touched
this same protein family and found real signal each time it came up — or
(b) a broader project decision to add E. coli/yeast as permanent BLAST
references for reasons beyond just this feature, which would lower this
feature's marginal cost enough to revisit the calculus above.
