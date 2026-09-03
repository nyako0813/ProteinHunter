# Experimental protein-interaction data: curation & locus-tag mapping

Status: **curation only, no calibration logic implemented**. This document
records how the 27 rows in
`Methanosarcina_acetivorans_experimental_protein_interactions.xlsx`
(`Experimental_interactions` sheet) were cleaned, expanded, and mapped to
this project's `old_locus_tag` (GCF_000007345.1, strain C2A) before any use
as a calibration/training set. Every locus-tag identification below was
verified against `data/input/genome.gff` directly (gene symbol, or
`product=` description when no gene symbol was assigned), not guessed from
memory.

## Step 1: excluded rows (non-protein partner)

Three rows removed per instruction -- all involve MmcA (MA_0658) and a
small molecule / metal ion, not a second protein:

| Row (original) | Partner | Reason |
|---|---|---|
| MmcA -- Methanophenazine / 2-hydroxyphenazine | redox cofactor | not a protein |
| MmcA -- AQDS | synthetic redox mediator | not a protein |
| MmcA -- Fe3+ | metal ion | not a protein |

## Step 2: a major finding from species verification (instruction 4)

The Notes sheet's caveat turned out to matter a great deal for one whole
cluster of rows. Checked each cited paper's actual organism (web search,
not assumed from the spreadsheet's own "Organism context" column, which
turned out to be inconsistent with the papers themselves in places):

- **"Structure of the ATP-driven MCR activation complex" (*Nature*, 2025)
  was performed in *Methanococcus maripaludis*, not *M. acetivorans*.**
  ([Nature](https://www.nature.com/articles/s41586-025-08890-7),
  [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12176620/)) The
  spreadsheet's own "Organism context" column already hedged this
  ("Methanogen complex; M. acetivorans homologs" rather than a plain "M.
  acetivorans"), and the web search confirms the hedge was necessary: McrC,
  Mmp7, Mmp17, Mmp3, and A2 were purified and structurally resolved in *M.
  maripaludis*. Every row sourced from this paper (the coarse
  McrA/B/C-Mmp3/7/17/AtwA row, and the five specific "direct structural"
  pairs) is about a **different species**, not experimental evidence for
  strain C2A. *M. acetivorans* orthologs were still identified (see Step 4)
  since the mmp-nomenclature genes are real and present in this genome, but
  the physical-contact claims themselves were not tested in C2A.
- The nitrogenase-PII supercomplex cryo-EM paper (bioRxiv 2025.09,
  Lessner lab) **is** confirmed *M. acetivorans*, native purification.
  ([bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.09.675011v2))
  No caveat needed for rows 11-14.
- The Lieber et al. 2014 PLoS ONE HdrD1 paper **is** confirmed *M.
  acetivorans*.
  ([PLoS ONE](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0107563))
  No caveat needed for rows 1-5.
- The "Interface swapping orchestrates carbon transfer in archaeal ACDS"
  (2026) paper **could not be located** by search at all. A closely related
  but distinct 2024 PNAS cryo-EM ACDS paper exists for **Methanosarcina
  *thermophila*** (a different species in the same genus), which raises the
  possibility the spreadsheet's citation is for a similar non-C2A study, or
  a very recent (2026) paper not yet indexed. Either way, this could not be
  confirmed as C2A data -- **flagged 要確認**, on top of a second, unrelated
  problem found in Step 4 (this genome carries two paralogous CODH/ACS gene
  clusters, and it is unclear which one -- if either -- the paper studied).
- The 2022 "Dissertation experimental work" source for the last two rows
  could not be identified/verified independently (no author or title given,
  no search hit). Its own "Organism context" column plainly states "M.
  acetivorans" (unlike the Nature-2025 rows' hedged wording), so it is kept
  as probably-C2A, but flagged 要確認 for source verifiability.

## Step 3: slash-separated rows expanded into individual pairs

Per instruction 2. Two different situations turned up, handled
differently:

- **Genuinely distinct binding partners** (e.g. `Mmp3/Mmp7/Mmp17`,
  `NifI1/I2`): expanded into one row per partner.
- **One coarse "complex forms" summary row duplicating more specific rows
  already in the sheet** (`McrA/McrB/McrC` x `Mmp3/Mmp7/Mmp17/AtwA etc.`;
  `HdrD1` x `ACDS complex`; `ACDS subunits` x `CoFeSP scaffold`): **not**
  cross-product-expanded. Blindly expanding e.g. the McrABC x Mmp3/7/17/AtwA
  row would fabricate 12 specific pairwise claims (e.g. "McrB-Mmp17") that
  no cited source actually states -- the sheet's own more specific rows
  (McrC-Mmp7, McrC-McrA, Mmp3-McrA, ...) already cover what was actually
  demonstrated. These three coarse rows are dropped as redundant/non-specific,
  not converted into pairs.
- `NifDK` x `NifI1/I2` (a structural-level confirmation of the same
  relationship rows 11-12 already report at the NifD-only level) was
  expanded to add the two **new** pairs it contributes (NifK-NifI1,
  NifK-NifI2) that aren't already covered by the NifD-level rows.

## Step 4: locus-tag mapping for gene/complex names (instruction 3)

Verified directly against `data/input/genome.gff`. Two search strategies were
needed: gene symbol (`gene=`) for well-annotated genes, and `product=`
description for genes RefSeq left with only a generic `MA_RS#####` name.

| Name in source | old_locus_tag | How found | Notes |
|---|---|---|---|
| RNAP subunit D | **MA_1111** | GFF `product=DNA-directed RNA polymerase subunit D` | **This locus is annotated as a pseudogene (frameshifted) in the current RefSeq (GCF_000007345.1, 2024-10-21) annotation.** The 2012 paper demonstrated a real, purifiable protein, so either the modern frameshift call is an assembly/annotation artifact, or something has changed since 2012. Locus identity itself is unambiguous (only one gene fits); the pseudogene status is worth a second look before using this pair. |
| RNAP subunit L | **MA_0721** | GFF `product=DNA-directed RNA polymerase subunit L` | Unambiguous, real protein-coding gene. |
| NifD | **MA_3898** | GFF `gene=nifD` | Unambiguous. |
| NifK | **MA_3899** | GFF `gene=nifK` | Unambiguous. |
| NifI1, NifI2 | **MA_3896**, **MA_3897** | GFF: two unnamed `P-II family nitrogen regulator` genes sitting directly between `nifH` (MA_3895) and `nifD` (MA_3898) | No gene symbol confirms which is I1 vs I2 -- assignment (I1=MA_3896, I2=MA_3897) follows operon order and the published NifI1,2 heterotrimer stoichiometry (2 x NifI1 + 1 x NifI2), but the I1/I2 label-to-locus assignment specifically should be treated as a reasonable inference, not confirmed. The genome has **6 additional** "P-II family nitrogen regulator" genes elsewhere (paralogs for the Fe-only/V nitrogenase systems) -- only these two, positionally inside the Mo-nitrogenase operon, are plausible NifI1/I2 candidates. |
| CdhC | **MA_1014** *and* **MA_3862** | GFF `gene=cdhC` | **This genome has two complete, paralogous CODH/ACS (cdh) gene clusters** (MA_1011-1016 and MA_3860-3865), each with its own cdhA/B/C/D + acsC. Both loci are genuinely annotated `cdhC`; a proteomics study cannot distinguish near-identical paralogs by shared peptides alone. This is very likely *why* the source spreadsheet itself already lists two candidate loci for this one partner -- kept as a genuine two-way ambiguity, not resolved to one. |
| CdhA | **MA_1016** / **MA_3860** | GFF `gene=cdhA` (two paralogs, same cluster issue) | |
| CdhD | **MA_1012** / **MA_3864** | GFF `gene=cdhD` (two paralogs, same cluster issue) | |
| CdhE | **no gene named `cdhE` exists in this genome** | -- | Best guess is `acsC` (product `...subunit gamma`, MA_1011/MA_3865), reasoning by elimination against classic CODH/ACS subunit nomenclature (alpha/CdhA, "epsilon"/CdhB, beta/CdhC, delta/CdhD, gamma/CdhE~acsC) -- **not confirmed**, flagged 要確認. |
| Mmp7 | **MA_3992** | GFF `product=methanogenesis marker 7 protein` | High confidence -- found in the mmp operon (MA_3990-3999) flanking `atwA`. |
| Mmp17 | **MA_3993** | GFF `product=methanogenesis marker 17 protein` | Same operon, high confidence. |
| Mmp3 | **MA_3997** | GFF `product=methyl-coenzyme M reductase-associated protein Mmp3` | Same operon, product string names it explicitly -- high confidence. |
| A2 / AtwA | **MA_3998** | GFF `gene=atwA`, `product=...component A2` | Both aliases confirmed by the same product line -- high confidence. |
| McrA | **MA_4546** | GFF `gene=mcrA` | Unambiguous. |
| McrB | **MA_4550** | GFF `gene=mcrB` | Unambiguous. |
| McrC | **MA_4548** | GFF `gene=mcrC` | Unambiguous. |
| Rnf complex membrane proteins | **not resolved to one locus** | GFF confirms a real, adjacent `rnfABCEG` cluster (MA_0659-0663), with MmcA (MA_0658) genomically adjacent to `rnfC` (MA_0659) | The cited paper's own text says "membrane proteins of Rnf" (plural, unspecified) -- there is no single named partner to map, only a plausible candidate set (RnfA/MA_0663, RnfD/MA_0660, RnfE/MA_0662 are the membrane-embedded subunits). Left 要確認 rather than picking one. |
| ACDS subunits / CoFeSP scaffold | not resolved to a specific pair | -- | Too coarse a description to map to one locus pair beyond what row CdhD-CdhC/CdhA/CdhE already covers; treated as redundant/non-specific, dropped (same reasoning as Step 3's coarse-row rule). |

Rows that already had explicit `(MA_XXXX)` in the source (HdrD1, Mer, DnaK,
Hsp20, MtpA, MtpC, MtsF, MmcA) were used as given, per instruction 3 --
not re-verified against the GFF, since the source already provided them.

## Final curated pair table

Columns: pair (both `old_locus_tag`), strict/soft (instruction 5),
confidence, source. "Confidence" values: **明記済み** (locus given directly
in the source spreadsheet), **GFF特定** (resolved here by searching
`genome.gff`), **要確認** (could not be resolved with confidence, or the
species/paper itself is in question -- do not use as a calibration label
without a closer look).

| # | Protein A | Protein B | Class | Confidence | Source |
|---|---|---|---|---|---|
| 1 | MA_0688 (HdrD1) | MA_1014 (CdhC, cluster 1) | strict | 明記済み + GFF特定(パラログ2候補の一方、要確認) | Lieber et al., PLoS ONE 2014 |
| 2 | MA_0688 (HdrD1) | MA_3862 (CdhC, cluster 2) | strict | 明記済み + GFF特定(パラログ2候補の一方、要確認) | Lieber et al., PLoS ONE 2014 |
| 3 | MA_0688 (HdrD1) | MA_3733 (Mer) | strict | 明記済み | Lieber et al., PLoS ONE 2014 |
| 4 | MA_0688 (HdrD1) | MA_1478 (DnaK) | soft | 明記済み | Lieber et al., PLoS ONE 2014 |
| 5 | MA_0688 (HdrD1) | MA_4574 (Hsp20) | soft | 明記済み | Lieber et al., PLoS ONE 2014 |
| 6 | MA_4165 (MtpA) | MA_4164 (MtpC) | strict | 明記済み | Reichlen/Fu et al., J. Bacteriol. 2019 |
| 7 | MA_4165 (MtpA) | MA_4384 (MtsF) | soft | 明記済み | Biochemical Characterization of MtpA/MtpC, 2019 |
| 8 | MA_1111 (RNAP subunit D, incl. DΔD3 construct) | MA_0721 (RNAP subunit L) | strict | GFF特定(D=現行RefSeq注釈でpseudogene、要確認) | Lessner et al., J. Biol. Chem. 2012 |
| 9 | MA_3898 (NifD) | MA_3899 (NifK) | strict | 明記済み | Cryo-EM nitrogenase-PII supercomplex, 2025 |
| 10 | MA_3898 (NifD) | MA_3896 (NifI1) | strict | GFF特定(遺伝子名なし、オペロン順序による推定) | 同上 |
| 11 | MA_3898 (NifD) | MA_3897 (NifI2) | strict | GFF特定(同上) | 同上 |
| 12 | MA_3899 (NifK) | MA_3896 (NifI1) | strict | GFF特定(NifDK複合体としての構造確認、同上) | 同上 |
| 13 | MA_3899 (NifK) | MA_3897 (NifI2) | strict | GFF特定(同上) | 同上 |
| 14 | MA_1012 (CdhD, cluster 1) | MA_1014 (CdhC, cluster 1) | strict | 要確認(論文自体を特定できず、クラスター不明) | "Interface swapping..." archaeal ACDS, 2026(未特定) |
| 15 | MA_1012 (CdhD, cluster 1) | MA_1016 (CdhA, cluster 1) | strict | 要確認(同上) | 同上 |
| 16 | MA_1012 (CdhD, cluster 1) | MA_1011 (acsC≈CdhE?, cluster 1) | strict | 要確認(遺伝子名不在+同上) | 同上 |
| 17 | MA_3864 (CdhD, cluster 2) | MA_3862 (CdhC, cluster 2) | strict | 要確認(同上、クラスター2側の可能性) | 同上 |
| 18 | MA_3864 (CdhD, cluster 2) | MA_3860 (CdhA, cluster 2) | strict | 要確認(同上) | 同上 |
| 19 | MA_3864 (CdhD, cluster 2) | MA_3865 (acsC≈CdhE?, cluster 2) | strict | 要確認(同上) | 同上 |
| 20 | MA_0658 (MmcA) | Rnf膜サブユニット(MA_0660/0662/0663のいずれか、未特定) | soft | 要確認(論文が単一パートナーを特定せず) | MmcA is an electron conduit..., Nat. Commun. 2024 |
| 21 | MA_4548 (McrC) | MA_3992 (Mmp7) | strict | GFF特定(高確信度)、**要確認(生物種: 実験はM. maripaludis)** | Structure of the ATP-driven MCR activation complex, Nature 2025 |
| 22 | MA_4548 (McrC) | MA_3993 (Mmp17) | strict | GFF特定、**要確認(生物種)** | 同上 |
| 23 | MA_4548 (McrC) | MA_4546 (McrA) | strict | GFF特定(MCR本体サブユニット同士のため保存性は高い)、**要確認(生物種)** | 同上 |
| 24 | MA_3997 (Mmp3) | MA_4546 (McrA) | strict | GFF特定、**要確認(生物種)** | 同上 |
| 25 | MA_3998 (A2/AtwA) | MA_4546 (McrA) | strict | GFF特定、**要確認(生物種)** | 同上 |
| 26 | MA_4548 (McrC) | MA_3997 (Mmp3) | soft | GFF特定。Organism context欄は"M. acetivorans"と明記だが出典(学位論文)自体は未検証 | Dissertation experimental work, 2022(未特定) |
| 27 | MA_3992 (Mmp7) | MA_4546 (McrA) | soft | GFF特定、出典未検証(同上) | 同上 |
| 28 | MA_3992 (Mmp7) | MA_4550 (McrB) | soft | GFF特定、出典未検証(同上) | 同上 |
| 29 | MA_3992 (Mmp7) | MA_3997 (Mmp3) | soft | GFF特定、出典未検証(同上) | 同上 |

Pairs 21/26 and 22 share the same locus identity as pairs already listed
(McrC-Mmp7 appears as both #21 [strict, M. maripaludis structural] and #26
[soft, M. acetivorans-context co-purification]) -- kept as separate rows
deliberately, since they come from different sources with different
species-confidence, not because they are independent biological findings.
A future calibration step should decide whether to merge them (e.g. take
the max evidence level) or keep them distinct.

## Rows dropped as too coarse to map (not in the table above)

- `HdrD1` x `ACDS complex` (Lieber et al. 2014) -- redundant with pairs 1-2.
- `ACDS subunits` x `CoFeSP scaffold` ("Interface swapping..." 2026) --
  redundant with pairs 14-19, and shares the same organism/cluster
  uncertainty.
- `McrA/McrB/McrC` x `Mmp3/Mmp7/Mmp17/AtwA etc.` (Nature 2025) -- superseded
  by the specific pairs 21-25; cross-product-expanding it would fabricate
  claims (e.g. "McrB-Mmp17") no source actually makes.

## Summary counts

- 27 source rows -> 3 excluded (non-protein) + 3 dropped (too coarse) + 21
  expanded/kept, yielding **29 curated pairs**.
- Strict-positive: **21 pairs** (rows 1-3, 6, 8-19, 21-25), of which 14
  carry a 要確認/caveat flag (2 CdhC paralog-ambiguity, 6 CdhD-cluster/CdhE
  pairs, 1 RNAP-D pseudogene caveat, 5 McrC/Mmp/McrA species caveat) and
  **7 are fully clean** (rows 3, 6, 9, 10, 11, 12, 13).
- Soft-association: **8 pairs** (rows 4, 5, 7, 20, 26-29), of which 5 carry
  a flag (1 MmcA-Rnf unspecified-partner, 4 dissertation-source-unverified)
  and **3 are fully clean** (rows 4, 5, 7).
- **Fully clean overall (no 要確認/caveat at all): 10 of 29 pairs** -- rows
  3, 4, 5, 6, 7, 9, 10, 11, 12, 13. Recommend starting calibration with this
  clean subset (which already includes the full NifD/NifK/NifI1/NifI2
  supercomplex) and treating the other 19 as a secondary, lower-trust set
  pending the specific follow-ups noted above (paper identification for the
  2026 ACDS paper, a closer look at the RNAP-D pseudogene call, and a
  decision on how to weigh the M. maripaludis-sourced Mcr-activation-complex
  pairs given C2A has confirmed orthologs but the interaction itself was
  not tested in this strain).
