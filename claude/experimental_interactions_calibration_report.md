# Calibration diagnostic report: Tier A/B positive pairs vs. AlphaFold3 negatives

Status: **diagnostic report only, no cap/weight changes implemented**. This
records a single, consistent pipeline run (`scoring_model: v2_evidence_based`,
STRING + GEO coexpression both enabled, all BLAST classification buckets
enabled, `max_candidates_per_query: 5000` so nothing is rank-truncated) with
11 queries: the 4 Tier A/B query proteins (HdrD1, MtpA, NifD, NifK), 6 more
Tier B query proteins (MmcA, McrC, Mmp3, AtwA, Mmp7 -- RNAP subunit D could
not be queried, see below), and the existing AlphaFold3-calibration query
MA_4115, all under identical configuration so results are directly
comparable. Output workbook:
`data/output/ProteinHunter_results_calibration_check.xlsx` (not committed --
regenerable from the config below). Full per-pair data:
`claude/experimental_interactions_calibration_report_pairs.csv` and
`claude/experimental_interactions_calibration_report_negatives.csv`.

## Headline finding: known true positives mostly do not land in the Candidates bucket

Every single Tier A/B partner pair's candidate protein was classified by
this pipeline's own BLAST-based positive/negative logic into
**Negative_hit** (or a stricter sub-bucket) or, less often,
**Candidates_relaxed** -- never into the plain **Candidates** bucket. This
is because HdrD1/Mer/CdhC/MtpA/MtpC/Nif*/Mcr*/Mmp* etc. are ancient,
broadly-conserved central-metabolism/methanogenesis proteins: they BLAST-hit
strongly against this project's *negative* reference genomes too (other
archaea/bacteria that also carry homologs of these core enzymes), which is
exactly what routes a protein into Negative_hit rather than Candidates
under the current ortholog-filter logic. This is a structural
characteristic of the candidate-classification step, not a bug, but it
means: **for well-conserved proteins like these, `candidate_rank` and
`interaction_priority_score` are computed within Negative_hit/Negative_strong_hit
(bucket sizes 2088/780), not within the much smaller Candidates bucket
(151)** -- the numbers below should be read with that in mind.

Despite that, within their actual bucket, most Tier A pairs rank very
highly: e.g. McrC-McrA at rank 1-2 of 780-2088, NifD/NifK-related pairs at
rank 1-9 in several buckets. See the full table below.

## RNAP subunit D (MA_1111) cannot be tested at all

Confirmed directly: `MA_1111` resolved status is **`unresolved`** in the
`Interaction_query` sheet. Checked why -- `data/input/target.faa` (RefSeq's
translated CDS set, which this pipeline's entire protein universe is built
from) has **no entry at all** for MA_1111, consistent with its pseudogene
annotation (frameshifted, no `protein_id` in the GFF's CDS record). This
means MA_1111 cannot be used as a query **and cannot appear as a candidate
under any other query either** -- it is not a low-ranked candidate, it is
completely absent from the protein universe this pipeline works with. This
fully answers the diagnostic question about pseudogene impact: the effect
is total exclusion, not a ranking penalty.

## Tier A pairs (primary calibration basis, 8 query-candidate rows / 6 unique physical pairs)

| Query | Candidate | Label | BLAST buckets (candidate) | Primary sheet | Rank in bucket | interaction_score | interaction_priority_score | evidence_tier | coexpr_GSE77738 | coexpr_GSE64349 | STRING cooccurrence | STRING neighborhood |
|---|---|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|
| MA_0688 (HdrD1) | MA_3733 (Mer) | strict | Negative_hit;Negative_strong_hit | Interaction_Neg_hit | 1433/2088 | **0.000** | 0.000 | Tier4_Weak | 0.009 | MISSING | 0.000 | 0.000 |
| MA_4165 (MtpA) | MA_4164 (MtpC) | strict | Negative_hit;Negative_strong_hit | Interaction_Neg_hit | 106/2088 (2/780 in Neg_strong) | **17.739** | 44.168 | Tier3_Moderate | 0.999 | 1.000 | 0.000 | 0.585 |
| MA_3898 (NifD) | MA_3899 (NifK) | strict | Negative_hit;Negative_strong_hit | Interaction_Neg_hit | 7/2088 (2/780) | **24.958** | 50.907 | Tier2_Strong | 0.360 | 0.787 | 0.340 | 0.756 |
| MA_3899 (NifK) | MA_3898 (NifD) | strict | Candidates_relaxed;Negative_hit;Negative_medium_hit | Interaction_Candidates_relaxed | **1**/1288 | **51.298** | 51.483 | Tier2_Strong | 0.416 | 0.758 | 0.340 | 0.756 |
| MA_3898 (NifD) | MA_3896 (NifI1) | GFF特定 | Negative_hit;Negative_strong_hit | Interaction_Neg_hit | 76/2088 (5/780) | **19.926** | 50.447 | Tier2_Strong | 0.947 | MISSING | 0.229 | 0.722 |
| MA_3898 (NifD) | MA_3897 (NifI2) | GFF特定 | Negative_hit;Negative_strong_hit | Interaction_Neg_hit | 287/2088 (7/780) | **13.333** | 43.529 | Tier3_Moderate | 0.735 | MISSING | 0.000 | 0.802 |
| MA_3899 (NifK) | MA_3896 (NifI1) | GFF特定 | Negative_hit;Negative_strong_hit | Interaction_Neg_hit | 424/2088 (5/780) | **10.940** | 37.717 | Tier3_Moderate | 0.221 | MISSING | 0.258 | 0.651 |
| MA_3899 (NifK) | MA_3897 (NifI2) | GFF特定 | Negative_hit;Negative_strong_hit | Interaction_Neg_hit | 785/2088 (9/780) | **5.528** | 32.473 | Tier3_Moderate | 0.157 | MISSING | 0.000 | 0.720 |

**Tier A summary (n=8):** `interaction_score` mean **39.55**, median **43.85**.
`interaction_priority_score` mean 17.97, median 15.54.
`coexpression_gse77738` mean 0.48, median 0.39.
`coexpression_gse64349` mean 0.85, median 0.79 (n=4, the other 4 are MISSING -- GSE64349 doesn't cover NifI1/NifI2/Mer).
`string_cooccurrence` mean 0.15, median 0.12.
`string_neighborhood` mean 0.62, median 0.72.

## AlphaFold3-confirmed-negative comparison (MA_4115 query, 28 calibration entries, same run)

**n=28:** `interaction_score` mean **12.95**, median **12.27**.
`interaction_priority_score` mean 20.55, median 23.66.
`coexpression_gse77738` mean 0.65, median 0.73.
`coexpression_gse64349` mean 0.60, median 0.74.
`string_cooccurrence` mean 0.05, median 0.00.
`string_neighborhood` mean 0.04, median 0.00.

(Full per-entry table: `claude/experimental_interactions_calibration_report_negatives.csv`.
Note: with the wider bucket/rank settings used for this run, all 28 entries
now land somewhere -- including the 8 that the original calibration file
recorded as "not listed in any sheet" under the narrower default settings.)

## Reading the two tables together (observations, not conclusions -- 判断は保留)

- **`interaction_score` separates the two sets far better than
  `interaction_priority_score` does.** Tier A positives: mean 39.55 vs.
  AF3 negatives: mean 12.95 (~3x). `interaction_priority_score` shows
  almost no separation (17.97 vs. 20.55 -- negatives *slightly* higher).
  This is consistent with `interaction_priority_score` including
  candidate-quality-only components (`source_classification`,
  `co_occurrence`) that reflect "how BLAST-good is this candidate", not
  "does it interact with this query" -- and the AF3-negative set was
  originally drawn from the Candidates/Candidates_relaxed buckets (already
  BLAST-favorable), while most Tier A positives are not. This is an
  argument for `ranking_metric: interaction_score` (already an existing
  config option, Phase 5 M5) over the default `interaction_priority_score`
  when the goal is identifying true interacting partners specifically.
- **`string_neighborhood` is the cleanest single separator found here**
  (0.62 vs. 0.04, ~16x), though this is expected and somewhat circular for
  half of Tier A -- NifI1/NifI2 sit in the same operon as NifD/NifK, so
  genomic proximity alone nearly guarantees a high neighborhood score. It
  says less about HdrD1-Mer or MtpA-MtpC, where the neighborhood score
  (0.000, 0.585) is more informative and still favors the positive
  direction for MtpA-MtpC.
- **`coexpression_gse77738` does *not* separate the two sets in this small
  sample -- if anything it runs slightly backwards** (0.48 positives vs.
  0.65 negatives). This is driven substantially by HdrD1-Mer's very low
  value (0.009, a single low outlier in an n=8 set) and is consistent with
  Phase 6b's own finding that AF3-confirmed negatives already sit at a
  fairly high coexpression percentile in this candidate pool (median 0.73)
  -- i.e., this project's own prior finding that coexpression alone is a
  weak standalone discriminator continues to hold up against a real,
  independent positive set, not just the negative set it was originally
  checked against.
- **`coexpression_gse64349` shows some separation** (0.85 vs. 0.60) but
  only 4 of 8 Tier A pairs have any value at all (GSE64349 doesn't cover
  Mer/NifI1/NifI2), so this is a very small sample to draw anything from.
- **n=8 for Tier A is small.** None of the above should be read as a
  statistically solid signal -- it is a first real look, not a calibration
  fit. HdrD1-Mer in particular pulls several Tier A averages down
  (`interaction_score=0.000`, `coexpression_gse77738=0.009`) despite being
  one of the strongest "明記済み" entries in the curated set; it may be
  worth a closer individual look (its `Interaction_Evidence_Detail` row
  shows every category at 0 except `co_occurrence=1.0`, i.e. essentially
  no positive evidence of any kind was found for this specific pair despite
  it being real, published, direct experimental evidence).

## Tier B pairs (diagnostic only, not used for cap/weight judgment)

17 query-candidate rows across CdhC paralogs, DnaK/Hsp20/MtsF (demoted from
Tier A per your instruction), MmcA-Rnf (not run -- no single locus to
query), and the Mcr-activation-complex cluster (McrC/Mmp3/Mmp7/Mmp17/AtwA/McrA/McrB).
Full table: `claude/experimental_interactions_calibration_report_pairs.csv`
(`tier_final` column). Two points worth flagging without over-interpreting
(these are 要確認-flagged pairs, species-uncertain or source-unverified):

- The Mcr-activation-complex pairs (M. maripaludis-sourced, species caveat)
  score comparably to or higher than Tier A on `interaction_score` in
  several cases (e.g. McrC-McrA: 38.006, Mmp7-Mmp3: 33.667), and their
  coexpression percentiles are mostly very high (0.8-1.0 on GSE77738 for
  several). Whether that reflects real conserved biology or the general
  "many things land at a high coexpression percentile in this candidate
  pool" pattern already seen in the AF3-negative set is not something this
  report resolves.
- The CdhC-paralog pairs (HdrD1-MA_1014, HdrD1-MA_3862) both show high
  coexpression (0.905, 0.807) but zero STRING evidence and a near-zero
  `interaction_score` (2.362, 1.190) -- consistent with Phase 6a's own
  finding that STRING has essentially no data for this organism's
  metabolic-complex proteins beyond genomic neighbors.

## Configuration used for this run

`interaction_scoring.enabled: true`, `scoring_model: v2_evidence_based`,
`string_ppi_ncbi_taxon_id: 188937`, `geo_coexpression_enabled: true`, all 9
`candidate_sources` enabled, `max_candidates_per_query: 5000`,
`evidence_detail_sheet.include_no_hit: true`, 11 `query_proteins`
(MA_0688, MA_4165, MA_3898, MA_3899, MA_1111, MA_0658, MA_4548, MA_3997,
MA_3998, MA_3992, MA_4115). This was a **temporary, local-only** config
change for this diagnostic run -- not committed, and `config.yaml` in the
repository is unchanged.
