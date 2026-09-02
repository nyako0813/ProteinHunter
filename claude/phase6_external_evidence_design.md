# Phase 6a: STRING PPI evidence — investigation & design record

Status: investigation complete, M1-M5 approved, implementing on
`feature/string-ppi-evidence-phase6a`. This document records what was
actually verified against live data before writing any code, so the
reasoning behind each decision below survives independently of chat
history.

## Background

Phase 5 (`protein_hunter_score` / `interaction_score` split, PR #6) is
merged. Phase 6 adds evidence from external knowledge bases to
`interaction_score`. The original plan had four legs:

1. STRING (known/predicted PPI database)
2. Public coexpression data for *M. acetivorans*
3. Gene fusion ("Rosetta Stone") evidence via STRING's `fusion` channel
4. Sequence coevolution / direct coupling analysis (DCA/EVcomplex/GREMLIN)

(4) was decided against without further investigation: per-pair coupling
analysis needs its own homolog collection (hundreds of sequences) and
statistical fitting per pair, comparable in cost to AlphaFold3 or worse for
non-model organisms. Treated the same as AlphaFold3 — manual, out of
pipeline scope.

(1)-(3) were investigated against live STRING/NCBI data before designing
anything. Findings below.

## Investigation: STRING

### Taxonomy ID pitfall

The species-level NCBI taxid for *M. acetivorans* (**2214**) returns
nothing from STRING. STRING only recognizes the **strain-level taxid
188937** ("Methanosarcina acetivorans C2A"), which matches this project's
target assembly (GCF_000007345.1). Verified via `get_string_ids` — species
2214 returns `[]` for every identifier tried; species 188937 works.

### Identifier mapping

Tested `get_string_ids` with several identifier forms for MA_4115:

| identifier tried | species | result |
|---|---|---|
| `WP_011024006.1` (RefSeq, versioned) | 188937 | no match |
| `WP_011024006` (RefSeq, unversioned) | 188937 | no match |
| `MA4115` (no underscore) | 188937 | no match |
| `MA_4115` (old locus tag, underscore) | 188937 | **matches**, `stringId: "188937.MA_4115"` |

STRING's data for this organism predates the RefSeq WP_ accession
reannotation and is keyed by the classic locus tag — exactly the format
already stored in `ProteinRecord.old_locus_tag` throughout this pipeline.
No new identifier-resolution logic is needed; `old_locus_tag` is the join
key.

Batch mapping (148 old_locus_tags from the current Candidates bucket, one
POST call): **142/148 (96%) mapped**.

### Live network query for MA_4115

`interaction_partners` for `188937.MA_4115` (default confidence): top hits
are immediate genomic neighbors (MA_4114/4116/4117) and ribosomal proteins,
scored almost entirely from `nscore` (neighborhood) and `ascore`
(coexpression); `escore` (experiments) and `dscore` (databases) are **0**
for every one of MA_4115's top partners.

`network` for MA_4115 against all 142 mapped Candidates-bucket proteins:

- At STRING's default/medium confidence (required_score ≥ 400): **0/142**
  candidates have any recorded association with MA_4115.
- Lowered to STRING's own low-confidence floor (required_score ≥ 150):
  **18/142 (~13%)** show a nonzero score (0.15-0.32, all weak). Every one
  of those 18 has `escore = dscore = fscore = 0`; the signal is almost
  entirely `pscore` (cooccurrence) plus a little `nscore`.

Channel availability across the full candidate×candidate network among
those 142 proteins (1394 edges at required_score ≥ 150):

| channel | edges with signal | share |
|---|---:|---:|
| cooccurrence (pscore) | 1217 | 87% |
| neighborhood (nscore) | 811 | 58% |
| textmining (tscore) | 304 | 22% |
| coexpression (ascore) | 152 | 11% |
| experiments (escore) | 67 | 4.8% |
| databases (dscore) | 19 | 1.4% |
| **fusion (fscore)** | **0** | **0%** |

**Fusion evidence is completely absent for this organism** (0 of 1394
edges, 0 among MA_4115's own top partners). The original plan to use
STRING's fusion channel as a Rosetta-Stone substitute is **retracted** —
there is nothing to substitute with, for this species. Revisit only when
working with a query protein in a better-studied lineage.

STRING's own documentation confirms all evidence types, including
`experiments`/`databases`/`coexpression`, can be "transferred" from other
organisms via homology rather than observed directly in the organism of
interest. For an archaeon with essentially no direct experimental PPI
study, the already-tiny non-zero `escore`/`dscore`/`ascore` values found
above are most likely dominated by transferred evidence, not primary data
for *M. acetivorans* itself.

### Practical implication for MA_4115 specifically

None of the 142 mapped Candidates-bucket proteins have a
default-confidence STRING association with MA_4115 at all. Implementing
this will very likely not change MA_4115's own `interaction_score` results
much, if at all — MA_4115 is annotated by STRING itself as a "conserved
hypothetical protein" with no strong characterized partners. Decided to
implement anyway, since the value is in the general capability for future,
better-characterized query proteins, not a fix targeted at MA_4115.

### Terms of use / access method

- `caller_identity`: not strictly mandatory, "strongly recommended." Should
  identify the tool (e.g. an app name/domain), **not a contact email** —
  the original brief's assumption was wrong on this point.
- No enforced numeric rate limit; documented courtesy guidance is "wait
  one second between calls" and "avoid running scripts in parallel."
- Data license: CC BY 4.0 (attribution required).
- STRING's own guidance: **for anything beyond occasional/limited access,
  download the full dataset instead of repeated API calls.**

### Bulk download (the actual integration point)

Per-organism full-network dumps exist and were downloaded and inspected
directly:

- `https://stringdb-downloads.org/download/protein.links.detailed.v12.0/188937.protein.links.detailed.v12.0.txt.gz`
  — 12,196,138 bytes gzipped (~76 MB decompressed), 1,553,475 rows.
  Space-separated: `protein1 protein2 neighborhood fusion cooccurence
  coexpression experimental database textmining combined_score`, all
  score columns on STRING's raw 0-1000 integer scale. `protein1`/`protein2`
  are already `"188937.MA_####"` — splitting on the first `.` yields the
  `old_locus_tag` directly, no alias table needed.
- `https://stringdb-downloads.org/download/protein.info.v12.0/188937.protein.info.v12.0.txt.gz`
  — 4541 rows, tab-separated `#string_protein_id preferred_name
  protein_size annotation`. Used only to know the full universe of
  old_locus_tags STRING has data for, so a pair absent from the links file
  can be told apart from a protein STRING has never heard of (see M2).
  84 KB.

Only 1246 of the 1,553,475 rows involve MA_4115 specifically (genome-wide,
not restricted to the Candidates bucket) — confirms that a per-query
filtered cache entry stays small (tens to low hundreds of KB) even though
the source file is large.

## Investigation: public coexpression data

Searched NCBI GEO (via eutils, `db=gds`) for *M. acetivorans* expression
profiling. Three real datasets found, all strain C2A (matches this
project's genome):

| Accession | Samples | Conditions | Strain | Supplementary files | Reference |
|---|---:|---|---|---|---|
| GSE77738 | 61 | acetate / methanol / trimethylamine growth + RNA half-life | wild-type | XLS/XLSX (processed) + SRA raw reads (SRP069835) | PMID 27852217 |
| GSE64349 | 15 | methylated sulfur compounds (DMS/MMPA/MeSH) | wild-type | XLSX (processed) | — |
| GSE66445 | 6 | methane vs. methanol growth | **metabolically engineered** (ANME-1 Mcr inserted) | XLSX (processed) | — |

Up to 82 samples total, reasonable condition diversity, processed
(not raw-only) files available for all three. This is a usable resource,
not a "nothing found" result. GSE66445 uses an engineered strain and
should not be pooled with the two wild-type datasets without care.

Not downloaded or processed yet — deferred to **Phase 6b**, a separate
request, to keep this phase's blast radius small.

## Decisions

1. New v2 evidence category **`external_ppi_evidence`**, containing one
   component: `string_cooccurrence` (STRING's `pscore`/cooccurrence
   channel).
2. `string_neighborhood` (STRING's `nscore`/neighborhood channel) is added
   as a **second component under the existing `genomic_context` category**
   (shares its cap), not a new category — same pattern already used for
   `source_classification` + `sequence_evidence`.
3. `fscore` (fusion) — not implemented, evidence retracted (see above).
   `escore`/`dscore` (experiments/databases) — not implemented this phase,
   too sparse and likely homology-transferred for this organism.
   `tscore` (textmining) — not implemented, standard practice to exclude as
   noisy. `ascore` (STRING's own coexpression) — not implemented; if Phase
   6b produces a real *M. acetivorans*-specific coexpression score from the
   GEO data above, that will be strictly more specific than STRING's
   (likely transferred) `ascore`, making it redundant.
4. Coexpression (Phase 6b) is a separate phase — different kind of work
   (data engineering: download, normalize, correlate) from STRING (file
   parse + lookup).
5. Missing vs. evaluated-zero: a pair is `MISSING` when either protein's
   `old_locus_tag` is absent from STRING's `protein.info` for the species
   (STRING has never heard of that protein). A pair where both proteins are
   known to STRING but the pair itself is absent from the links file is
   treated as **evaluated, zero score** (`AVAILABLE`, `normalized_value =
   0.0`) — STRING's cooccurrence/neighborhood methods are computed
   systematically across the whole proteome, so a pair's absence from the
   sparse dump most likely reflects "STRING computed it and got zero
   everywhere," not "STRING never touched this pair." Documented as a
   pragmatic interpretation, not a guarantee from STRING's own docs.
6. Caching: bulk-download-once, not per-pair live calls, per STRING's own
   guidance. The raw `.txt.gz` (and `.info.txt.gz`) files are kept as
   ordinary local files (not through `JsonCache`, which is JSON-only) so
   they are downloaded once and re-scanned locally on demand. A **per-query**
   filtered result (`{candidate_old_locus_tag: {channel: raw_score, ...}}`)
   is cached via the existing `core/cache.py::JsonCache` pattern under a
   `string_ppi` namespace, keyed by `f"{ncbi_taxon_id}:{query_old_locus_tag}"`
   — small (a query typically has a few hundred to ~1200 partner rows in
   the raw file), consistent with how CDD already caches per-protein.
7. A small live-API fallback (`get_string_ids` + `network`, POST,
   `caller_identity` set to an app identifier, self-throttled to ≥1s
   between calls) is implemented for species with no local bulk file yet —
   not the primary path.

## M1-M5 plan (Phase 6a; each step its own commit, tests green)

- **M1**: bulk download + parse + per-query cache (`analysis/string_ppi_bridge.py`), plus the live-API fallback.
- **M2**: `external_ppi_evidence` category, `string_cooccurrence` component. Cap value provisional (documented as such).
- **M3**: `string_neighborhood` component, sharing `genomic_context`'s cap.
- **M4**: `legacy_additive` follow-up (small, mirrors Phase 5's M4).
- **M5**: `Interaction_Evidence_Detail` reflects the new components; STRING CC BY 4.0 attribution added to the Index sheet/docs; real-data verification against the 28-candidate AlphaFold3 calibration set.
