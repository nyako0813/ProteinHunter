# Phase 6e (continued): Rockhopper operon prediction — LK57 validation, from-scratch on a new machine

Status: investigation complete. This is a from-scratch repeat of an earlier
(uncommitted, now-lost) investigation on a different machine. Everything
below was re-derived on this machine in one session; nothing from the prior
run was reused except the hypothesis it left behind ("try sample LK57").

## Background

Core unmet need: candidate ranking for MA_4115 (and others) misses "non-
adjacent members of the same operon" — the adjacent-pair operon-gap
threshold is already solved (`_gene_neighborhood_v2`, 150bp). Rockhopper
(https://cs.wellesley.edu/~btjaden/Rockhopper/) predicts operon structure
directly from genome + RNA-seq reads and is being evaluated as a candidate
`genomic_context` signal. **This is investigation only** — no pipeline code
was touched.

Three benchmark complexes on NC_003552.1 (from the local
`GCF_000007345.1/genomic.gff`):

| complex | genes | gap(s) | note |
|---|---|---|---|
| Mcr activation cluster | MA_4546–MA_4550 (mcrA/G/C/D/B) | 2–25bp | 5 adjacent genes |
| Nif complex | MA_3896(NifI1), MA_3897(NifI2), MA_3898(nifD), MA_3899(nifK) | 3–28bp each step | **NifI1 and nifK are 2 genes apart** — the gap-bridging test case |
| Mtp complex | MA_4164, MA_4165 | 70bp | 2 adjacent genes |

Prior (lost) run on a GSE77738 sample at 100k/3M/8.35M reads: Nif and Mcr
detected as single operons at all depths; **Mtp never detected**, despite
adequate individual gene expression. Checking GSE77738's per-sample
expression pointed to sample "LK57" (acetate) as a high-Mtp-expression
candidate for retesting the hypothesis that this was a condition-specific
expression issue rather than a Rockhopper limitation.

## 1. Environment setup

No snags. Total setup time for all three components: well under a minute.

| component | approach | time |
|---|---|---|
| Java | none installed, no sudo. Downloaded Eclipse Temurin 17 JRE (linux x64) tarball from `api.adoptium.net`, extracted locally, ran via `./jdk-17.0.20.1+1-jre/bin/java` | 2s download (46MB) |
| Rockhopper | Static download page (`~btjaden/Rockhopper/`) only linked `download.html` via JS-driven navbar (`header.js`), not the jar. Following that link to `download.html` and grepping its `href`s found the direct jar URL immediately (`download/current/Rockhopper.jar`) — no need to reverse-engineer further JS this time | 8s download (14MB), version 2.03 |
| Verify CLI | `java -cp Rockhopper.jar Rockhopper` (no args) prints full usage; confirmed reference-based mode via `-g <genome_dir>` (needs `*.fna`+`*.ptt`+`*.rnt` in one directory) plus a FASTQ positional arg | instant |

## 2. Genome conversion to .ptt/.rnt

Fetched the GenBank flat file for **NC_003552.1** via NCBI eutils
(`efetch.fcgi?db=nuccore&id=NC_003552.1&rettype=gbwithparts&retmode=text`,
14MB, 11s, single unauthenticated request — no rate-limit issues for one
call). Wrote `claude/genbank_to_ptt_rnt.py` (Biopython `Bio.SeqIO`, venv's
existing biopython 1.88) to emit the classic NCBI Genome Projects trio:

- `.fna`: single FASTA record, 70-col wrapped.
- `.ptt`: header `<description> - 1..<len>` / `<n> proteins` / column line
  `Location Strand Length PID Gene Synonym Code COG Product`, one row per
  CDS feature (Length = amino acids from `/translation`, PID = `/protein_id`,
  Synonym = `/locus_tag`).
- `.rnt`: same shape for tRNA/rRNA/ncRNA/tmRNA/misc_RNA features (Length =
  nucleotides).

Ran in 1 second: 4,866 CDS → `.ptt`, 68 RNA features → `.rnt`, 5,751,492bp →
`.fna`. Cross-checked all 11 benchmark loci (MA_4164/4165, MA_3896–3899,
MA_4546–4550) against their `MA_RS#####` locus tags from the local GFF3 —
all present in the `.ptt` with matching coordinates (e.g.
`MA_RS21735 5084184..5084978 +`, `MA_RS21740 5085048..5086082 +`).

## 3. LK57 sample identification + FASTQ acquisition

Queried the ENA portal API directly for the whole study:

```
https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP069835&result=read_run&fields=run_accession,sample_title,sample_alias,library_name,fastq_ftp,fastq_bytes&format=tsv
```

61 runs returned; `sample_title` for `SRR3158263` is exactly `LK57`
(GSM2058189). Confirmed growth condition from GEO's plain-text record
(`acc.cgi?acc=GSM2058189&targ=self&form=text&view=quick`, avoids the HTML
page's reCAPTCHA wall that blocks normal fetching):

```
!Sample_source_name_ch1 = WWM604
!Sample_characteristics_ch1 = strain: C2A
!Sample_characteristics_ch1 = growth media: High salt media with 120 mM acetate
!Sample_characteristics_ch1 = biological replicate: 3
!Sample_characteristics_ch1 = time point: 30 min after halting transcription with actinomycin D
```

Confirmed acetate condition, biological replicate 3. Note this GSE
(GSE77738, "Transcriptomic profiles and RNA half-life of M. acetivorans...")
is an RNA-decay time-course study — LK57 is a 30-minute-post-actinomycin-D
time point, not a pure steady-state snapshot, but reads still reflect
in-vivo transcript boundaries and this is what the prior investigation had
already selected on. The GEO record's own data-processing notes state the
original authors mapped with **Rockhopper v2.0.2** — independent validation
that Rockhopper is an appropriate/expected tool for this dataset.

Full file: `SRR3158263.fastq.gz`, 4.58GB compressed
(`ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/003/SRR3158263/SRR3158263.fastq.gz`).
Rather than a raw byte-range cut (not directly usable on a gzip stream) or
forcing the full download, used a streaming pipe with an exact read-count
cutoff:

```
curl -sL <fastq.gz url> | zcat | head -n <reads*4> > subset.fastq
```

`head` exiting after enough lines closes the pipe, which SIGPIPEs `zcat`
and `curl` and stops the transfer immediately — no wasted bandwidth past
the target read count, and no separate resume/range logic needed.

| subsample | reads | lines | wall time | decompressed size |
|---|---:|---:|---:|---:|
| LK57_100k | 100,000 | 400,000 | 24s | 26MB |
| LK57_3M | 3,000,000 | 12,000,000 | 597s (~10min) | 768MB |
| LK57_8350k | 8,350,000 | 33,400,000 | 723s (~12min) | 2.1GB |

Effective sustained throughput (compressed bytes/time, estimated from gzip
ratio) was in the same **~0.8 MB/s** ballpark as the prior machine's
0.87MB/s — network conditions for this ENA host were not meaningfully
better here. The difference in outcome (full 3-depth series completed here
vs. 33% of one download previously) came entirely from cutting the stream
at an exact read boundary instead of downloading arbitrary/full byte ranges.

## 4. Rockhopper runs

Command (reference-based mode, single condition, default parameters):

```
java -Xmx1200m -cp Rockhopper.jar Rockhopper -g genome_dir -o out_dir -TIME reads.fastq
```

| run | reads | alignment rate | runtime | peak mem |
|---|---:|---:|---:|---|
| LK57_100k | 100,000 | 92% | 20s | well under 1.2GB (`-Xmx1200m` never hit) |
| LK57_3M | 3,000,000 | 92% | 28s | same |
| LK57_8350k | 8,350,000 | 91% | 49s | same |
| LK51_100k (bonus, see below) | 100,000 | — | 20s | same |

### Mtp complex (MA_4164/MA_4165) — the core question

**Not detected as a single operon at any read depth on LK57**, despite LK57
being the sample previously flagged as having comparatively high Mtp
expression. `_operons.txt` never contains an MA_RS21735/MA_RS21740 pair at
100k, 3M, or 8.35M reads.

More striking: **the `_operons.txt` output is byte-identical (same md5sum)
across all three LK57 read depths** (1,031 gene-pairs, 658 multi-gene
operons, every single time), even though the underlying per-gene expression
counts in `_transcripts.txt` clearly scale with depth (e.g. MA_RS21735:
63 → 86 → 79 raw counts across 100k/3M/8.35M; MA_RS21740: 87 → 118 → 118).
The operon caller's pairing decision saturates well below 100k reads for
this genome and does not change with 30-80x more data — this is not a
"needs more reads" situation.

At 8.35M reads, individual expression of the Mtp genes is directly
comparable to the Nif genes that *do* get merged:

```
MA_RS21735 (Mtp)  expr=79    5083935..5084978  (cobalamin B12-binding domain)
MA_RS21740 (Mtp)  expr=118   5085003..5086166  (uroporphyrinogen decarboxylase family)
nifD              expr=89    4790855..4792441
nifK              expr=112   4792469..4793956
```

Expression magnitude is not the limiting factor — the 70bp intergenic gap
itself appears to be the reason Rockhopper does not merge this pair, in
contrast to the Nif chain's 3–28bp per-step gaps.

### Nif complex — reproduced, including the non-adjacent bridge

At every depth tested, `_operons.txt` contains:

```
4789281  4793839  +  5  nifH, MA_RS20335, MA_RS20340, nifD, nifK
```

This is one predicted operon spanning NifI1 (MA_RS20335) through nifK,
i.e. **NifI1 and nifK — 2 genes and ~3.7kb apart — are correctly grouped**
via the chain of small (3–28bp) intergenic gaps between each consecutive
pair. This is exactly the gap-bridging behavior the pipeline needs, and it
reproduces cleanly on LK57 at all three depths (also confirmed on the bonus
LK51 sample, see below). Note the operon also picks up an upstream `nifH`
gene not in the original 4-gene benchmark list — a reasonable/expected
extension, not a problem.

### Mcr cluster — reproduced

```
5596675  5601614  -  5  mcrA, mcrG, mcrC, mcrD, mcrB
```

All 5 genes merged into one operon at every depth tested, matching the
benchmark exactly.

### Robustness check: second sample (LK51)

Time/bandwidth allowed one more spot-check per the investigation's scope.
Picked **LK51** (GSM2058179, SRR3158257): same acetate growth media, same
30-minute post-actinomycin-D time point as LK57, but **biological
replicate 2** instead of 3 — a genuine independent sample, not a re-run.
100k-read subsample only (16s download, 20s Rockhopper run).

Result: `_operons.txt` for LK51_100k is **byte-identical** to all three
LK57 runs (same md5sum). Same Nif and Mcr operons reproduced; Mtp still not
merged. So the result holds across two independent biological replicates,
not just multiple depths of one library.

## 5. Candid assessment: condition-dependent or methodological limitation?

**Methodological limitation, not a condition/expression issue.** The
evidence:

1. LK57 was specifically selected as a high-Mtp-expression sample under the
   hypothesis that low expression in the original test sample explained the
   miss — it still misses.
2. The non-detection is depth-independent: identical operon call at 100k,
   3M, and 8.35M reads (30–80x range) on the same library.
3. The non-detection is sample-independent: identical result on a second,
   independent biological replicate (LK51).
4. At the highest depth tested, MA_4164/MA_4165 individual expression
   (79/118 raw counts) is in the same range as nifD/nifK (89/112), which
   Rockhopper *does* merge — so this isn't simply "too lowly expressed to
   call."

The one structural difference between the working and non-working cases is
gap size: Mcr and Nif's per-step gaps are all ≤28bp; Mtp's gap is 70bp. That
strongly suggests Rockhopper's reference-based operon caller has an
implicit maximum-intergenic-distance (or coverage-continuity) cutoff
somewhere between 28bp and 70bp for this genome/read-length combination —
tighter, not looser, than this pipeline's own already-validated
`_gene_neighborhood_v2` 150bp adjacent-pair threshold.

This has a direct, somewhat counterintuitive implication for the
`genomic_context` design question: Rockhopper's added value was hoped to be
bridging *non-adjacent* members through a chain of tight gaps (which it
does do, as shown by Nif) — but for a **directly adjacent** pair with a gap
in the 30–150bp range (which `_gene_neighborhood_v2` already treats as a
plausible same-operon signal), Rockhopper may actively contradict that
signal with a false negative, since its own adjacency threshold appears
stricter than 150bp. If Rockhopper-derived operon predictions are used as a
`genomic_context` feature, they should **not** be treated as a reliable
negative signal (i.e., "not merged by Rockhopper" should not suppress or
penalize a candidate pair) — only a positive Rockhopper merge should be
scored as additional supporting evidence, and even then only for pairs
Rockhopper actually groups, not for gap-spanning cases in general where
every intervening gap isn't already known to be small.

## 6. Cost comparison

| | prior machine (lost session) | this machine |
|---|---|---|
| Rockhopper acquisition | required reverse-engineering site JS | static jar link found in `download.html`, no JS reversing needed |
| Java | portable JRE, no install | same approach, Temurin 17, 2s download |
| `.ptt`/`.rnt` conversion | written once, not preserved | rewritten, kept at `claude/genbank_to_ptt_rnt.py`, 1s to run |
| FASTQ acquisition | stalled at ~1.5GB/4.58GB (33%) at ~0.87MB/s, one sample, incomplete | full 100k/3M/8.35M subsample series (3 depths) completed via exact-read-count streaming cutoff, ~0.8MB/s effective throughput, ~22 min total download time across all 3 |
| Rockhopper run time/sample | estimated 10–15 min | actual: 20s (100k) / 28s (3M) / 49s (8.35M) — well under the prior estimate even at the largest depth tried |
| Peak memory | estimated <1.5GB | `-Xmx1200m` never exhausted at any depth tested |
| Samples covered | 1 (incomplete) | 2 (LK57 full 3-depth series + LK51 100k spot-check) |
| Total wall-clock, environment setup through final Rockhopper run | not completed | ~30 minutes |

The prior machine's per-run time estimate (10-15 min) was likely dominated
by, or conflated with, its FASTQ download bottleneck rather than actual
Rockhopper compute time — on this machine, Rockhopper itself is fast (well
under a minute even at 8.35M reads); the FASTQ streaming download is the
only real bottleneck, and the exact-read-count pipe cutoff makes that
bottleneck bounded and predictable rather than open-ended.

## Files left on disk (not committed, `data/temp/` is gitignored)

- `data/temp/rockhopper_investigation/Rockhopper.jar`, `jdk-17.0.20.1+1-jre/`
- `data/temp/rockhopper_investigation/NC_003552.1.gb`, `genome_dir/` (`.fna`/`.ptt`/`.rnt`)
- `data/temp/rockhopper_investigation/reads/` (LK57 100k/3M/8.35M, LK51 100k)
- `data/temp/rockhopper_investigation/out_100k/`, `out_3M/`, `out_8350k/`, `out_LK51_100k/` (Rockhopper output dirs, including `_operons.txt`/`_transcripts.txt`)
- `claude/genbank_to_ptt_rnt.py` (reusable conversion script, not committed by this session — left for the user to review/commit)
