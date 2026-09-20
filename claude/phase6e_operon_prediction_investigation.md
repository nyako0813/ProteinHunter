# Phase 6e candidate: direct operon prediction (Rockhopper) — investigation (pre-design)

Status: **investigation only, nothing approved or implemented**. Same
process as Phase 6c/6d. Written per
`patches/claude_code_instructions_operon_investigation.md`. Unlike 6c/6d
this one required actually running an external tool against real RNA-seq
reads, not just querying an API — details below.

## Background

The instruction document's premise, verified rather than assumed:

- `analysis/interaction_scoring.py:139` (`_OPERON_TIGHT_MAX_BP = 150`)
  already solves adjacent-pair operon detection from GFF coordinates
  alone — confirmed directly in the code (the comment at lines 113-138
  documents the same three reference complexes the instruction document
  cites: Mcr cluster, nifI1-nifI2-nifK-nifD, mtpA-mtpC).
- **`claude/genomic_distance_weight_finding.md`, the document the
  instruction file cites as the source of that decision, does not exist
  in this repository** (`git log --all` finds no commit ever adding it,
  on any branch including remotes). This doesn't block the investigation
  — the decision's substance is preserved directly in the code comment
  above — but the user should know the write-up itself appears to be
  lost, if it ever existed as a separate file.
- The real open gap is non-adjacent same-operon members (e.g. NifI1 to
  NifD, with NifI2/NifK in between) — pure pairwise distance scoring
  can't represent "these two are part of one larger transcript" when
  they aren't each other's nearest neighbor.

## Method

1. Confirmed DOOR/MicrobesOnline are unreachable (curl, 15s timeout).
2. Installed a portable JRE (no `apt`/`sudo` available or used), found
   Rockhopper's real download URL (the one in the instruction document
   404s — see below), ran it against a real *M. acetivorans* genome.
3. Built the legacy PTT/RNT genome format Rockhopper's command line
   requires (see finding 2) from the NCBI GenBank flat file via
   Biopython, since NCBI no longer distributes that pre-2016 format
   directly.
4. Pulled real GSE77738 reads directly from ENA (`ftp.sra.ebi.ac.uk`
   over HTTPS) rather than the SRA toolkit the instruction document
   suggested — simpler, no extra tool needed, and produces the same
   FASTQ.
5. Ran Rockhopper at three read-depths from one sample (SRR3158229 /
   "LK2") to both validate against the known operons and get a grounded
   runtime/memory scaling estimate: 100K reads (quick smoke test), 3M
   reads, and 8.35M reads (the largest chunk obtainable before this
   sandbox's network became the limiting factor — see finding 5).

## 1. DOOR / MicrobesOnline: confirmed dead

```
$ curl -sI --max-time 15 http://csbl.bmb.uga.edu/DOOR/     → curl exit 28 (timeout)
$ curl -sI --max-time 15 https://www.microbesonline.org/   → curl exit 60 (SSL cert failure)
```

Matches Cowork's pre-survey exactly — both confirmed dead from this
environment too. No further investigation attempted, per instruction
("5分程度で切り上げてよい"). The pivot to Rockhopper (an actively
maintained, single-organism, run-it-yourself tool) rather than relying on
an external operon database was the right call independent of anything
else in this report.

## 2. Rockhopper: works, with two real gotchas not in the instruction doc

- **Java**: not installed (`command not found`), and this environment has
  no passwordless `sudo`, so `apt install` was not attempted (would have
  needed the user's password/consent, and it's an unnecessary system
  change anyway). Instead: downloaded a portable Eclipse Temurin JRE 17
  tarball (`api.adoptium.net`, ~47 MB), extracted it locally under the
  scratchpad, and ran Rockhopper against that `java` binary directly —
  no root, no system modification, fully reversible. Confirms "Java is
  needed" doesn't have to mean "install Java system-wide."
- **The download URL in the instruction document is stale.**
  `cs.wellesley.edu/~btjaden/Rockhopper/Rockhopper.jar` 404s — the site's
  navigation is JS-rendered (`header.js` → `writeNavbar()`), so a plain
  page fetch (what Cowork's WebFetch pre-survey would have hit too, same
  class of problem seen with BioGRID's downloads page in Phase 6d) can't
  see the real link. Found via `download.html`: the actual path is
  `.../Rockhopper/download/current/Rockhopper.jar` (14 MB, version 2.03).
- **Genome input format**: the command-line tool (not the GUI) requires a
  directory with a `*.fna` (genome sequence), `*.ptt` (protein table),
  and `*.rnt` (RNA table) — **the pre-2016 legacy NCBI format**, not GFF3
  and not the `data/input/genome.gff` this project already has. NCBI
  stopped distributing `.ptt`/`.rnt` files directly years ago (confirmed
  by listing the current `GCF_000007345.1_ASM734v1` FTP directory — no
  such files, only `.gff.gz`/`.gtf.gz`/`.gbff.gz`). Built a ~50-line
  Biopython script converting the GenBank flat file (`.gbff`, which NCBI
  still provides) into `.ptt`/`.rnt` — straightforward, one-time, <1
  second to run, but it is new code a real implementation would need
  (`analysis/operon_bridge.py` or a small preprocessing helper), not
  something obtainable for free from existing project inputs.
  - **Identifier gotcha inside that conversion**: the GenBank file
    carries *two* `old_locus_tag` qualifiers per gene — one without the
    underscore (`MA0001`, the historical PTT-native form) and one with
    it (`MA_0001`, this project's canonical `old_locus_tag` format used
    everywhere else in the pipeline). Biopython's parser returns both as
    a list; naively taking the first one grabs the wrong (no-underscore)
    form. The converter explicitly filters for the underscored variant.
    Rockhopper itself doesn't care which form is in its `Synonym` column,
    but downstream code joining Rockhopper's output back to this
    project's `old_locus_tag` field would silently break without this.
- The Rockhopper GUI can reportedly auto-download genome+annotation data
  from GenBank by organism name (user guide, "if genomic sequence
  information... is available from GenBank, Rockhopper will automatically
  download..."), but the command-line tool's own `-g` help text only
  documents the manual `.fna`/`.ptt`/`.rnt` directory form — not verified
  further since this pipeline needs the headless/scriptable path anyway.

## 3. Input data: Phase 6b's files are not enough, raw reads are on ENA (no SRA toolkit needed)

Confirmed: `.cache/coexpression/` holds only `GSE77738_ReadCounts.xls`,
`GSE64349_TableS1/S2_GEO.xlsx` — processed count/RPKM tables, no raw
reads, matching the instruction document's expectation.

- SRA toolkit (`prefetch`/`fasterq-dump`) is **not** the simplest path.
  ENA mirrors SRA and serves FASTQ directly over HTTPS:
  `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP069835&
  result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,...&
  format=tsv` returns one row per of the 61 runs with a direct
  `ftp.sra.ebi.ac.uk/.../<run>.fastq.gz` path (works over `https://` too)
  — no toolkit install, no two-step prefetch+dump, just `curl`.
- Reads are single-end, 100 bp (`library_layout=SINGLE`), confirming the
  instruction document's "possibly paired-end" was not the case for this
  dataset — simpler than expected on that axis.
- **This sandbox's network is the actual bottleneck, not compute.**
  Sample SRR3158229 (the smallest of the 61 at 1.80 GB gzipped) hit a
  600s curl timeout at 524 MB (≈0.87 MB/s effective). A second sample
  (SRR3158228) stalled at 41 MB and made no further progress. This may be
  an artifact of this specific sandboxed environment rather than a
  general finding — worth the user re-checking actual throughput on
  their own machine before taking the extrapolated 61-sample download
  time below at face value.

## 4. Validation against the known operons — the core question

Ran at three depths from SRR3158229 ("LK2"): 100K / 3M / 8.35M reads (the
last being as much of one sample as this sandbox's network allowed in
the time available — see finding 3). Genome: `NC_003552.1`
(GCF_000007345.1, matches this project's target assembly), PTT/RNT built
per finding 2.

| complex | genes | result (all 3 depths) |
|---|---|---|
| Mcr activation cluster | MA_4546-4550 | ✅ called as **one 5-gene operon**, exact coordinates `5596675-5601614` matching the annotated gene boundaries |
| Nif complex (**including the non-adjacent pair**) | MA_3896(NifI1)-3897(NifI2)-3898-3899, plus upstream `nifH` | ✅ called as **one 5-gene operon**, `4789281-4793839`, spanning `nifH, MA_3896, MA_3897, nifD, nifK` — **NifI1 and the distal gene are correctly grouped into the same transcript despite not being adjacent.** This is the exact gap the instruction document asked Rockhopper to fill, and it filled it, reproducibly at all three read depths (identical operon count — 657 multi-gene operons, 1028 grouped gene-pairs — at 100K, 3M, and 8.35M reads). |
| Mtp pair | MA_4164-MA_4165, 70 bp gap | ❌ **not called as an operon at any depth tested**, despite both genes independently showing stable, comparable, nonzero expression once coverage was adequate (61-64 and 73-76 reads at 3M/8.35M — confirms it wasn't purely a coverage artifact; something about this specific pair or this specific sample/condition didn't cross Rockhopper's boundary-merging threshold). |

**A gene-identity correction, independent of Rockhopper's result**: the
instruction document's table (and the code comment in
`interaction_scoring.py:121`, "nifI1-nifI2-nifK-nifD") label MA_3898 as
NifK and MA_3899 as NifD. Both this project's own
`data/input/genome.gff` and the GenBank record independently say the
opposite — **MA_3898 = nifD, MA_3899 = nifK**. This doesn't affect
Rockhopper's result (which operates on coordinates, not gene-symbol
assumptions) and doesn't affect `_OPERON_TIGHT_MAX_BP`'s correctness (a
generic distance rule, not gene-specific), but the gene-symbol labels in
that source comment and in the (missing) design doc this instruction
carried forward are swapped and worth a one-line doc fix.

**Bottom line on the central question**: Rockhopper delivered exactly the
value the instruction document was hoping for on the flagship example
(NifI1↔NifD/NifK), and did so with a trivial fraction (0.35%) of one
sample's reads — but it is not uniformly reliable: it missed a second,
independently-confirmed real complex (Mtp) with an even tighter gap than
Nif's, on the same sample, at every depth tested. One test sample and one
condition is not enough to know whether Mtp fails to call in general or
only for that sample/condition — a real implementation could not treat
Rockhopper's operon list as ground truth on its own, only as one more
evidence component (same posture as STRING/coexpression already have).

## 5. Cost evaluation

- **Runtime** (this environment, 8 cores, single Rockhopper process per
  run): 100K reads → 25s; 3M reads → 3m26s; 8.35M reads → 5m19s. Scaling
  is clearly sub-linear (large fixed cost — JVM startup, genome indexing,
  and the transcript-boundary/operon-calling stages, which showed
  *identical* output at every depth tested, i.e. they're dominated by
  genome size not read count). Extrapolating from the two largest,
  least-fixed-cost-dominated points (3M→8.35M, ≈21s per additional
  million reads) to a full ≈28.6M-read sample: **roughly 10-15 minutes
  per sample**.
- **Memory**: peak RSS 594 MB (8.35M-read run) to 1.5 GB (3M-read run,
  likely GC-timing noise rather than a real inverse relationship) —
  single-digit GB in the worst case, comfortable even in this sandbox's
  3.7 GB total RAM, not a constraint on any reasonable dev machine.
- **Disk**: per-sample Rockhopper output (`summary.txt` +
  `_transcripts.txt` + `_operons.txt`) is small — **≈356 KB** for the 3M-
  read run (transcript/operon *counts* were identical across all three
  depths, so this size is a reasonable estimate at full scale too, not
  just a subset artifact). Raw FASTQ is the only large footprint (1.8-4.5
  GB gzipped per sample, ≈218 GB gzipped for all 61 GSE77738 samples,
  computed from ENA's `fastq_bytes` field) and does not need to be kept
  after each sample's operons are extracted and cached.
- **Full 61-sample estimate**: compute ≈61 × 10-15 min ≈ **10-15 hours
  serial** (trivially parallelizable across samples given 8 cores and
  ≤600 MB-1.5 GB per process — 4-6 concurrent runs would fit comfortably,
  cutting wall time to roughly 2-4 hours). **Download is the larger cost
  in this environment**: 218 GB at the ≈0.87 MB/s observed here is ≈70
  hours — almost certainly a sandbox-specific network ceiling rather than
  a property of Rockhopper or the dataset itself, but the user should
  verify actual throughput before committing to a 61-sample run.
- **Caching**: yes, realistic, same pattern as STRING (Phase 6a decision
  6) and coexpression (Phase 6b) — run once per sample (or once per
  condition if samples are pooled/collapsed, likely fewer than 61
  distinct biological conditions), keep only the small text outputs
  (`_operons.txt`/`_transcripts.txt`, cacheable via the existing
  `JsonCache`-adjacent pattern after a light parse step), discard the
  large raw FASTQ. No live-query pattern is needed or possible here —
  Rockhopper is a batch tool, not an API.

## 6. Is this worth implementing?

**Cautiously yes — the strongest case of the four Phase 6c/6d/6e
investigations, but with a real reliability caveat the other three
didn't have.**

Unlike KEGG pathways (Phase 6c, rejected: the query genes themselves had
no data) and interolog inference (Phase 6d, rejected: needs a new BLAST
reference genome, and even where data existed it was either empty or
noisy), Rockhopper answered exactly the question this project has an
existing, documented, unfilled gap for — non-adjacent operon members —
using data (GSE77738) and infrastructure (Phase 6b's coexpression
pipeline already normalizes this same dataset) this project already
invested in. It required no new BLAST reference genome and no new
external account/API key (unlike BioGRID). The Java/PTT-RNT friction
points are real but one-time and already solved by this investigation
(portable JRE + a small conversion script), not open problems.

The caveat that keeps this from being an unqualified "implement it now":
one sample, one condition, two real complexes tested — one succeeded
cleanly (Nif, the flagship case), one failed outright (Mtp) despite
adequate coverage. That is not "weak signal, useful someday" like KEGG;
it's "works, but not reliably enough to trust unverified" on the exact
kind of pair this feature exists to catch. A real implementation should
not treat a Rockhopper-predicted operon boundary as ground truth — it
should be one more `EvidenceComponent`, most naturally as a **new
component under the existing `genomic_context` category** (parallel to
how `string_neighborhood` was added alongside `same_gene_neighborhood` in
Phase 6a decision 2), not a replacement for `_gene_neighborhood_v2`, and
not scored as strongly as a same-strand ≤150bp pair already is.

**Recommendation**: worth a scoped follow-up phase, but scope it as
*evidence*, not *ground truth*, and validate on more than one
sample/condition before trusting a "not called as an operon" result as a
negative (Mtp's failure here came from a single sample — a different
growth condition, or pooling multiple GSE77738 samples the way Phase 6b
already does for coexpression, might change that result). Running the
full 61-sample GSE77738 set (or a representative subset across its
distinct growth conditions, not all 61) and cross-checking Mtp and a few
other independently-known pairs across multiple samples before writing
any scoring code would be the natural next step — explicitly out of
scope for this investigation per the instruction document, but the clear
next question its own findings raise.
