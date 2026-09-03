"""Optional bridge to public GEO coexpression evidence for *M. acetivorans*.

Reads two published, processed RNA-seq supplementary files from NCBI GEO
(GSE77738, GSE64349) -- never anything from this pipeline's own scoring --
and turns them into per-query-gene coexpression evidence. See
``claude/phase6b_coexpression_design.md`` for the investigation this module
is built from: in particular, why GSE77738's 61 GEO-listed samples are
**not** 61 independent replicates (most are an actinomycin-D RNA-decay time
course; only 13 are true steady-state samples), why GSE64349's Delta-msrH
mutant subset is excluded while its "WWM82 (parental strain)" subset is kept
as additional wild-type replicates, why gene IDs need only a trivial
transform (not a real mapping project), and why a fixed linear map from
Pearson r to ``normalized_value`` is not used (GSE64349's small sample count
badly inflates its background gene-pair correlation, so ``normalized_value``
is instead each pair's percentile rank within its own query gene's
background correlation distribution).

Data license: GEO is an NCBI/NIH public database with no login, fee, or
reuse restriction; there is no license text comparable to STRING's CC BY 4.0
to enforce programmatically. Crediting the original studies by PMID
(27852217 for GSE77738, 25691524 for GSE64349) is scientific courtesy, not a
license requirement -- see ``output/excel.py``'s Index sheet.
"""

from __future__ import annotations

import gzip
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from core.cache import JsonCache
from core.exceptions import CoexpressionAnnotationError

GEO_FTP_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"

GSE77738_PMID = "27852217"
GSE64349_PMID = "25691524"

#: GSE77738 is primarily an actinomycin-D RNA-decay time course (cells
#: sampled at 0/5/10/20/30/60/120/240 min after halting transcription, to
#: measure RNA half-life), not 61 independent condition replicates. Using
#: every sample in a correlation matrix would mix real coexpression signal
#: with shared decay-kinetics correlation (nearly all transcripts decline
#: together after actinomycin D). These 13 RPKM-sheet column names are the
#: ones NOT part of that decay chase: 9 explicit "time point: 0 min"
#: replicates (3 methanol, 3 trimethylamine, 3 acetate) plus 4 more samples
#: with no time-point annotation at all (2 methanol, 2 trimethylamine).
#: Verified directly against GSE77738_series_matrix.txt's
#: growth-media/time-point sample characteristics -- see
#: claude/phase6b_coexpression_design.md. Hardcoded (not re-derived from the
#: series matrix at runtime) because GSE77738 is a fixed, already-published
#: dataset that will not change, the same reasoning
#: analysis/string_ppi_bridge.py uses to hardcode STRING's strain taxid.
#: Ambiguous cases -- e.g. "Metcalf_C2AM1_R1.PF.fastq" is GEO sample title
#: "C2AM1" (the trailing "_R1" is a lane/replicate suffix in the filename,
#: not part of the title) vs. "Metcalf2_C2AT_R1.PF.fastq" is GEO sample
#: title "C2AT_R1" (there the trailing "_R1" IS part of the title) -- cannot
#: be resolved by a general regex, which is why exact column names are
#: listed rather than a pattern.
GSE77738_STEADY_STATE_RPKM_COLUMNS: frozenset[str] = frozenset(
    {
        "LK1_ATCACG_L007_R1_001.fastq",
        "LK9_TTAGGC_L003_R1_001.fastq",
        "LK17_GGCTAC_L004_R1_001.fastq",
        "LK25_ATCACG_L003_R1_001.fastq",
        "LK31_CAGATC_L004_R1_001.fastq",
        "LK37_ATCACG_L005_R1_001.fastq",
        "LK43_CAGATC_L006_R1_001.fastq",
        "LK49_ATCACG_L007_R1_001.fastq",
        "LK55_CAGATC_L008_R1_001.fastq",
        "Metcalf_C2AM1_R1.PF.fastq",
        "Metcalf_C2AM3_R1.PF.fastq",
        "Metcalf2_C2AT_R1.PF.fastq",
        "Metcalf2_C2AT_R2.PF.fastq",
    }
)

#: GSE64349's TableS2 bundles a *different* comparison (Delta-msrH deletion
#: mutant vs. its "WWM82 (parental strain)" control) into the same GEO
#: series as TableS1's wild-type DMS/MMPA/MeSH/MeOH comparison. The mutant
#: subset is excluded for the same reason GSE66445 (a metabolically
#: engineered strain) was excluded in Phase 6a's scope; "WWM82 (parental
#: strain)" is genetically wild-type and is kept as additional replicates,
#: since GSE64349 is already sample-starved (see phase6b design doc).
_GSE64349_MUTANT_COLUMN_PREFIX = "delta-msrH -"
_GSE64349_RPKM_COLUMN_SUFFIX = " - RPKM"

_MA_LOCUS_PATTERN = re.compile(r"^MA(\d{4})$")


def _to_old_locus_tag(gene_locus: object) -> str | None:
    """Convert a bare GEO gene-locus id ('MA0001') to this project's 'MA_0001' form.

    Returns None for anything that is not a 4-digit protein-coding locus
    (e.g. 'MAt4684' tRNA features, or a missing/blank cell).
    """
    match = _MA_LOCUS_PATTERN.match(str(gene_locus).strip())
    if match is None:
        return None
    return f"MA_{match.group(1)}"


@dataclass(slots=True, frozen=True)
class CoexpressionPairValue:
    """One query/candidate pair's coexpression evidence from one GEO dataset."""

    correlation: float
    """Pearson correlation of log2(RPKM+1) across the dataset's retained samples."""

    percentile: float
    """This pair's correlation, as a 0.0-1.0 rank within the query gene's own
    background correlation distribution (against every other known gene in
    this dataset). See the module docstring for why a fixed linear map of
    ``correlation`` is not used instead."""


@dataclass(slots=True, frozen=True)
class CoexpressionBundle:
    """A parsed, per-query-gene index of coexpression evidence for one GEO dataset."""

    dataset_id: str
    known_tags: frozenset[str]
    pairs_by_query: dict[str, dict[str, CoexpressionPairValue]]
    warnings: tuple[str, ...]
    n_samples: int = 0
    """Number of retained (post-filtering) samples the correlations were
    computed from -- exposed for explanation text and testing, since this
    dataset's small, hand-curated sample counts are exactly what motivated
    the percentile-based normalization (see the module docstring)."""

    def lookup(self, query_old_locus_tag: str, candidate_old_locus_tag: str) -> CoexpressionPairValue | None:
        """Return coexpression evidence for one pair, or None when unavailable.

        Unlike STRING's sparse links file (Phase 6a), this dataset's
        correlation is *dense* -- once a query gene's expression is known,
        its correlation with every other known gene can be computed. So
        there is no separate "known to the dataset but this specific pair
        was never evaluated" case the way STRING has one: None here means
        either the query or the candidate gene is absent from this
        dataset's gene list entirely (MISSING), or one of the two has
        zero-variance expression across the retained samples (correlation
        is mathematically undefined for it, treated the same as MISSING).
        """
        if not query_old_locus_tag or not candidate_old_locus_tag:
            return None
        if query_old_locus_tag not in self.known_tags or candidate_old_locus_tag not in self.known_tags:
            return None
        return self.pairs_by_query.get(query_old_locus_tag, {}).get(candidate_old_locus_tag)


def _empty_bundle(dataset_id: str) -> CoexpressionBundle:
    return CoexpressionBundle(dataset_id=dataset_id, known_tags=frozenset(), pairs_by_query={}, warnings=())


def _download_gz_and_decompress(url: str, destination: Path, timeout: int) -> None:
    """Download a .gz file and leave its decompressed contents at ``destination``."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp_gz = destination.with_suffix(destination.suffix + ".gz.part")
    try:
        with requests.get(url, timeout=timeout, stream=True) as response:
            response.raise_for_status()
            with tmp_gz.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    handle.write(chunk)
        tmp_final = destination.with_suffix(destination.suffix + ".part")
        with gzip.open(tmp_gz, "rb") as src, tmp_final.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        tmp_final.replace(destination)
    except (requests.RequestException, OSError) as exc:
        raise CoexpressionAnnotationError(
            f"Could not download/decompress GEO supplementary file from '{url}'. "
            "Please check the network connection."
        ) from exc
    finally:
        tmp_gz.unlink(missing_ok=True)


def _build_symbol_to_locus_table(gse77738_readcounts_path: Path) -> dict[str, str]:
    """Build a gene-symbol -> old_locus_tag lookup from GSE77738's own Gene Name column.

    GSE64349's 'Feature ID' column mixes gene symbols ('cdc6_1', 'repA')
    with bare locus ids for genes that have no common name -- this table
    resolves the symbol rows without any external resource, using data
    already downloaded for the GSE77738 bundle. See the module docstring.
    """
    sheet = pd.read_excel(gse77738_readcounts_path, sheet_name="RPKM Normalized Read Counts")
    table: dict[str, str] = {}
    for gene_locus, gene_name in zip(sheet["Gene Locus"], sheet["Gene Name"]):
        locus = _to_old_locus_tag(gene_locus)
        name = str(gene_name).strip()
        if locus and name and name != "-":
            table.setdefault(name.lower(), locus)
    return table


def _feature_id_to_locus(feature_id: object, symbol_table: dict[str, str]) -> str | None:
    """Resolve one GSE64349 'Feature ID' cell to an old_locus_tag, or None."""
    direct = _to_old_locus_tag(feature_id)
    if direct:
        return direct
    key = str(feature_id).strip().lower()
    if key in symbol_table:
        return symbol_table[key]
    # A trailing "_<digit>" disambiguates duplicate symbols (e.g. 'cdc6_1');
    # GSE77738's own Gene Name column never carries this suffix.
    stripped = re.sub(r"_\d+$", "", key)
    return symbol_table.get(stripped)


def _log2_matrix(expression: pd.DataFrame) -> pd.DataFrame:
    """log2(RPKM + 1), the standard transform for count-derived expression data."""
    return np.log2(expression.astype(float) + 1.0)


def _build_bundle_from_matrix(
    dataset_id: str,
    log_matrix: pd.DataFrame,
    query_old_locus_tags: list[str],
    cache: JsonCache,
    cache_namespace: str,
    warnings: list[str],
) -> CoexpressionBundle:
    """Shared correlation/percentile/caching logic for both GEO datasets.

    ``log_matrix``: rows indexed by old_locus_tag (already deduplicated),
    columns are the dataset's retained samples, values are log2(RPKM+1).
    """
    known_tags: set[str] = set(log_matrix.index)

    gene_variance = log_matrix.var(axis=1)
    zero_variance_genes = set(gene_variance.index[gene_variance == 0.0])
    if zero_variance_genes:
        warnings.append(
            f"{dataset_id}: {len(zero_variance_genes)} gene(s) have zero-variance "
            "expression across the retained samples and cannot be correlated "
            "with anything; treated as unavailable for those specific pairs."
        )

    usable = log_matrix.drop(index=zero_variance_genes) if zero_variance_genes else log_matrix
    gene_order = list(usable.index)
    values = usable.to_numpy(dtype=float)
    means = values.mean(axis=1, keepdims=True)
    stds = values.std(axis=1, keepdims=True)
    z_scores = (values - means) / stds
    n_samples = values.shape[1]
    gene_index = {gene: i for i, gene in enumerate(gene_order)}

    pairs_by_query: dict[str, dict[str, CoexpressionPairValue]] = {}
    for query_tag in dict.fromkeys(tag for tag in query_old_locus_tags if tag):
        cached = cache.get(cache_namespace, query_tag)
        if isinstance(cached, dict):
            pairs_by_query[query_tag] = _decode_cached_pairs(cached)
            known_tags.add(query_tag)
            known_tags.update(pairs_by_query[query_tag])
            continue

        if query_tag not in gene_index:
            # Absent from the dataset entirely (MISSING), or one of the
            # zero-variance genes (cannot be correlated) -- either way,
            # nothing to compute or cache for this query.
            continue

        query_row = z_scores[gene_index[query_tag]]
        correlations = (z_scores @ query_row) / n_samples
        self_idx = gene_index[query_tag]
        other_mask = np.ones(len(gene_order), dtype=bool)
        other_mask[self_idx] = False
        other_correlations = correlations[other_mask]
        sorted_others = np.sort(other_correlations)

        pairs: dict[str, CoexpressionPairValue] = {}
        for gene, idx in gene_index.items():
            if gene == query_tag:
                continue
            r = float(correlations[idx])
            percentile = float(np.searchsorted(sorted_others, r, side="right")) / len(sorted_others)
            pairs[gene] = CoexpressionPairValue(correlation=r, percentile=percentile)

        cache.set(cache_namespace, query_tag, _encode_pairs_for_cache(pairs))
        pairs_by_query[query_tag] = pairs

    return CoexpressionBundle(
        dataset_id=dataset_id,
        known_tags=frozenset(known_tags),
        pairs_by_query=pairs_by_query,
        warnings=tuple(warnings),
        n_samples=n_samples,
    )


def _encode_pairs_for_cache(pairs: dict[str, CoexpressionPairValue]) -> dict[str, dict[str, float]]:
    return {tag: {"correlation": v.correlation, "percentile": v.percentile} for tag, v in pairs.items()}


def _decode_cached_pairs(cached: dict) -> dict[str, CoexpressionPairValue]:
    decoded: dict[str, CoexpressionPairValue] = {}
    for tag, value in cached.items():
        if not isinstance(value, dict):
            continue
        try:
            decoded[tag] = CoexpressionPairValue(
                correlation=float(value.get("correlation", 0.0)),
                percentile=float(value.get("percentile", 0.0)),
            )
        except (TypeError, ValueError):
            continue
    return decoded


def load_gse77738_coexpression_bundle(
    enabled: bool,
    query_old_locus_tags: list[str],
    cache: JsonCache,
    cache_dir: Path,
) -> CoexpressionBundle:
    """Load GSE77738 (acetate/methanol/TMA growth) coexpression evidence.

    Returns an empty bundle (every lookup returns None) when ``enabled`` is
    False. Never raises: a download/parse failure with nothing cached
    degrades to an empty bundle with a warning, the same "optional,
    best-effort evidence" behavior as the STRING and PIH bridges.
    """
    dataset_id = "gse77738"
    if not enabled:
        return _empty_bundle(dataset_id)

    warnings: list[str] = []
    readcounts_path = Path(cache_dir) / "coexpression" / "GSE77738_ReadCounts.xls"
    try:
        if not readcounts_path.exists():
            _download_gz_and_decompress(
                f"{GEO_FTP_BASE}/GSE77nnn/GSE77738/suppl/GSE77738_ReadCounts.xls.gz",
                readcounts_path,
                timeout=120,
            )
        rpkm_sheet = pd.read_excel(readcounts_path, sheet_name="RPKM Normalized Read Counts")
    except (CoexpressionAnnotationError, OSError, ValueError) as exc:
        warnings.append(f"gse77738: could not obtain/parse GSE77738_ReadCounts.xls: {exc}")
        return CoexpressionBundle(dataset_id=dataset_id, known_tags=frozenset(), pairs_by_query={}, warnings=tuple(warnings))

    steady_state_columns = [c for c in rpkm_sheet.columns if c in GSE77738_STEADY_STATE_RPKM_COLUMNS]
    missing_columns = GSE77738_STEADY_STATE_RPKM_COLUMNS - set(rpkm_sheet.columns)
    if missing_columns:
        warnings.append(
            f"gse77738: expected steady-state sample column(s) not found in the "
            f"downloaded file, continuing with what is present: {sorted(missing_columns)}"
        )

    rpkm_sheet = rpkm_sheet.copy()
    rpkm_sheet["old_locus_tag"] = rpkm_sheet["Gene Locus"].map(_to_old_locus_tag)
    expression = (
        rpkm_sheet.dropna(subset=["old_locus_tag"])
        .drop_duplicates(subset=["old_locus_tag"])
        .set_index("old_locus_tag")[steady_state_columns]
    )
    log_matrix = _log2_matrix(expression)

    return _build_bundle_from_matrix(
        dataset_id, log_matrix, query_old_locus_tags, cache, "coexpression_gse77738", warnings
    )


def load_gse64349_coexpression_bundle(
    enabled: bool,
    query_old_locus_tags: list[str],
    cache: JsonCache,
    cache_dir: Path,
) -> CoexpressionBundle:
    """Load GSE64349 (methylated sulfur compounds) coexpression evidence.

    Pools TableS1 (9 wild-type DMS/MMPA/MeSH/MeOH samples) with TableS2's
    "WWM82 (parental strain)" subset (3 more wild-type replicates); TableS2's
    Delta-msrH mutant subset is excluded. See the module docstring and
    claude/phase6b_coexpression_design.md.

    Returns an empty bundle (every lookup returns None) when ``enabled`` is
    False. Never raises: a download/parse failure with nothing cached
    degrades to an empty bundle with a warning, matching the other bridges.
    """
    dataset_id = "gse64349"
    if not enabled:
        return _empty_bundle(dataset_id)

    warnings: list[str] = []
    coexpression_dir = Path(cache_dir) / "coexpression"
    table1_path = coexpression_dir / "GSE64349_TableS1_GEO.xlsx"
    table2_path = coexpression_dir / "GSE64349_TableS2_GEO.xlsx"
    # The symbol->locus table is built from GSE77738's own data (see
    # _build_symbol_to_locus_table); this bundle downloads that file too if
    # it is not already present, so GSE64349 can be enabled independently of
    # GSE77738. In the common case both are enabled together and this reuses
    # whatever load_gse77738_coexpression_bundle already fetched.
    readcounts_path = coexpression_dir / "GSE77738_ReadCounts.xls"

    try:
        if not table1_path.exists():
            _download_gz_and_decompress(
                f"{GEO_FTP_BASE}/GSE64nnn/GSE64349/suppl/GSE64349_TableS1_GEO.xlsx.gz",
                table1_path,
                timeout=180,
            )
        if not table2_path.exists():
            _download_gz_and_decompress(
                f"{GEO_FTP_BASE}/GSE64nnn/GSE64349/suppl/GSE64349_TableS2_GEO.xlsx.gz",
                table2_path,
                timeout=120,
            )
        if not readcounts_path.exists():
            _download_gz_and_decompress(
                f"{GEO_FTP_BASE}/GSE77nnn/GSE77738/suppl/GSE77738_ReadCounts.xls.gz",
                readcounts_path,
                timeout=120,
            )
        symbol_table = _build_symbol_to_locus_table(readcounts_path)
        table1 = pd.read_excel(table1_path)
        table2 = pd.read_excel(table2_path)
    except (CoexpressionAnnotationError, OSError, ValueError) as exc:
        warnings.append(f"gse64349: could not obtain/parse GEO supplementary files: {exc}")
        return CoexpressionBundle(dataset_id=dataset_id, known_tags=frozenset(), pairs_by_query={}, warnings=tuple(warnings))

    table1_rpkm_cols = [c for c in table1.columns if c.endswith(_GSE64349_RPKM_COLUMN_SUFFIX)]
    table2_wt_cols = [
        c
        for c in table2.columns
        if c.endswith(_GSE64349_RPKM_COLUMN_SUFFIX) and not c.startswith(_GSE64349_MUTANT_COLUMN_PREFIX)
    ]

    table1 = table1.copy()
    table1["old_locus_tag"] = table1["Feature ID"].map(lambda f: _feature_id_to_locus(f, symbol_table))
    table2 = table2.copy()
    table2["old_locus_tag"] = table2["Feature ID"].map(lambda f: _feature_id_to_locus(f, symbol_table))

    expr1 = (
        table1.dropna(subset=["old_locus_tag"]).drop_duplicates(subset=["old_locus_tag"]).set_index("old_locus_tag")[
            table1_rpkm_cols
        ]
    )
    expr2 = (
        table2.dropna(subset=["old_locus_tag"]).drop_duplicates(subset=["old_locus_tag"]).set_index("old_locus_tag")[
            table2_wt_cols
        ]
    )
    combined = expr1.join(expr2, how="outer")
    log_matrix = _log2_matrix(combined)

    return _build_bundle_from_matrix(
        dataset_id, log_matrix, query_old_locus_tags, cache, "coexpression_gse64349", warnings
    )


__all__: tuple[str, ...] = (
    "GSE64349_PMID",
    "GSE77738_PMID",
    "CoexpressionAnnotationError",
    "CoexpressionBundle",
    "CoexpressionPairValue",
    "load_gse64349_coexpression_bundle",
    "load_gse77738_coexpression_bundle",
)
