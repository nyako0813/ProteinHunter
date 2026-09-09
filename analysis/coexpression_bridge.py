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

See ``claude/coexpression_mutual_rank_normalization.md`` for a second,
later finding along the same lines: even the percentile-rank fix above is
one-sided (ranked only against the *query* gene's background), so a
candidate gene that is simply highly and precisely measured -- a
coexpression "hub" broadly correlated with most of the genome, for reasons
that have nothing to do with the query -- still lands at a high percentile.
A live run against MA_4115 found roughly half of all 351 candidates scoring
above the 0.7 percentile on GSE77738 alone (a properly discriminating
signal should center near 0.5 with most pairs far from the extremes), and a
Pearson correlation of 0.55 between GSE77738's and GSE64349's *independent*
percentiles for the same 291 candidates -- too much cross-dataset agreement
to be query-specific signal, and consistent with a shared "hub" property of
those genes rather than a real relationship to the query. ``normalized_value``
is therefore now each pair's **mutual rank**: the geometric mean of the
correlation's percentile from the query's side and from the candidate's own
side (a gene that is highly correlated with nearly everything has its own
background shifted high too, which deflates exactly this case) -- the same
idea CoNekT/ATTED-II-style coexpression network tools call "Mutual Rank",
adapted here to already-percentile-ranked values. This is a normalization,
not an exclusion: GSE77738 stays in the score (see the module-level note
above about the *other*, earlier finding that led here), it is just no
longer scored one-sided.

See ``claude/coexpression_mutual_rank_normalization.md`` also for a third,
related fix: both loader functions used to check the on-disk per-query
cache only *after* unconditionally checking for/downloading/parsing the
GEO source file, so a query that was already fully cached would still fail
(or silently degrade to an empty bundle) whenever that source file was
missing -- e.g. deleted after first use, or an offline environment. Each
dataset's whole correlation matrix (a single, small, per-dataset object --
not one per candidate, unlike the per-gene backgrounds behind Mutual Rank
above) is now cached once alongside the existing per-query pairs, so a run
whose query genes have all been scored before can be answered entirely
from cache; see ``_try_fully_cached_bundle``.

Data license: GEO is an NCBI/NIH public database with no login, fee, or
reuse restriction; there is no license text comparable to STRING's CC BY 4.0
to enforce programmatically. Crediting the original studies by PMID
(27852217 for GSE77738, 25691524 for GSE64349) is scientific courtesy, not a
license requirement -- see ``output/excel.py``'s Index sheet.
"""

from __future__ import annotations

import gzip
import math
import re
import shutil
from dataclasses import dataclass, field
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
class _OneSidedCoexpression:
    """Internal: one gene's correlation with another, ranked only within the
    FIRST gene's own background distribution. This is the shape actually
    persisted to the on-disk JsonCache (keyed by that first gene) -- see
    CoexpressionPairValue, the public, symmetrized value ``lookup()``
    actually returns, built from two of these (query-side and
    candidate-side) on demand."""

    correlation: float
    percentile: float


@dataclass(slots=True, frozen=True)
class CoexpressionPairValue:
    """One query/candidate pair's coexpression evidence from one GEO dataset."""

    correlation: float
    """Pearson correlation of log2(RPKM+1) across the dataset's retained samples."""

    query_percentile: float
    """This pair's correlation, as a 0.0-1.0 rank within the QUERY gene's own
    background correlation distribution (against every other known gene)."""

    candidate_percentile: float
    """The same correlation value, as a 0.0-1.0 rank within the CANDIDATE
    gene's own background distribution instead. A gene that is simply
    correlated with most of the genome (a coexpression "hub", often just a
    highly and precisely measured gene, for reasons unrelated to the query)
    has its own background shifted high too, so this side deflates exactly
    that case rather than rewarding it."""

    percentile: float
    """The actual score: the geometric mean of ``query_percentile`` and
    ``candidate_percentile`` (a "Mutual Rank", see the module docstring).
    High only when both sides agree this pair is unusually correlated
    relative to their own typical partners -- this is what
    ``normalized_value``/the Excel explanation text use, not either
    one-sided percentile alone."""


@dataclass(slots=True, frozen=True)
class CoexpressionBundle:
    """A parsed, per-query-gene index of coexpression evidence for one GEO dataset."""

    dataset_id: str
    known_tags: frozenset[str]
    pairs_by_query: dict[str, dict[str, _OneSidedCoexpression]]
    warnings: tuple[str, ...]
    n_samples: int = 0
    """Number of retained (post-filtering) samples the correlations were
    computed from -- exposed for explanation text and testing, since this
    dataset's small, hand-curated sample counts are exactly what motivated
    the percentile-based normalization (see the module docstring)."""
    _gene_order: tuple[str, ...] = ()
    _gene_index: dict[str, int] = field(default_factory=dict)
    _z_scores: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    """Retained purely so ``lookup()`` can derive a candidate gene's own
    background distribution on demand (see ``_percentile_from_gene_side``).
    Empty for ``_empty_bundle`` and for datasets that failed to load."""
    _background_cache: dict[str, np.ndarray] = field(default_factory=dict)
    """Per-gene sorted background correlations, memoized in-process only
    (never written to the on-disk JsonCache): the matrix-vector product
    behind it is microseconds, and persisting one full background per
    *candidate* gene -- there can be hundreds per run, unlike the handful of
    configured query genes the on-disk cache is sized for -- would bloat
    that cache file several-hundred-fold for no real benefit."""

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
        one_sided = self.pairs_by_query.get(query_old_locus_tag, {}).get(candidate_old_locus_tag)
        if one_sided is None:
            return None

        candidate_percentile = self._percentile_from_gene_side(candidate_old_locus_tag, one_sided.correlation)
        if candidate_percentile is None:
            # Not expected in practice -- candidate_old_locus_tag was just
            # confirmed correlatable above, from the query's side -- but
            # degrade to the one-sided percentile rather than losing the
            # pair if it ever happens.
            candidate_percentile = one_sided.percentile

        mutual = math.sqrt(max(one_sided.percentile, 0.0) * max(candidate_percentile, 0.0))
        return CoexpressionPairValue(
            correlation=one_sided.correlation,
            query_percentile=one_sided.percentile,
            candidate_percentile=candidate_percentile,
            percentile=mutual,
        )

    def _percentile_from_gene_side(self, gene_tag: str, value: float) -> float | None:
        """Percentile rank of ``value`` within ``gene_tag``'s own background
        correlation distribution (every other known gene, computed from
        ``gene_tag``'s side instead of the query's). Returns None when
        ``gene_tag`` has no usable expression profile (absent from the
        dataset, or zero-variance)."""
        idx = self._gene_index.get(gene_tag)
        if idx is None:
            return None
        sorted_background = self._background_cache.get(gene_tag)
        if sorted_background is None:
            row_correlations = (self._z_scores @ self._z_scores[idx]) / self.n_samples
            mask = np.ones(len(self._gene_order), dtype=bool)
            mask[idx] = False
            sorted_background = np.sort(row_correlations[mask])
            self._background_cache[gene_tag] = sorted_background
        if len(sorted_background) == 0:
            return None
        return float(np.searchsorted(sorted_background, value, side="right")) / len(sorted_background)


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

    # Cache the dataset-level matrix itself (once per dataset, not per
    # query/candidate) so a future run can be answered entirely from cache
    # -- see _try_fully_cached_bundle.
    _cache_matrix(cache, cache_namespace, gene_order, z_scores, n_samples)

    pairs_by_query: dict[str, dict[str, _OneSidedCoexpression]] = {}
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

        pairs: dict[str, _OneSidedCoexpression] = {}
        for gene, idx in gene_index.items():
            if gene == query_tag:
                continue
            r = float(correlations[idx])
            percentile = float(np.searchsorted(sorted_others, r, side="right")) / len(sorted_others)
            pairs[gene] = _OneSidedCoexpression(correlation=r, percentile=percentile)

        cache.set(cache_namespace, query_tag, _encode_pairs_for_cache(pairs))
        pairs_by_query[query_tag] = pairs

    return CoexpressionBundle(
        dataset_id=dataset_id,
        known_tags=frozenset(known_tags),
        pairs_by_query=pairs_by_query,
        warnings=tuple(warnings),
        n_samples=n_samples,
        _gene_order=tuple(gene_order),
        _gene_index=dict(gene_index),
        _z_scores=z_scores,
    )


def _encode_pairs_for_cache(pairs: dict[str, _OneSidedCoexpression]) -> dict[str, dict[str, float]]:
    return {tag: {"correlation": v.correlation, "percentile": v.percentile} for tag, v in pairs.items()}


def _decode_cached_pairs(cached: dict) -> dict[str, _OneSidedCoexpression]:
    decoded: dict[str, _OneSidedCoexpression] = {}
    for tag, value in cached.items():
        if not isinstance(value, dict):
            continue
        try:
            decoded[tag] = _OneSidedCoexpression(
                correlation=float(value.get("correlation", 0.0)),
                percentile=float(value.get("percentile", 0.0)),
            )
        except (TypeError, ValueError):
            continue
    return decoded


def _matrix_cache_namespace(cache_namespace: str) -> str:
    """A dedicated namespace for the dataset-level correlation matrix.

    Kept separate from ``cache_namespace`` (which holds one entry per
    *query* gene) because this holds exactly one entry for the whole
    dataset -- a single ~1-2MB blob for GSE77738/GSE64349's gene counts,
    unlike the per-candidate backgrounds in ``CoexpressionBundle._background_cache``
    which are deliberately NOT persisted (see that field's docstring).
    """
    return f"{cache_namespace}_matrix"


def _cache_matrix(cache: JsonCache, cache_namespace: str, gene_order: list[str], z_scores: np.ndarray, n_samples: int) -> None:
    """Persist the dataset's full correlation matrix so a future run can
    rebuild a bundle -- including candidate-side Mutual Rank percentiles --
    without re-reading the GEO source file at all, as long as every
    requested query gene's pairs are already cached too (see
    ``_try_fully_cached_bundle``)."""
    cache.set(
        _matrix_cache_namespace(cache_namespace),
        "matrix",
        {"gene_order": gene_order, "z_scores": z_scores.tolist(), "n_samples": n_samples},
    )


def _load_cached_matrix(cache: JsonCache, cache_namespace: str) -> tuple[list[str], np.ndarray, int] | None:
    """Load a previously cached correlation matrix, or None if absent/malformed."""
    cached = cache.get(_matrix_cache_namespace(cache_namespace), "matrix")
    if not isinstance(cached, dict):
        return None
    try:
        gene_order = [str(g) for g in cached["gene_order"]]
        z_scores = np.array(cached["z_scores"], dtype=float)
        n_samples = int(cached["n_samples"])
    except (KeyError, TypeError, ValueError):
        return None
    if z_scores.ndim != 2 or z_scores.shape[0] != len(gene_order) or z_scores.shape[1] != n_samples:
        return None
    return gene_order, z_scores, n_samples


def _try_fully_cached_bundle(
    dataset_id: str,
    query_old_locus_tags: list[str],
    cache: JsonCache,
    cache_namespace: str,
) -> CoexpressionBundle | None:
    """Build a complete bundle purely from the on-disk cache, without ever
    touching the GEO source file (no existence check, no download attempt,
    no ``pd.read_excel``).

    Returns None -- meaning "fall through to the normal file-read path" --
    when the dataset's correlation matrix has not been cached yet, or when
    any of the currently requested query genes has not been scored (and
    thus cached) before. This is deliberately conservative: it never
    returns a bundle missing evidence for a query the caller asked about.

    This is the fix for the bug where a query already fully answered by the
    cache would still trigger a source-file existence check/re-download
    attempt first -- and, if that file was unavailable (e.g. deleted, or an
    offline environment), fail or silently degrade to an empty bundle even
    though the cache had everything needed. See
    ``claude/coexpression_mutual_rank_normalization.md``.
    """
    matrix = _load_cached_matrix(cache, cache_namespace)
    if matrix is None:
        return None
    gene_order, z_scores, n_samples = matrix
    gene_index = {gene: i for i, gene in enumerate(gene_order)}

    known_tags: set[str] = set(gene_order)
    pairs_by_query: dict[str, dict[str, _OneSidedCoexpression]] = {}
    for query_tag in dict.fromkeys(tag for tag in query_old_locus_tags if tag):
        cached_pairs = cache.get(cache_namespace, query_tag)
        if not isinstance(cached_pairs, dict):
            # A query never scored (and thus never cached) before -- the
            # matrix alone cannot answer it, so give up on the shortcut
            # entirely rather than return a bundle silently missing this
            # query's evidence.
            return None
        decoded = _decode_cached_pairs(cached_pairs)
        pairs_by_query[query_tag] = decoded
        known_tags.add(query_tag)
        known_tags.update(decoded)

    return CoexpressionBundle(
        dataset_id=dataset_id,
        known_tags=frozenset(known_tags),
        pairs_by_query=pairs_by_query,
        warnings=(),
        n_samples=n_samples,
        _gene_order=tuple(gene_order),
        _gene_index=gene_index,
        _z_scores=z_scores,
    )


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

    fully_cached = _try_fully_cached_bundle(dataset_id, query_old_locus_tags, cache, "coexpression_gse77738")
    if fully_cached is not None:
        return fully_cached

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

    fully_cached = _try_fully_cached_bundle(dataset_id, query_old_locus_tags, cache, "coexpression_gse64349")
    if fully_cached is not None:
        return fully_cached

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
