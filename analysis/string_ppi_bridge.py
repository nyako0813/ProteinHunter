"""Optional bridge to STRING (string-db.org) known/predicted PPI evidence.

STRING is a separate, independent, publicly funded database; this module
only reads its published bulk-download files (and, as a fallback, its
public REST API) -- never anything from this pipeline's own scoring, and
STRING is never imported as code. See
``claude/phase6_external_evidence_design.md`` for the investigation this
module is built from: in particular, why the strain-level NCBI taxid must
be used (species-level taxids return nothing), why ``old_locus_tag`` (not
the RefSeq protein_id) is the correct join key for this organism, why the
bulk-download files are used instead of live per-pair API calls (STRING's
own guidance: download the full dataset for anything beyond occasional
access), and why only the ``cooccurrence``/``neighborhood`` channels are
read here (``fusion`` evidence measured 0% for this organism;
``experiments``/``databases``/``textmining``/``coexpression`` were all
either too sparse, likely homology-transferred rather than direct, or
noisy for this project's target organism).

Data license: STRING data is CC BY 4.0 -- see
https://string-db.org/cgi/access?footer_active_subpage=licensing. Any
Excel output that includes STRING-derived evidence must credit STRING
(see ``output/excel.py``'s Index sheet).
"""

from __future__ import annotations

import gzip
import time
from dataclasses import dataclass
from pathlib import Path

import requests

from core.cache import JsonCache
from core.exceptions import StringPpiAnnotationError

STRING_API_BASE = "https://string-db.org/api"
STRING_DOWNLOAD_BASE = "https://stringdb-downloads.org/download"
STRING_VERSION = "12.0"
STRING_CALLER_IDENTITY = "proteinhunter_v5"
#: STRING's own courtesy guidance: wait at least this long between live API calls.
STRING_MIN_CALL_INTERVAL_SECONDS = 1.0

#: Bulk-file column order after protein1/protein2 (space-separated), and
#: the field names this module uses for them everywhere. STRING spells
#: "cooccurence" with one 'r' in its own files; normalized to
#: "cooccurrence" (correct spelling) as soon as it is parsed.
_BULK_FILE_CHANNELS: tuple[str, ...] = (
    "neighborhood",
    "fusion",
    "cooccurrence",
    "coexpression",
    "experimental",
    "database",
    "textmining",
)

#: Live API JSON field prefixes (nscore/fscore/pscore/ascore/escore/dscore/tscore)
#: mapped to the same channel names used above, so both sources produce
#: identical StringPairScores regardless of which one supplied the data.
_LIVE_API_CHANNEL_PREFIXES: dict[str, str] = {
    "n": "neighborhood",
    "f": "fusion",
    "p": "cooccurrence",
    "a": "coexpression",
    "e": "experimental",
    "d": "database",
    "t": "textmining",
}


@dataclass(slots=True, frozen=True)
class StringPairScores:
    """One query/candidate pair's STRING evidence, on STRING's native 0-1000 scale."""

    neighborhood: int = 0
    fusion: int = 0
    cooccurrence: int = 0
    coexpression: int = 0
    experimental: int = 0
    database: int = 0
    textmining: int = 0
    combined_score: int = 0


_ZERO_SCORES = StringPairScores()


@dataclass(slots=True, frozen=True)
class StringPpiBundle:
    """A parsed, per-query-protein index of STRING evidence for one species."""

    ncbi_taxon_id: int
    known_tags: frozenset[str]
    pairs_by_query: dict[str, dict[str, StringPairScores]]
    warnings: tuple[str, ...]

    def lookup(self, query_old_locus_tag: str, candidate_old_locus_tag: str) -> StringPairScores | None:
        """Return STRING evidence for one pair, or None when either side is unknown to STRING.

        A pair where both proteins are known to STRING but the pair itself
        never appears in STRING's data returns an all-zero
        ``StringPairScores`` (evaluated, no signal) rather than None --
        STRING's cooccurrence/neighborhood methods are computed
        systematically across the whole proteome, so absence from the
        sparse bulk file is treated as "computed, scored zero," not "never
        evaluated." See claude/phase6_external_evidence_design.md decision 5.
        """
        if not query_old_locus_tag or not candidate_old_locus_tag:
            return None
        if query_old_locus_tag not in self.known_tags or candidate_old_locus_tag not in self.known_tags:
            return None
        return self.pairs_by_query.get(query_old_locus_tag, {}).get(candidate_old_locus_tag, _ZERO_SCORES)


def load_string_ppi_bundle(
    ncbi_taxon_id: int | None,
    query_old_locus_tags: list[str],
    cache: JsonCache,
    cache_dir: Path,
) -> StringPpiBundle:
    """Load STRING evidence for the given query proteins, from cache/local bulk files/live API.

    Returns an empty bundle (every lookup returns None) when
    ``ncbi_taxon_id`` is None -- STRING evidence is opt-in via
    ``interaction_scoring.string_ppi_ncbi_taxon_id``. Never raises: a
    download/parse failure with nothing cached degrades to an empty bundle
    with a warning, the same "optional, best-effort evidence" behavior as
    the PIH bridge (analysis/pih_evidence_bridge.py) and consistent with
    the design specification's rule that an unavailable external source
    must not stop a local run.
    """
    if ncbi_taxon_id is None:
        return StringPpiBundle(ncbi_taxon_id=0, known_tags=frozenset(), pairs_by_query={}, warnings=())

    warnings: list[str] = []
    string_dir = Path(cache_dir) / "string_ppi_network"

    try:
        known_tags: set[str] = set(_load_known_tags(ncbi_taxon_id, string_dir, warnings))
    except StringPpiAnnotationError as exc:
        warnings.append(str(exc))
        known_tags = set()

    pairs_by_query: dict[str, dict[str, StringPairScores]] = {}
    last_live_call = 0.0
    for query_tag in dict.fromkeys(tag for tag in query_old_locus_tags if tag):
        cache_key = f"{ncbi_taxon_id}:{query_tag}"
        cached = cache.get("string_ppi", cache_key)
        if isinstance(cached, dict):
            decoded = _decode_cached_pairs(cached)
            pairs_by_query[query_tag] = decoded
            # A cached result proves this query (and whatever partners it
            # has) were resolved before, via whichever path (bulk scan or
            # live fallback) originally populated it -- restore that
            # "known" status here too, since protein.info alone cannot
            # reconstruct it when the live-fallback path was used
            # (see the live-fallback branch below).
            known_tags.add(query_tag)
            known_tags.update(decoded)
            continue

        try:
            pairs = _scan_bulk_file_for_query(ncbi_taxon_id, query_tag, string_dir, warnings)
        except StringPpiAnnotationError as exc:
            warnings.append(str(exc))
            pairs = None

        if pairs is None:
            # No local bulk file for this species -- fall back to a single
            # live call. protein.info (the normal source of known_tags)
            # was never fetched either in this branch, so a successful
            # live response is used as a stand-in "this protein is known
            # to STRING" signal for the query and every partner it
            # returned; candidates outside that already-filtered live
            # result stay indistinguishable from MISSING here. This is a
            # deliberate simplification for the secondary fallback path --
            # see claude/phase6_external_evidence_design.md.
            last_live_call = _throttle(last_live_call)
            try:
                pairs = _fetch_string_partners_live(ncbi_taxon_id, query_tag)
                known_tags.add(query_tag)
                known_tags.update(pairs)
            except StringPpiAnnotationError as exc:
                warnings.append(str(exc))
                pairs = {}

        cache.set("string_ppi", cache_key, _encode_pairs_for_cache(pairs))
        pairs_by_query[query_tag] = pairs

    return StringPpiBundle(
        ncbi_taxon_id=ncbi_taxon_id,
        known_tags=frozenset(known_tags),
        pairs_by_query=pairs_by_query,
        warnings=tuple(warnings),
    )


def _throttle(last_call: float) -> float:
    """Sleep if needed so live STRING calls stay at least STRING_MIN_CALL_INTERVAL_SECONDS apart."""
    elapsed = time.monotonic() - last_call
    if last_call > 0.0 and elapsed < STRING_MIN_CALL_INTERVAL_SECONDS:
        time.sleep(STRING_MIN_CALL_INTERVAL_SECONDS - elapsed)
    return time.monotonic()


def _bulk_file_paths(ncbi_taxon_id: int, string_dir: Path) -> tuple[Path, Path]:
    links_path = string_dir / f"{ncbi_taxon_id}.protein.links.detailed.v{STRING_VERSION}.txt.gz"
    info_path = string_dir / f"{ncbi_taxon_id}.protein.info.v{STRING_VERSION}.txt.gz"
    return links_path, info_path


def _download_bulk_file(url: str, destination: Path, timeout: int) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = destination.with_suffix(destination.suffix + ".part")
    try:
        with requests.get(url, timeout=timeout, stream=True) as response:
            response.raise_for_status()
            with tmp_path.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    handle.write(chunk)
    except requests.RequestException as exc:
        tmp_path.unlink(missing_ok=True)
        raise StringPpiAnnotationError(
            f"Could not download STRING bulk file from '{url}'. Please check the network connection."
        ) from exc
    tmp_path.replace(destination)


def _load_known_tags(ncbi_taxon_id: int, string_dir: Path, warnings: list[str]) -> frozenset[str]:
    """Return every old_locus_tag STRING has data for in this species (protein.info file)."""
    _links_path, info_path = _bulk_file_paths(ncbi_taxon_id, string_dir)
    if not info_path.exists():
        url = (
            f"{STRING_DOWNLOAD_BASE}/protein.info.v{STRING_VERSION}/"
            f"{ncbi_taxon_id}.protein.info.v{STRING_VERSION}.txt.gz"
        )
        _download_bulk_file(url, info_path, timeout=60)

    tags: set[str] = set()
    try:
        with gzip.open(info_path, "rt", encoding="utf-8") as handle:
            next(handle, None)  # header
            for line in handle:
                columns = line.rstrip("\n").split("\t")
                if not columns or not columns[0]:
                    continue
                tags.add(_strip_taxon_prefix(columns[0]))
    except OSError as exc:
        raise StringPpiAnnotationError(
            f"Could not read the cached STRING protein.info file '{info_path}': {exc}"
        ) from exc
    return frozenset(tags)


def _scan_bulk_file_for_query(
    ncbi_taxon_id: int, query_tag: str, string_dir: Path, warnings: list[str]
) -> dict[str, StringPairScores] | None:
    """Scan the local bulk links file for one query's partners.

    Returns None (triggering the live-API fallback) only when the bulk
    file itself could not be obtained -- an empty dict is a legitimate
    "found the file, this query has zero partners in it" result.
    """
    links_path, _info_path = _bulk_file_paths(ncbi_taxon_id, string_dir)
    if not links_path.exists():
        url = (
            f"{STRING_DOWNLOAD_BASE}/protein.links.detailed.v{STRING_VERSION}/"
            f"{ncbi_taxon_id}.protein.links.detailed.v{STRING_VERSION}.txt.gz"
        )
        try:
            _download_bulk_file(url, links_path, timeout=180)
        except StringPpiAnnotationError:
            return None

    query_id = f"{ncbi_taxon_id}.{query_tag}"
    pairs: dict[str, StringPairScores] = {}
    try:
        with gzip.open(links_path, "rt", encoding="utf-8") as handle:
            next(handle, None)  # header
            for line in handle:
                columns = line.rstrip("\n").split(" ")
                if len(columns) != 2 + len(_BULK_FILE_CHANNELS) + 1:
                    continue
                protein1, protein2 = columns[0], columns[1]
                if protein1 != query_id and protein2 != query_id:
                    continue
                other = protein2 if protein1 == query_id else protein1
                other_tag = _strip_taxon_prefix(other)
                channel_values = [_safe_int(value) for value in columns[2:-1]]
                pairs[other_tag] = StringPairScores(
                    **dict(zip(_BULK_FILE_CHANNELS, channel_values)),
                    combined_score=_safe_int(columns[-1]),
                )
    except OSError as exc:
        warnings.append(f"Could not read the cached STRING links file '{links_path}': {exc}")
        return None
    return pairs


def _fetch_string_partners_live(ncbi_taxon_id: int, query_tag: str) -> dict[str, StringPairScores]:
    """Fallback for a species with no local bulk file yet: one live interaction_partners call."""
    params = {
        "identifiers": f"{ncbi_taxon_id}.{query_tag}",
        "species": ncbi_taxon_id,
        "caller_identity": STRING_CALLER_IDENTITY,
        "required_score": 1,
        "limit": 5000,
    }
    try:
        response = requests.get(f"{STRING_API_BASE}/json/interaction_partners", params=params, timeout=30)
        response.raise_for_status()
        payload = response.json()
    except requests.RequestException as exc:
        raise StringPpiAnnotationError(
            f"STRING live API query failed for '{query_tag}'. Please check the network connection."
        ) from exc
    except ValueError as exc:
        raise StringPpiAnnotationError(f"STRING returned invalid JSON for '{query_tag}'.") from exc

    if not isinstance(payload, list):
        raise StringPpiAnnotationError(f"STRING returned an unexpected response for '{query_tag}'.")

    query_id = f"{ncbi_taxon_id}.{query_tag}"
    pairs: dict[str, StringPairScores] = {}
    for entry in payload:
        if not isinstance(entry, dict):
            continue
        id_a, id_b = entry.get("stringId_A"), entry.get("stringId_B")
        if id_a != query_id and id_b != query_id:
            continue
        other = id_b if id_a == query_id else id_a
        if not isinstance(other, str):
            continue
        other_tag = _strip_taxon_prefix(other)
        channel_values = {
            channel: round(_safe_float(entry.get(f"{prefix}score")) * 1000)
            for prefix, channel in _LIVE_API_CHANNEL_PREFIXES.items()
        }
        pairs[other_tag] = StringPairScores(
            **channel_values, combined_score=round(_safe_float(entry.get("score")) * 1000)
        )
    return pairs


def _strip_taxon_prefix(string_id: str) -> str:
    """Return the old_locus_tag part of a STRING id ('188937.MA_4115' -> 'MA_4115')."""
    return string_id.split(".", 1)[1] if "." in string_id else string_id


def _safe_int(value: object) -> int:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0


def _safe_float(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return 0.0


def _encode_pairs_for_cache(pairs: dict[str, StringPairScores]) -> dict[str, dict[str, int]]:
    return {
        tag: {
            "neighborhood": scores.neighborhood,
            "fusion": scores.fusion,
            "cooccurrence": scores.cooccurrence,
            "coexpression": scores.coexpression,
            "experimental": scores.experimental,
            "database": scores.database,
            "textmining": scores.textmining,
            "combined_score": scores.combined_score,
        }
        for tag, scores in pairs.items()
    }


def _decode_cached_pairs(cached: dict) -> dict[str, StringPairScores]:
    decoded: dict[str, StringPairScores] = {}
    for tag, value in cached.items():
        if not isinstance(value, dict):
            continue
        decoded[tag] = StringPairScores(
            neighborhood=_safe_int(value.get("neighborhood", 0)),
            fusion=_safe_int(value.get("fusion", 0)),
            cooccurrence=_safe_int(value.get("cooccurrence", 0)),
            coexpression=_safe_int(value.get("coexpression", 0)),
            experimental=_safe_int(value.get("experimental", 0)),
            database=_safe_int(value.get("database", 0)),
            textmining=_safe_int(value.get("textmining", 0)),
            combined_score=_safe_int(value.get("combined_score", 0)),
        )
    return decoded


__all__: tuple[str, ...] = (
    "STRING_CALLER_IDENTITY",
    "STRING_VERSION",
    "StringPairScores",
    "StringPpiAnnotationError",
    "StringPpiBundle",
    "load_string_ppi_bundle",
)
