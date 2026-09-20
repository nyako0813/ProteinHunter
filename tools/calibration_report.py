"""Regenerate the calibration comparison (curated true pairs vs. known negatives) from a results workbook.

An analysis-only tool: it reads three inputs and never touches scoring or the
pipeline's own output code (see claude/calibration_pipeline_expansion_design.md).

Inputs
  --curated    CSV of curated interacting pairs with a ``tier`` column
               (A / B / excluded); schema of claude/experimental_interactions_curated.csv.
  --negatives  CSV of known non-interacting candidates: ``old_locus_tag`` plus
               optional ``source``, ``af3_classification``, ``af3_ipTM`` and
               ``query_old_locus_tag`` (which query the negative was assessed
               against; ``--negatives-query`` supplies it when the column is absent).
  --results    A ProteinHunter results workbook (``04_Score_Breakdown``,
               ``11_Raw_Audit``, ``03_Candidate_Overview``). If
               ``<stem>.run_provenance.yaml`` sits next to it, its fingerprint is
               copied into the summary; otherwise that is skipped.

Outputs (in --out)
  pairs.xlsx (sheet Pairs), negatives_matched.xlsx (sheet Negatives_Matched),
  calibration_summary.md

The two tables are separate workbooks rather than two sheets of one: they have
different columns and are read independently (a reviewer opens one, a script
loads one), the file names stay the ones the design names, and the overwrite
guard stays a simple per-file check. Scores are written as numbers and missing
values as empty cells, so the tables filter and sort in Excel as they are.

Statistics use only numpy/pandas: Mann-Whitney U (exact when the smaller sample has
at most 8 values and there are no ties, otherwise the normal approximation with tie and continuity
correction -- the same rule as scipy.stats.mannwhitneyu's default) and
ROC-AUC = U / (n1 * n2).

Usage
  python tools/calibration_report.py --curated claude/experimental_interactions_curated.csv \
      --negatives claude/calibration/af3_negatives_MA_4115.csv \
      --results data/output/ProteinHunter_results_calibration_check.xlsx \
      --out claude/calibration/2026-09-19_run
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

SCORE_COLUMNS: tuple[str, ...] = (
    "interaction_score",
    "interaction_priority_score",
    "string_neighborhood",
    "string_cooccurrence",
    "coexpression_gse77738",
    "coexpression_gse64349",
    "final_score",
)

#: Components read from 11_Raw_Audit (normalized_value) and attached to every matched row.
COMPONENT_COLUMNS: tuple[str, ...] = (
    "source_classification",
    "sequence_evidence",
    "genomic_context",
    "string_neighborhood",
    "co_occurrence",
    "domain_complementarity",
    "string_cooccurrence",
    "coexpression_gse77738",
    "coexpression_gse64349",
)

#: Components whose AVAILABLE/MISSING status is kept as ``<name>_status`` (external data may be absent).
STATUS_COLUMNS: tuple[str, ...] = (
    "string_neighborhood",
    "string_cooccurrence",
    "coexpression_gse77738",
    "coexpression_gse64349",
)

VALID_TIERS = ("A", "B", "excluded")
EXACT_MAX_N = 8  # exact p-value when the smaller sample has at most this many values (and no ties)
EXACT_MAX_PRODUCT = 5000  # ...and n1 * n2 stays small enough for the exact counting to be cheap
SIDECAR_SUFFIX = ".run_provenance.yaml"

PAIRS_FILE = "pairs.xlsx"
NEGATIVES_FILE = "negatives_matched.xlsx"
SUMMARY_FILE = "calibration_summary.md"
PAIRS_SHEET = "Pairs"
NEGATIVES_SHEET = "Negatives_Matched"

#: Columns written as numbers (blank -> empty cell) rather than text.
NUMERIC_COLUMNS: tuple[str, ...] = (
    "candidate_rank",
    "final_score",
    "interaction_score",
    "interaction_priority_score",
    *COMPONENT_COLUMNS,
)
MAX_COLUMN_WIDTH = 45

PAIR_RESULT_COLUMNS: tuple[str, ...] = (
    "candidate_source",
    "candidate_sources",
    "candidate_rank",
    "negative_hit_strength",
    "final_score",
    "final_score_tier",
    "interaction_score",
    "interaction_priority_score",
    "interaction_evidence_tier",
    *COMPONENT_COLUMNS,
    *(f"{name}_status" for name in STATUS_COLUMNS),
)


class CalibrationInputError(Exception):
    """A user-facing problem with an input file."""


# ---------------------------------------------------------------------------
# Statistics (numpy/pandas only)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MannWhitney:
    n1: int
    n2: int
    u1: float  # U statistic of the first sample (ties count 1/2)
    auc: float  # u1 / (n1 * n2); 0.5 = no separation, 1.0 = first sample always higher
    p_value: float  # two-sided
    method: str  # "exact" or "asymptotic"


@lru_cache(maxsize=None)
def _u_counts(n1: int, n2: int) -> tuple[int, ...]:
    """Number of rank arrangements giving each U value 0..n1*n2 (no ties)."""
    if n1 == 0 or n2 == 0:
        return (1,)
    # c(n1, n2, u) = c(n1 - 1, n2, u - n2) + c(n1, n2 - 1, u)
    first = _u_counts(n1 - 1, n2)
    second = _u_counts(n1, n2 - 1)
    size = n1 * n2 + 1
    counts = [0] * size
    for u in range(size):
        if 0 <= u - n2 < len(first):
            counts[u] += first[u - n2]
        if u < len(second):
            counts[u] += second[u]
    return tuple(counts)


def mann_whitney(first: Any, second: Any) -> MannWhitney | None:
    """Two-sided Mann-Whitney U test of ``first`` vs ``second``; None when either sample is empty."""
    x = pd.to_numeric(pd.Series(first), errors="coerce").dropna().to_numpy(dtype=float)
    y = pd.to_numeric(pd.Series(second), errors="coerce").dropna().to_numpy(dtype=float)
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return None

    combined = pd.Series(np.concatenate([x, y]))
    ranks = combined.rank(method="average").to_numpy()
    u1 = float(ranks[:n1].sum() - n1 * (n1 + 1) / 2)
    u2 = n1 * n2 - u1
    auc = u1 / (n1 * n2)

    _, tie_counts = np.unique(combined.to_numpy(), return_counts=True)
    has_ties = bool((tie_counts > 1).any())
    u_max = max(u1, u2)

    if min(n1, n2) <= EXACT_MAX_N and n1 * n2 <= EXACT_MAX_PRODUCT and not has_ties:
        counts = _u_counts(n1, n2)
        total = sum(counts)
        tail = sum(counts[int(round(u_max)):])
        return MannWhitney(n1, n2, u1, auc, min(1.0, 2.0 * tail / total), "exact")

    n = n1 + n2
    tie_term = float(((tie_counts**3) - tie_counts).sum()) / (n * (n - 1))
    variance = n1 * n2 / 12.0 * ((n + 1) - tie_term)
    if variance <= 0:
        return MannWhitney(n1, n2, u1, auc, float("nan"), "asymptotic")
    z = (u_max - n1 * n2 / 2.0 - 0.5) / math.sqrt(variance)  # continuity-corrected
    p_value = min(1.0, 2.0 * 0.5 * math.erfc(z / math.sqrt(2.0)))
    return MannWhitney(n1, n2, u1, auc, p_value, "asymptotic")


# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------


def _read_csv(path: Path, what: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    except FileNotFoundError as exc:
        raise CalibrationInputError(f"{what} file not found: {path}") from exc


def read_curated(path: Path) -> pd.DataFrame:
    frame = _read_csv(path, "curated pairs")
    required = {"protein_a_old_locus_tag", "protein_b_old_locus_tag", "tier"}
    missing = required - set(frame.columns)
    if missing:
        hint = " (add a `tier` column with A / B / excluded)" if "tier" in missing else ""
        raise CalibrationInputError(f"{path}: missing column(s) {sorted(missing)}{hint}")
    frame = frame.copy()
    frame["tier"] = frame["tier"].str.strip().map(lambda value: "excluded" if value.lower() == "excluded" else value.upper())
    bad = sorted(set(frame["tier"]) - set(VALID_TIERS))
    if bad:
        raise CalibrationInputError(f"{path}: tier must be one of {VALID_TIERS}, found {bad}")
    for column in ("protein_a_label", "protein_b_label", "class", "confidence", "source"):
        if column not in frame.columns:
            frame[column] = ""
    return frame


def read_negatives(path: Path, default_query: str) -> pd.DataFrame:
    frame = _read_csv(path, "negatives")
    if "old_locus_tag" not in frame.columns:
        raise CalibrationInputError(f"{path}: missing column 'old_locus_tag'")
    frame = frame.copy()
    for column, default in (("source", ""), ("af3_classification", ""), ("af3_ipTM", "")):
        if column not in frame.columns:
            frame[column] = default
    if "query_old_locus_tag" not in frame.columns:
        frame["query_old_locus_tag"] = default_query
    else:
        frame["query_old_locus_tag"] = frame["query_old_locus_tag"].replace("", default_query)
    return frame


def _read_sheet(path: Path, sheet: str, required: set[str]) -> pd.DataFrame | None:
    """Read a workbook sheet whose header is on row 1 or 2 (row 1 holds the 'Back to Index' link)."""
    try:
        workbook = pd.ExcelFile(path)
    except FileNotFoundError as exc:
        raise CalibrationInputError(f"results workbook not found: {path}") from exc
    if sheet not in workbook.sheet_names:
        return None
    for header in (1, 0, 2):
        frame = workbook.parse(sheet, header=header)
        if required <= set(map(str, frame.columns)):
            return frame
    return None


@dataclass(frozen=True)
class Results:
    pairs: pd.DataFrame  # one row per (query, candidate, candidate_source) with scores and components
    tag_by_protein: dict[str, str]


def read_results(path: Path) -> Results:
    breakdown = _read_sheet(path, "04_Score_Breakdown", {"query_id", "candidate_protein_id", "interaction_score"})
    if breakdown is None:
        raise CalibrationInputError(
            f"{path}: sheet '04_Score_Breakdown' with query_id / candidate_protein_id / interaction_score "
            "was not found (is this a ProteinHunter results workbook?)"
        )
    audit = _read_sheet(path, "11_Raw_Audit", {"query_id", "candidate_protein_id", "component_name", "status"})
    overview = _read_sheet(path, "03_Candidate_Overview", {"protein_id", "old_locus_tag"})

    tag_by_protein: dict[str, str] = {}
    if overview is not None:
        for protein_id, tag in zip(overview["protein_id"], overview["old_locus_tag"]):
            if isinstance(tag, str) and tag.strip():
                tag_by_protein[str(protein_id)] = tag.strip()
    if audit is not None:
        for id_column, tag_column in (("query_id", "query_old_locus_tag"), ("candidate_protein_id", "candidate_old_locus_tag")):
            if tag_column in audit.columns:
                for protein_id, tag in zip(audit[id_column], audit[tag_column]):
                    if isinstance(tag, str) and tag.strip() and str(protein_id) not in tag_by_protein:
                        tag_by_protein[str(protein_id)] = tag.strip()

    pairs = breakdown.copy()
    pairs["query_id"] = pairs["query_id"].astype(str)
    pairs["candidate_protein_id"] = pairs["candidate_protein_id"].astype(str)
    pairs["query_tag"] = pairs["query_id"].map(tag_by_protein)
    candidate_tag = pairs["candidate_old_locus_tag"] if "candidate_old_locus_tag" in pairs.columns else pd.Series(index=pairs.index, dtype=object)
    pairs["candidate_tag"] = candidate_tag.where(candidate_tag.astype(str).str.strip().ne("") & candidate_tag.notna())
    pairs["candidate_tag"] = pairs["candidate_tag"].fillna(pairs["candidate_protein_id"].map(tag_by_protein))

    for column in ("candidate_rank", "final_score", "interaction_score", "interaction_priority_score"):
        if column in pairs.columns:
            pairs[column] = pd.to_numeric(pairs[column], errors="coerce")
        else:
            pairs[column] = np.nan
    for column in ("candidate_source", "negative_hit_strength", "final_score_tier", "interaction_evidence_tier"):
        if column not in pairs.columns:
            pairs[column] = ""

    components = _pivot_components(audit)
    pairs = pairs.merge(components, how="left", on=["query_id", "candidate_protein_id", "candidate_source"])
    for column in COMPONENT_COLUMNS:
        if column not in pairs.columns:
            pairs[column] = np.nan
    for column in STATUS_COLUMNS:
        if f"{column}_status" not in pairs.columns:
            pairs[f"{column}_status"] = ""
    return Results(pairs=pairs, tag_by_protein=tag_by_protein)


def _pivot_components(audit: pd.DataFrame | None) -> pd.DataFrame:
    """Wide component table from the long 11_Raw_Audit sheet (which also holds other stacked tables)."""
    keys = ["query_id", "candidate_protein_id", "candidate_source"]
    if audit is None or not set(keys + ["component_name", "status", "normalized_value"]) <= set(audit.columns):
        return pd.DataFrame(columns=keys)
    rows = audit[audit["component_name"].isin(COMPONENT_COLUMNS) & audit["status"].isin(["AVAILABLE", "MISSING", "NOT_APPLICABLE"])].copy()
    rows["query_id"] = rows["query_id"].astype(str)
    rows["candidate_protein_id"] = rows["candidate_protein_id"].astype(str)
    rows["value"] = pd.to_numeric(rows["normalized_value"], errors="coerce").where(rows["status"].eq("AVAILABLE"))
    values = rows.pivot_table(index=keys, columns="component_name", values="value", aggfunc="first", dropna=False)
    frame = values.reset_index()
    status = rows.pivot_table(index=keys, columns="component_name", values="status", aggfunc="first")
    status = status[[column for column in STATUS_COLUMNS if column in status.columns]]
    status.columns = [f"{column}_status" for column in status.columns]
    return frame.merge(status.reset_index(), how="left", on=keys)


def read_sidecar(results_path: Path) -> dict[str, Any] | None:
    """The run-provenance sidecar next to the results workbook, or None (absent or unreadable)."""
    sidecar = results_path.with_name(f"{results_path.stem}{SIDECAR_SUFFIX}")
    if not sidecar.is_file():
        return None
    try:
        payload = yaml.safe_load(sidecar.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError):
        return None
    return payload if isinstance(payload, dict) and isinstance(payload.get("run_provenance"), dict) else None


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------


def _best_row(rows: pd.DataFrame) -> pd.Series:
    """One row per (query, candidate): highest final_score, then lowest candidate_rank."""
    ordered = rows.sort_values(["final_score", "candidate_rank"], ascending=[False, True], na_position="last", kind="stable")
    best = ordered.iloc[0].copy()
    best["candidate_sources"] = ";".join(sorted({str(value) for value in rows["candidate_source"] if str(value)}))
    return best


def _result_fields(row: pd.Series | None) -> dict[str, Any]:
    if row is None:
        return {column: "" for column in PAIR_RESULT_COLUMNS}
    return {column: row.get(column, "") for column in PAIR_RESULT_COLUMNS}


def _lookup(results: Results, query: str, candidate: str) -> pd.Series | None:
    pairs = results.pairs
    hits = pairs[(pairs["query_tag"] == query) & (pairs["candidate_tag"] == candidate)]
    return None if hits.empty else _best_row(hits)


def match_curated(curated: pd.DataFrame, results: Results) -> pd.DataFrame:
    queries_present = set(results.pairs["query_tag"].dropna())
    rows: list[dict[str, Any]] = []
    for _, item in curated.iterrows():
        a, b = item["protein_a_old_locus_tag"].strip(), item["protein_b_old_locus_tag"].strip()
        base = {
            "curated_source": item["source"],
            "curated_class": item["class"],
            "curated_confidence": item["confidence"],
            "tier": item["tier"],
        }
        directions = [(a, item["protein_a_label"], b, item["protein_b_label"]), (b, item["protein_b_label"], a, item["protein_a_label"])]
        emitted = False
        for query, query_label, candidate, candidate_label in directions:
            if item["tier"] == "excluded" or "/" in candidate or query not in queries_present:
                continue
            row = _lookup(results, query, candidate)
            rows.append(
                {**base, "query": query, "query_label": query_label, "candidate": candidate, "candidate_label": candidate_label,
                 "status": "matched" if row is not None else "candidate_not_in_results", **_result_fields(row)}
            )
            emitted = True
        if not emitted:
            reason = "excluded" if item["tier"] == "excluded" else (
                "unresolved_multi_locus" if "/" in b or "/" in a else "no_query_in_results"
            )
            rows.append(
                {**base, "query": a, "query_label": item["protein_a_label"], "candidate": b, "candidate_label": item["protein_b_label"],
                 "status": reason, **_result_fields(None)}
            )
    return pd.DataFrame(rows)


def match_negatives(negatives: pd.DataFrame, results: Results) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, item in negatives.iterrows():
        query, candidate = item["query_old_locus_tag"].strip(), item["old_locus_tag"].strip()
        row = _lookup(results, query, candidate)
        rows.append(
            {
                "old_locus_tag": candidate,
                "query_old_locus_tag": query,
                "source": item["source"],
                "af3_classification": item["af3_classification"],
                "af3_ipTM": item["af3_ipTM"],
                "found": row is not None,
                **_result_fields(row),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(frame[column], errors="coerce").dropna()


def _fmt(value: float | None, digits: int = 2) -> str:
    return "n/a" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:.{digits}f}"


def _fmt_p(value: float) -> str:
    if math.isnan(value):
        return "n/a"
    return "<0.001" if value < 0.001 else f"{value:.3f}"


def render_summary(
    *,
    curated_path: Path,
    negatives_path: Path,
    results_path: Path,
    matched_pairs: pd.DataFrame,
    matched_negatives: pd.DataFrame,
    sidecar: dict[str, Any] | None,
    generated_at: datetime,
) -> str:
    matched = matched_pairs[matched_pairs["status"] == "matched"]
    groups = {"Tier A": matched[matched["tier"] == "A"], "Tier B": matched[matched["tier"] == "B"], "Negatives": matched_negatives[matched_negatives["found"].astype(bool)]}

    lines = ["# Calibration summary", ""]
    lines += [
        f"Generated: {generated_at:%Y-%m-%d %H:%M} by `tools/calibration_report.py`. "
        "Analysis only: nothing here changes scoring.",
        "",
        "## Inputs",
        "",
        f"- curated pairs: `{curated_path}`",
        f"- negatives: `{negatives_path}`",
        f"- results workbook: `{results_path}`",
    ]
    if sidecar is None:
        lines.append(f"- run provenance: none (no `{results_path.stem}{SIDECAR_SUFFIX}` next to the workbook)")
    else:
        run = sidecar["run_provenance"]
        dirty = "+dirty" if run.get("git_dirty") else ""
        check = run.get("reference_genome_check_passed")
        lines += [
            f"- run provenance (`{results_path.stem}{SIDECAR_SUFFIX}`):",
            f"  - config hash: `{run.get('config_hash', 'unknown')}`",
            f"  - code: {run.get('app_version', '?')} (git {run.get('git_commit') or 'unknown'}{dirty})",
            f"  - reference genome check: {'passed' if check is True else 'FAILED - see that run log' if check is False else 'not verified'}",
            f"  - generated: {run.get('generated_at', 'unknown')}",
        ]
        settings = _effective_settings(sidecar.get("effective_config"))
        if settings:
            lines.append(f"  - settings: {settings}")
    lines.append("")

    lines += ["## Match coverage", ""]
    tier_status = matched_pairs.groupby(["tier", "status"]).size()
    lines += ["| tier | status | curated rows/directions |", "|---|---|---:|"]
    for (tier, status), count in tier_status.items():
        lines.append(f"| {tier} | {status} | {count} |")
    n_neg = len(matched_negatives)
    lines += ["", f"Negatives found in the results: {int(matched_negatives['found'].astype(bool).sum())} of {n_neg}.", ""]

    lines += ["## Scores by group (mean / median)", ""]
    header = "| score | " + " | ".join(f"{name} (n, mean / median)" for name in groups) + " |"
    lines += [header, "|---|" + "---|" * len(groups)]
    for column in SCORE_COLUMNS:
        cells = []
        for frame in groups.values():
            values = _numeric(frame, column)
            cells.append(f"{len(values)}, {_fmt(values.mean() if len(values) else None)} / {_fmt(values.median() if len(values) else None)}")
        lines.append(f"| `{column}` | " + " | ".join(cells) + " |")
    lines.append("")

    lines += [
        "## Separation from the negatives",
        "",
        "AUC = probability that a random positive scores above a random negative "
        "(0.5 = no separation, ties count half); p is the two-sided Mann-Whitney U test.",
        "",
        "| score | Tier A vs neg: AUC | p | Tier B vs neg: AUC | p |",
        "|---|---:|---:|---:|---:|",
    ]
    for column in SCORE_COLUMNS:
        cells = []
        for tier in ("Tier A", "Tier B"):
            result = mann_whitney(_numeric(groups[tier], column), _numeric(groups["Negatives"], column))
            cells += ["n/a", "n/a"] if result is None else [_fmt(result.auc, 3), f"{_fmt_p(result.p_value)} ({result.method})"]
        lines.append(f"| `{column}` | " + " | ".join(cells) + " |")
    lines += [
        "",
        "Read with care: the samples are small, the negatives come from one query, and "
        "the positives from several, so a low p or a high AUC is a prompt to look, not a calibration fit.",
        "",
    ]
    return "\n".join(lines)


def _effective_settings(config: Any) -> str:
    """A one-line digest of the settings that most change what the numbers mean."""
    if not isinstance(config, dict):
        return ""
    scoring = config.get("interaction_scoring") or {}
    parts = []
    for label, value in (
        ("scoring_model", scoring.get("scoring_model")),
        ("ranking_metric", scoring.get("ranking_metric")),
        ("max_candidates_per_query", scoring.get("max_candidates_per_query")),
        ("geo_coexpression", scoring.get("geo_coexpression_enabled")),
        ("consider_cross_species_matches", config.get("consider_cross_species_matches")),
    ):
        if value is not None:
            parts.append(f"{label}={value}")
    sources = scoring.get("candidate_sources")
    if isinstance(sources, dict):
        parts.append("candidate_sources=" + ",".join(sorted(name for name, on in sources.items() if on)))
    return "; ".join(parts)


# ---------------------------------------------------------------------------
# Excel output
# ---------------------------------------------------------------------------


def _prepare_for_excel(frame: pd.DataFrame) -> pd.DataFrame:
    """Numbers as numbers, blanks as empty cells (the matching code uses "" for "no value")."""
    prepared = frame.copy()
    for column in NUMERIC_COLUMNS:
        if column in prepared.columns:
            prepared[column] = pd.to_numeric(prepared[column], errors="coerce")
    for column in prepared.columns:
        if prepared[column].dtype == object:
            prepared[column] = prepared[column].map(lambda value: None if isinstance(value, str) and value == "" else value)
    return prepared


def write_excel_table(frame: pd.DataFrame, path: Path, sheet_name: str) -> None:
    """Write one table as a single-sheet workbook: bold header, frozen header row, filter, readable widths."""
    prepared = _prepare_for_excel(frame)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        prepared.to_excel(writer, sheet_name=sheet_name, index=False)
        worksheet = writer.sheets[sheet_name]
        worksheet.freeze_panes = "A2"
        if len(prepared.columns) > 0:
            last_column = get_column_letter(len(prepared.columns))
            worksheet.auto_filter.ref = f"A1:{last_column}{max(1, len(prepared) + 1)}"
        for cell in worksheet[1]:
            cell.font = Font(bold=True)
        for index, column in enumerate(prepared.columns, start=1):
            if column in NUMERIC_COLUMNS or prepared[column].dtype.kind in "biuf":
                width = max(len(str(column)), 10)
            else:
                texts = [len(str(value)) for value in prepared[column].head(200) if value is not None and not pd.isna(value)]
                width = max([len(str(column)), *texts])
            worksheet.column_dimensions[get_column_letter(index)].width = min(width + 2, MAX_COLUMN_WIDTH)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def run(
    curated: Path,
    negatives: Path,
    results: Path,
    out: Path,
    *,
    negatives_query: str = "MA_4115",
    force: bool = False,
    generated_at: datetime | None = None,
) -> dict[str, Path]:
    curated_frame = read_curated(curated)
    negatives_frame = read_negatives(negatives, negatives_query)
    result_tables = read_results(results)

    targets = {name: out / name for name in (PAIRS_FILE, NEGATIVES_FILE, SUMMARY_FILE)}
    existing = [name for name, path in targets.items() if path.exists()]
    if existing and not force:
        raise CalibrationInputError(
            f"{out} already contains {', '.join(existing)}; use a new (e.g. dated) --out directory "
            "to keep the history, or --force to overwrite"
        )

    matched_pairs = match_curated(curated_frame, result_tables)
    matched_negatives = match_negatives(negatives_frame, result_tables)
    summary = render_summary(
        curated_path=curated,
        negatives_path=negatives,
        results_path=results,
        matched_pairs=matched_pairs,
        matched_negatives=matched_negatives,
        sidecar=read_sidecar(results),
        generated_at=generated_at or datetime.now(),
    )

    out.mkdir(parents=True, exist_ok=True)
    write_excel_table(matched_pairs, targets[PAIRS_FILE], PAIRS_SHEET)
    write_excel_table(matched_negatives, targets[NEGATIVES_FILE], NEGATIVES_SHEET)
    targets[SUMMARY_FILE].write_text(summary, encoding="utf-8")
    return targets


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Regenerate the calibration comparison from a ProteinHunter results workbook.")
    parser.add_argument("--curated", type=Path, required=True, help="curated pairs CSV with a tier column (A/B/excluded)")
    parser.add_argument("--negatives", type=Path, required=True, help="known-negative candidates CSV")
    parser.add_argument("--results", type=Path, required=True, help="ProteinHunter results workbook (.xlsx)")
    parser.add_argument("--out", type=Path, required=True, help="output directory (use a new, dated one to keep history)")
    parser.add_argument("--negatives-query", default="MA_4115", help="query old_locus_tag for negatives without a query_old_locus_tag column")
    parser.add_argument("--force", action="store_true", help="overwrite existing outputs in --out")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    try:
        targets = run(args.curated, args.negatives, args.results, args.out, negatives_query=args.negatives_query, force=args.force)
    except CalibrationInputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    for path in targets.values():
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
