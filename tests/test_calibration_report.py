"""Tests for tools/calibration_report.py."""

from __future__ import annotations

import csv
import math
from datetime import datetime
from pathlib import Path

import pandas as pd
import pytest
import yaml

from tools import calibration_report as cr

REPO_ROOT = Path(__file__).resolve().parent.parent
GENERATED_AT = datetime(2026, 9, 20, 12, 0)

# ---------------------------------------------------------------------------
# Mann-Whitney U / AUC
# ---------------------------------------------------------------------------
# Expected p-values below were produced with scipy.stats.mannwhitneyu (default
# method), which is NOT a dependency of this project; only the U-derived
# quantities are re-derived by hand in the comments.


def test_mann_whitney_perfect_separation_exact() -> None:
    # [3,4,5] vs [1,2]: every positive beats every negative -> U = 6 of 3*2, AUC = 1.
    # Exact two-sided p = 2 * (1 / C(5,2)) = 0.2.
    result = cr.mann_whitney([3, 4, 5], [1, 2])

    assert (result.n1, result.n2, result.u1, result.auc) == (3, 2, 6.0, 1.0)
    assert result.method == "exact"
    assert result.p_value == pytest.approx(0.2)


def test_mann_whitney_with_ties_counts_half_and_uses_normal_approximation() -> None:
    # [1,2,2,5] vs [2,3,3,4,6]: pairwise wins 0+.5+.5+... = U 5.0 of 20 -> AUC 0.25.
    result = cr.mann_whitney([1, 2, 2, 5], [2, 3, 3, 4, 6])

    assert result.u1 == 5.0
    assert result.auc == pytest.approx(0.25)
    assert result.method == "asymptotic"
    assert result.p_value == pytest.approx(0.26017490098354834, rel=1e-9)


def test_mann_whitney_larger_samples_use_normal_approximation() -> None:
    x = [0.9, 0.8, 0.75, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3]
    y = [0.7, 0.65, 0.52, 0.48, 0.42, 0.38, 0.33, 0.28, 0.25, 0.2, 0.15, 0.1]

    result = cr.mann_whitney(x, y)

    assert result.u1 == 91.0
    assert result.auc == pytest.approx(91 / 120)
    assert result.method == "asymptotic"
    assert result.p_value == pytest.approx(0.04431379228316605, rel=1e-9)


def test_auc_is_u_over_n1_n2_and_flips_with_the_order_of_the_samples() -> None:
    a = [0.2, 0.9, 0.5, 0.7]
    b = [0.1, 0.4, 0.8]

    forward, backward = cr.mann_whitney(a, b), cr.mann_whitney(b, a)

    assert forward.auc == pytest.approx(forward.u1 / (4 * 3))
    assert forward.auc + backward.auc == pytest.approx(1.0)
    assert forward.p_value == pytest.approx(backward.p_value)


def test_mann_whitney_edge_cases() -> None:
    assert cr.mann_whitney([], [1.0]) is None
    assert cr.mann_whitney([float("nan")], [1.0]) is None  # missing values are dropped
    all_tied = cr.mann_whitney([1, 1, 1], [1, 1])
    assert all_tied.auc == 0.5 and math.isnan(all_tied.p_value)


# ---------------------------------------------------------------------------
# Synthetic workbook
# ---------------------------------------------------------------------------

COMPONENTS = {
    "string_neighborhood": None,
    "string_cooccurrence": None,
    "coexpression_gse77738": None,
    "coexpression_gse64349": None,
    "genomic_context": None,
}


def pair_row(query: str, candidate: str, source: str, interaction: float, priority: float, final: float, rank: int = 1, **components: float | None) -> dict:
    return {
        "query": query,
        "candidate": candidate,
        "source": source,
        "interaction_score": interaction,
        "interaction_priority_score": priority,
        "final_score": final,
        "rank": rank,
        "components": {**COMPONENTS, **components},
    }


def write_results_workbook(path: Path, rows: list[dict], *, back_link_row: bool = True) -> None:
    """A minimal results workbook in the real layout: header on row 2 (row 1 = 'Back to Index')."""
    tags = sorted({r["query"] for r in rows} | {r["candidate"] for r in rows})
    overview = pd.DataFrame({"protein_id": [f"WP_{t}" for t in tags], "old_locus_tag": tags})
    breakdown = pd.DataFrame(
        [
            {
                "query_id": f"WP_{r['query']}",
                "query_protein_id": f"WP_{r['query']}",
                "candidate_rank": r["rank"],
                "candidate_protein_id": f"WP_{r['candidate']}",
                "candidate_old_locus_tag": r["candidate"],
                "candidate_source": r["source"],
                "negative_hit_strength": "none",
                "interaction_priority_score": r["interaction_priority_score"],
                "interaction_score": r["interaction_score"],
                "final_score": r["final_score"],
                "interaction_evidence_tier": "Tier3_Moderate",
                "final_score_tier": "Tier3_Moderate",
            }
            for r in rows
        ]
    )
    audit_rows = []
    for r in rows:
        for name, value in r["components"].items():
            audit_rows.append(
                {
                    "query_id": f"WP_{r['query']}",
                    "query_old_locus_tag": r["query"],
                    "candidate_protein_id": f"WP_{r['candidate']}",
                    "candidate_old_locus_tag": r["candidate"],
                    "candidate_source": r["source"],
                    "component_name": name,
                    "status": "AVAILABLE" if value is not None else "MISSING",
                    "raw_value": value,
                    "normalized_value": value,
                }
            )
    audit = pd.DataFrame(audit_rows)
    # The real 11_Raw_Audit stacks other tables underneath; they must be ignored.
    junk = pd.DataFrame([{"query_id": "resolved from target records", "candidate_protein_id": 1, "component_name": "notes", "status": "x"}])
    startrow = 1 if back_link_row else 0
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        for name, frame in (("03_Candidate_Overview", overview), ("04_Score_Breakdown", breakdown), ("11_Raw_Audit", pd.concat([audit, junk], ignore_index=True))):
            frame.to_excel(writer, sheet_name=name, index=False, startrow=startrow)
            if back_link_row:
                writer.sheets[name]["A1"] = "Back to Index"


def write_csv(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        csv.writer(handle).writerows([header, *rows])
    return path


CURATED_HEADER = ["protein_a_old_locus_tag", "protein_a_label", "protein_b_old_locus_tag", "protein_b_label", "class", "confidence", "source", "tier"]


@pytest.fixture
def small_case(tmp_path: Path) -> dict[str, Path]:
    """Two queries (Q1, Q2), Tier A/B pairs, three negatives for query Q1. Values chosen for hand calculation."""
    rows = [
        # Tier A: Q1->P1 (both directions exist for P2<->Q2) ...
        pair_row("MA_Q1", "MA_P1", "Candidates", 60, 40, 70, string_neighborhood=0.9, coexpression_gse77738=0.8),
        pair_row("MA_Q1", "MA_P2", "Candidates", 40, 20, 50, string_neighborhood=0.5, coexpression_gse77738=0.6),
        pair_row("MA_Q2", "MA_Q1", "Candidates", 50, 30, 60, string_neighborhood=0.7),
        # Tier B
        pair_row("MA_Q1", "MA_P3", "Candidates", 20, 10, 30, string_neighborhood=0.2),
        # negatives (query MA_Q1)
        pair_row("MA_Q1", "MA_N1", "Candidates", 10, 15, 25, string_neighborhood=0.0, coexpression_gse77738=0.9),
        pair_row("MA_Q1", "MA_N2", "Candidates", 30, 25, 45, string_neighborhood=0.1, coexpression_gse77738=0.7),
        # a candidate found in two buckets: best final_score wins, both sources are listed
        pair_row("MA_Q1", "MA_N3", "Candidates_relaxed", 5, 5, 10, string_neighborhood=0.0),
        pair_row("MA_Q1", "MA_N3", "Negative_hit", 20, 5, 35, string_neighborhood=0.3),
    ]
    results = tmp_path / "results.xlsx"
    write_results_workbook(results, rows)
    curated = write_csv(
        tmp_path / "curated.csv",
        CURATED_HEADER,
        [
            ["MA_Q1", "Q1", "MA_P1", "P1", "strict", "ok", "paper 1", "A"],
            ["MA_Q1", "Q1", "MA_P2", "P2", "strict", "ok", "paper 1", "A"],
            ["MA_Q1", "Q1", "MA_Q2", "Q2", "strict", "ok", "paper 2", "A"],  # both directions: Q1->Q2 not in results as candidate, Q2->Q1 is
            ["MA_Q1", "Q1", "MA_P3", "P3", "soft", "check", "paper 3", "B"],
            ["MA_Q1", "Q1", "MA_X9", "X9", "soft", "check", "paper 4", "B"],  # candidate absent from results
            ["MA_Z1", "Z1", "MA_Z2", "Z2", "strict", "-", "paper 5", "B"],  # neither protein is a query
            ["MA_Q1", "Q1", "MA_P1/MA_P2", "multi", "soft", "-", "paper 6", "B"],  # unresolved multi-locus partner
            ["MA_Q1", "Q1", "MA_P9", "P9", "strict", "-", "paper 7", "excluded"],
        ],
    )
    negatives = write_csv(
        tmp_path / "neg.csv",
        ["old_locus_tag", "af3_classification", "af3_ipTM", "source"],
        [["MA_N1", "neg", "0.2", "alphafold3"], ["MA_N2", "neg", "0.3", "alphafold3"], ["MA_N3", "neg", "0.1", "alphafold3"], ["MA_N4", "neg", "0.1", "alphafold3"]],
    )
    return {"results": results, "curated": curated, "negatives": negatives, "out": tmp_path / "out"}


def run_small(case: dict[str, Path], **kwargs) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    cr.run(case["curated"], case["negatives"], case["results"], case["out"], negatives_query="MA_Q1", generated_at=GENERATED_AT, **kwargs)
    pairs = pd.read_excel(case["out"] / "pairs.xlsx", sheet_name="Pairs")
    negatives = pd.read_excel(case["out"] / "negatives_matched.xlsx", sheet_name="Negatives_Matched")
    return pairs, negatives, (case["out"] / "calibration_summary.md").read_text(encoding="utf-8")


def test_pairs_are_matched_in_the_directions_whose_query_is_in_the_results(small_case: dict[str, Path]) -> None:
    pairs, _, _ = run_small(small_case)

    by_key = {(r.query, r.candidate): r for r in pairs.itertuples()}
    assert by_key[("MA_Q1", "MA_P1")].status == "matched"
    assert float(by_key[("MA_Q1", "MA_P1")].interaction_score) == 60
    assert float(by_key[("MA_Q1", "MA_P1")].coexpression_gse77738) == pytest.approx(0.8)
    # Q1-Q2 is curated once but exists in the results only as Q2 -> Q1 (Q1 -> Q2 has no row)
    assert by_key[("MA_Q2", "MA_Q1")].status == "matched"
    assert by_key[("MA_Q1", "MA_Q2")].status == "candidate_not_in_results"
    assert by_key[("MA_Q1", "MA_X9")].status == "candidate_not_in_results"
    # neither protein queried / unresolved partner / excluded tier are reported, never scored
    assert by_key[("MA_Z1", "MA_Z2")].status == "no_query_in_results"
    assert by_key[("MA_Q1", "MA_P1/MA_P2")].status == "unresolved_multi_locus"
    assert by_key[("MA_Q1", "MA_P9")].status == "excluded"
    assert pairs[pairs.status != "matched"]["interaction_score"].isna().all()


def test_missing_components_stay_blank_not_zero(small_case: dict[str, Path]) -> None:
    pairs, _, _ = run_small(small_case)

    row = pairs[(pairs["query"] == "MA_Q2") & (pairs["candidate"] == "MA_Q1")].iloc[0]
    assert float(row.string_neighborhood) == pytest.approx(0.7)
    assert pd.isna(row.coexpression_gse77738)
    assert row.coexpression_gse77738_status == "MISSING"


def test_negatives_are_matched_and_multi_bucket_candidates_take_the_best_row(small_case: dict[str, Path]) -> None:
    _, negatives, _ = run_small(small_case)

    by_tag = {r.old_locus_tag: r for r in negatives.itertuples()}
    assert [str(r.found) for r in negatives.itertuples()] == ["True", "True", "True", "False"]
    assert float(by_tag["MA_N3"].final_score) == 35  # Negative_hit row (35) beats Candidates_relaxed (10)
    assert by_tag["MA_N3"].candidate_source == "Negative_hit"
    assert by_tag["MA_N3"].candidate_sources == "Candidates_relaxed;Negative_hit"
    assert pd.isna(by_tag["MA_N4"].interaction_score)
    assert set(negatives.query_old_locus_tag) == {"MA_Q1"}


def test_summary_statistics_match_hand_calculation(small_case: dict[str, Path]) -> None:
    _, _, summary = run_small(small_case)

    # Tier A = P1 (60), P2 (40), Q2->Q1 (50): mean 50, median 50. Negatives = N1 (10), N2 (30), N3 (20): mean 20, median 20.
    assert "| `interaction_score` | 3, 50.00 / 50.00 | 1, 20.00 / 20.00 | 3, 20.00 / 20.00 |" in summary
    # Tier B has only P3 (20).
    assert "| `interaction_priority_score` | 3, 30.00 / 30.00 | 1, 10.00 / 10.00 | 3, 15.00 / 15.00 |" in summary
    # Every Tier A value beats every negative -> U = 9 of 3*3, AUC 1.000, exact two-sided p = 2 / C(6,3) = 0.1.
    # Tier B (20) vs negatives (10, 20, 30): 1 win + 1 tie -> U = 1.5, AUC 0.500; the tie forces the normal approximation (p = 1).
    assert "| `interaction_score` | 1.000 | 0.100 (exact) | 0.500 | 1.000 (asymptotic) |" in summary
    # Component present on only some rows: coexpression_gse77738 has Tier A (0.8, 0.6) and negatives (0.9, 0.7); N3 has none.
    assert "| `coexpression_gse77738` | 2, 0.70 / 0.70 | 0, n/a / n/a | 2, 0.80 / 0.80 |" in summary
    # A = [0.8, 0.6] vs N = [0.9, 0.7]: only 0.8 > 0.7 -> U = 1 of 4, AUC 0.250 (exact p = 2/3)
    assert "| `coexpression_gse77738` | 0.250 | 0.667 (exact) |" in summary
    assert "Negatives found in the results: 3 of 4." in summary


def test_coverage_table_lists_every_status(small_case: dict[str, Path]) -> None:
    _, _, summary = run_small(small_case)

    for line in (
        "| A | matched | 3 |",
        "| A | candidate_not_in_results | 1 |",
        "| B | matched | 1 |",
        "| B | candidate_not_in_results | 1 |",
        "| B | no_query_in_results | 1 |",
        "| B | unresolved_multi_locus | 1 |",
        "| excluded | excluded | 1 |",
    ):
        assert line in summary


# ---------------------------------------------------------------------------
# Provenance sidecar
# ---------------------------------------------------------------------------


def write_sidecar(results: Path, *, check_passed: bool | None = True) -> None:
    payload = {
        "run_provenance": {
            "app_version": "5.0",
            "git_commit": "1a2b3c4",
            "git_dirty": True,
            "config_hash": "deadbeef",
            "reference_genome_check_passed": check_passed,
            "generated_at": "2026-09-20T12:00:00",
        },
        "effective_config": {
            "consider_cross_species_matches": True,
            "interaction_scoring": {
                "scoring_model": "v2_evidence_based",
                "max_candidates_per_query": 200,
                "candidate_sources": {"candidates": True, "no_hit": True, "negative_hit": False},
            },
        },
    }
    results.with_name(f"{results.stem}.run_provenance.yaml").write_text(yaml.safe_dump(payload), encoding="utf-8")


def test_summary_copies_provenance_when_the_sidecar_exists(small_case: dict[str, Path]) -> None:
    write_sidecar(small_case["results"])

    _, _, summary = run_small(small_case)

    assert "config hash: `deadbeef`" in summary
    assert "code: 5.0 (git 1a2b3c4+dirty)" in summary
    assert "reference genome check: passed" in summary
    assert "scoring_model=v2_evidence_based" in summary
    assert "candidate_sources=candidates,no_hit" in summary


def test_summary_flags_a_failed_reference_genome_check(small_case: dict[str, Path]) -> None:
    write_sidecar(small_case["results"], check_passed=False)

    _, _, summary = run_small(small_case)

    assert "reference genome check: FAILED" in summary


def test_summary_runs_without_a_sidecar_and_says_so(small_case: dict[str, Path]) -> None:
    _, _, summary = run_small(small_case)

    assert "run provenance: none" in summary
    assert "config hash" not in summary


def test_an_unreadable_sidecar_is_treated_as_absent(small_case: dict[str, Path]) -> None:
    small_case["results"].with_name("results.run_provenance.yaml").write_text("{ not: [valid", encoding="utf-8")

    _, _, summary = run_small(small_case)

    assert "run provenance: none" in summary


# ---------------------------------------------------------------------------
# Inputs, layout, CLI
# ---------------------------------------------------------------------------


def test_workbook_without_the_back_link_row_is_read_too(tmp_path: Path) -> None:
    results = tmp_path / "hand_made.xlsx"
    write_results_workbook(results, [pair_row("MA_Q1", "MA_P1", "Candidates", 60, 40, 70)], back_link_row=False)

    tables = cr.read_results(results)

    assert list(tables.pairs["candidate_tag"]) == ["MA_P1"]
    assert list(tables.pairs["query_tag"]) == ["MA_Q1"]


def test_missing_tier_column_is_a_clear_error(tmp_path: Path) -> None:
    curated = write_csv(tmp_path / "c.csv", CURATED_HEADER[:-1], [["MA_A", "a", "MA_B", "b", "strict", "", ""]])

    with pytest.raises(cr.CalibrationInputError, match="tier"):
        cr.read_curated(curated)


def test_invalid_tier_value_is_rejected(tmp_path: Path) -> None:
    curated = write_csv(tmp_path / "c.csv", CURATED_HEADER, [["MA_A", "a", "MA_B", "b", "strict", "", "", "C"]])

    with pytest.raises(cr.CalibrationInputError, match="tier must be"):
        cr.read_curated(curated)


def test_results_workbook_without_the_score_sheet_is_a_clear_error(tmp_path: Path) -> None:
    path = tmp_path / "other.xlsx"
    pd.DataFrame({"a": [1]}).to_excel(path, sheet_name="Sheet1", index=False)

    with pytest.raises(cr.CalibrationInputError, match="04_Score_Breakdown"):
        cr.read_results(path)


def test_existing_outputs_are_not_overwritten_without_force(small_case: dict[str, Path]) -> None:
    run_small(small_case)

    with pytest.raises(cr.CalibrationInputError, match="already contains"):
        cr.run(small_case["curated"], small_case["negatives"], small_case["results"], small_case["out"], negatives_query="MA_Q1")

    run_small(small_case, force=True)  # explicit overwrite is fine


def test_only_the_new_outputs_count_for_the_overwrite_guard(small_case: dict[str, Path]) -> None:
    """Leftover pairs.csv / negatives_matched.csv from before the Excel switch neither block a run nor get touched."""
    out = small_case["out"]
    out.mkdir()
    legacy = {name: out / name for name in ("pairs.csv", "negatives_matched.csv")}
    for path in legacy.values():
        path.write_text("old,csv\n1,2\n", encoding="utf-8")

    run_small(small_case)

    assert (out / "pairs.xlsx").exists() and (out / "negatives_matched.xlsx").exists()
    assert all(path.read_text(encoding="utf-8") == "old,csv\n1,2\n" for path in legacy.values())


def test_excel_tables_are_typed_formatted_and_single_sheet(small_case: dict[str, Path]) -> None:
    from openpyxl import load_workbook

    run_small(small_case)

    for name, sheet in (("pairs.xlsx", "Pairs"), ("negatives_matched.xlsx", "Negatives_Matched")):
        workbook = load_workbook(small_case["out"] / name)
        assert workbook.sheetnames == [sheet]
        worksheet = workbook[sheet]
        headers = {cell.value: cell for cell in worksheet[1]}
        assert all(cell.font.bold for cell in headers.values())
        assert worksheet.freeze_panes == "A2"
        assert worksheet.auto_filter.ref.startswith("A1:")

    worksheet = load_workbook(small_case["out"] / "pairs.xlsx")["Pairs"]
    columns = {cell.value: index for index, cell in enumerate(worksheet[1], start=1)}
    matched_row = next(r for r in range(2, worksheet.max_row + 1) if worksheet.cell(r, columns["status"]).value == "matched")
    score = worksheet.cell(matched_row, columns["interaction_score"])
    assert score.data_type == "n" and isinstance(score.value, (int, float))  # a number, not text
    unmatched_row = next(r for r in range(2, worksheet.max_row + 1) if worksheet.cell(r, columns["status"]).value == "excluded")
    assert worksheet.cell(unmatched_row, columns["interaction_score"]).value is None  # blank, not ""
    assert worksheet.cell(unmatched_row, columns["candidate_source"]).value is None

    negatives = load_workbook(small_case["out"] / "negatives_matched.xlsx")["Negatives_Matched"]
    neg_columns = {cell.value: index for index, cell in enumerate(negatives[1], start=1)}
    assert negatives.cell(2, neg_columns["found"]).data_type == "b"  # booleans stay booleans


def test_empty_tables_still_produce_valid_workbooks(tmp_path: Path) -> None:
    cr.write_excel_table(pd.DataFrame(columns=["a", "b"]), tmp_path / "empty.xlsx", "Pairs")

    frame = pd.read_excel(tmp_path / "empty.xlsx", sheet_name="Pairs")
    assert list(frame.columns) == ["a", "b"] and frame.empty


def test_cli_main_writes_the_three_outputs_and_returns_zero(small_case: dict[str, Path], capsys: pytest.CaptureFixture[str]) -> None:
    code = cr.main(
        [
            "--curated", str(small_case["curated"]),
            "--negatives", str(small_case["negatives"]),
            "--results", str(small_case["results"]),
            "--out", str(small_case["out"]),
            "--negatives-query", "MA_Q1",
        ]
    )

    assert code == 0
    assert {p.name for p in small_case["out"].iterdir()} == {"pairs.xlsx", "negatives_matched.xlsx", "calibration_summary.md"}
    assert "calibration_summary.md" in capsys.readouterr().out


def test_cli_main_reports_input_errors_without_a_traceback(small_case: dict[str, Path], capsys: pytest.CaptureFixture[str]) -> None:
    code = cr.main(
        ["--curated", str(small_case["curated"]), "--negatives", str(small_case["negatives"]),
         "--results", str(small_case["results"].with_name("missing.xlsx")), "--out", str(small_case["out"])]
    )

    assert code == 2
    assert "error:" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# Regression against the committed calibration report
# ---------------------------------------------------------------------------
# The original workbook (ProteinHunter_results_calibration_check.xlsx) was never
# committed, so it is rebuilt here from the report's own per-pair CSVs: the tool
# must reproduce the figures printed in
# claude/experimental_interactions_calibration_report.md from those inputs.


def _float_or_none(value: str) -> float | None:
    return float(value) if value not in ("", None) else None


def _components_from_row(row: dict[str, str]) -> dict[str, float | None]:
    return {name: _float_or_none(row.get(name, "")) for name in cr.COMPONENT_COLUMNS if name in row}


def build_workbook_from_report_csvs(path: Path) -> None:
    rows = []
    with (REPO_ROOT / "claude/experimental_interactions_calibration_report_pairs.csv").open(encoding="utf-8-sig", newline="") as handle:
        for r in csv.DictReader(handle):
            rows.append(
                pair_row(r["query"], r["candidate"], "Candidates", float(r["interaction_score"]), float(r["interaction_priority_score"]), float("nan"), **_components_from_row(r))
            )
    with (REPO_ROOT / "claude/experimental_interactions_calibration_report_negatives.csv").open(encoding="utf-8-sig", newline="") as handle:
        for r in csv.DictReader(handle):
            rows.append(
                pair_row("MA_4115", r["old_locus_tag"], "Candidates", float(r["interaction_score"]), float(r["interaction_priority_score"]), float("nan"), **_components_from_row(r))
            )
    write_results_workbook(path, rows)


def test_reproduces_the_figures_in_the_committed_calibration_report(tmp_path: Path) -> None:
    results = tmp_path / "calibration_check.xlsx"
    build_workbook_from_report_csvs(results)

    # The curated pairs as they were when the report was written (new confirmed pairs have been added to
    # claude/experimental_interactions_curated.csv since); the report's figures only follow from that snapshot.
    cr.run(
        REPO_ROOT / "tests/fixtures/curated_pairs_calibration_report.csv",
        REPO_ROOT / "claude/calibration/af3_negatives_MA_4115.csv",
        results,
        tmp_path / "out",
        generated_at=GENERATED_AT,
    )
    summary = (tmp_path / "out/calibration_summary.md").read_text(encoding="utf-8")

    # "Tier A summary (n=8): interaction_score mean 39.55, median 43.85"; AF3 negatives (n=28): mean 12.95, median 12.27
    assert "| `interaction_score` | 8, 39.55 / 43.85 | 17, " in summary
    assert "| 28, 12.95 / 12.27 |" in summary
    # "interaction_priority_score ... 17.97 vs. 20.55"
    assert "| `interaction_priority_score` | 8, 17.97 /" in summary
    assert "| 28, 20.55 /" in summary
    # "string_neighborhood ... (0.62 vs. 0.04)", "coexpression_gse77738 ... (0.48 positives vs. 0.65 negatives)"
    assert "| `string_neighborhood` | 8, 0.62 /" in summary and "| 28, 0.04 /" in summary
    assert "| `coexpression_gse77738` | 8, 0.48 /" in summary and "| 28, 0.65 /" in summary
    # "coexpression_gse64349 ... (0.85 vs. 0.60)", only the Tier A pairs with a value counted
    assert "| `coexpression_gse64349` | 3, 0.85 /" in summary and "| 28, 0.60 /" in summary

    pairs = pd.read_excel(tmp_path / "out/pairs.xlsx", sheet_name="Pairs")
    matched = pairs[pairs.status == "matched"]
    assert (matched.tier == "A").sum() == 8
    assert (matched.tier == "B").sum() == 17  # 12 original Tier B rows + 5 demoted from Tier A (CdhC x2, DnaK, Hsp20, MtsF)
    # interaction_score separates Tier A from the negatives: report says ~3x on the means
    tier_a = pd.to_numeric(matched[matched.tier == "A"].interaction_score)
    negatives = pd.read_excel(tmp_path / "out/negatives_matched.xlsx", sheet_name="Negatives_Matched")
    result = cr.mann_whitney(tier_a, pd.to_numeric(negatives.interaction_score))
    assert result.n1 == 8 and result.n2 == 28
    assert result.auc > 0.8


# ---------------------------------------------------------------------------
# The committed curated pairs
# ---------------------------------------------------------------------------


def test_committed_curated_pairs_are_valid_and_have_no_duplicate_pairs() -> None:
    frame = cr.read_curated(REPO_ROOT / "claude/experimental_interactions_curated.csv")

    pairs = [frozenset((a, b)) for a, b in zip(frame["protein_a_old_locus_tag"], frame["protein_b_old_locus_tag"])]
    assert len(pairs) == len(set(pairs)), "a pair listed twice would be counted twice (once per tier)"
    assert set(frame["tier"]) <= {"A", "B", "excluded"}


def test_component_a2_mcr_pairs_are_tier_a() -> None:
    """Current Biology 2026 (PMC13432166): A2 (atwA, MA_3998) x McrA/McrB/McrG (MA_4546 / MA_4550 / MA_4547), M. acetivorans WWM73."""
    frame = cr.read_curated(REPO_ROOT / "claude/experimental_interactions_curated.csv")

    a2 = frame[frame["protein_a_old_locus_tag"] == "MA_3998"].to_dict("records")  # "class" is not a valid attribute name

    assert {(r["protein_b_old_locus_tag"], r["protein_b_label"], r["tier"], r["class"]) for r in a2} == {
        ("MA_4546", "McrA", "A", "strict"),
        ("MA_4550", "McrB", "A", "strict"),
        ("MA_4547", "McrG", "A", "strict"),
    }
    assert all("PMC13432166" in r["source"] and "WWM73" in r["confidence"] for r in a2)

