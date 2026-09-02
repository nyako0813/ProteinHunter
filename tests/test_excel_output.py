"""Tests for minimal Excel output helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from openpyxl import load_workbook

from core.exceptions import ExcelOutputError
from core.models import BlastHit, CandidateScore, DomainHit, ProteinRecord
from analysis.interaction_scoring import (
    INTERACTION_EVIDENCE_DETAIL_LEGACY_COLUMNS,
    INTERACTION_EVIDENCE_DETAIL_V2_COLUMNS,
    InteractionScoringResult,
    interaction_pair_columns,
)
from output.excel import (
    EXCEL_COLUMNS,
    INDEX_ROWS,
    records_to_dataframe,
    write_classification_workbook,
    write_records_to_excel,
)


def make_hit(
    subject_id: str,
    bitscore: float,
    evalue: float,
    source: str = "positive",
) -> BlastHit:
    """Create a BLAST hit for Excel output tests."""
    return BlastHit(
        query_id="protein_1",
        subject_id=subject_id,
        percent_identity=80.0,
        alignment_length=100,
        evalue=evalue,
        bitscore=bitscore,
        source=source,
    )


def make_record() -> ProteinRecord:
    """Create a protein record with enough fields for Excel tests."""
    record = ProteinRecord(
        protein_id="protein_1",
        description="candidate protein",
        old_locus_tag="MA_1234",
        sequence="MSTNPKPQR",
        positive_hits=[
            make_hit("positive_low", 20.0, 1e-50),
            make_hit("positive_best", 50.0, 1e-5),
            make_hit("positive_tie_best", 50.0, 1e-20),
        ],
        negative_hits=[make_hit("negative_best", 10.0, 1e-3, source="negative")],
        domains=[
            DomainHit(
                source="CDD",
                accession="cd12345",
                name="Thioredoxin_like",
                description="redox domain",
            ),
            DomainHit(
                source="Pfam",
                accession="PF00001",
                name="Domain",
                description="",
            ),
            DomainHit(
                source="CDD",
                accession="cd67890",
                name="Second_domain",
                description="second domain",
            ),
        ],
        motifs=["CXXC", "HXXH"],
        uniprot_accession="P12345",
        alphafold_url="https://example.test/model",
        notes=["reviewed", "export ready"],
    )
    score = CandidateScore(protein_id="protein_1")
    score.add_component("positive_hit", 5.0, "Positive BLAST hit found")
    score.add_component("domain_hit", 4.0, "Domain annotation found")
    record.score = score

    return record


def test_records_to_dataframe_column_order() -> None:
    """DataFrame columns should match the required Excel column order."""
    dataframe = records_to_dataframe({"protein_1": make_record()})

    assert list(dataframe.columns) == list(EXCEL_COLUMNS)
    assert list(dataframe.columns)[2] == "old_locus_tag"
    domain_count_index = list(dataframe.columns).index("domain_count")
    assert list(dataframe.columns)[domain_count_index + 1 : domain_count_index + 4] == [
        "unique_domain_count",
        "unique_domain_accessions",
        "unique_domain_names",
    ]


def test_empty_records_returns_expected_columns() -> None:
    """Empty records should return an empty DataFrame with stable columns."""
    dataframe = records_to_dataframe({})

    assert dataframe.empty
    assert list(dataframe.columns) == list(EXCEL_COLUMNS)


def test_records_to_dataframe_includes_score_fields() -> None:
    """Score details should be included in stable Excel columns."""
    dataframe = records_to_dataframe({"protein_1": make_record()})
    row = dataframe.iloc[0]

    assert row["total_score"] == 9.0
    assert row["score_components"] == "positive_hit=5.0; domain_hit=4.0"
    assert row["score_reasons"] == "Positive BLAST hit found; Domain annotation found"


def test_records_to_dataframe_without_score_uses_empty_score_fields() -> None:
    """Records without scores should use zero and blank score fields."""
    record = ProteinRecord(protein_id="protein_1")

    dataframe = records_to_dataframe({"protein_1": record})
    row = dataframe.iloc[0]

    assert row["total_score"] == 0
    assert row["score_components"] == ""
    assert row["score_reasons"] == ""
    assert row["old_locus_tag"] == ""


def test_records_to_dataframe_includes_old_locus_tag() -> None:
    """old_locus_tag should be exported immediately after description."""
    dataframe = records_to_dataframe({"protein_1": make_record()})
    row = dataframe.iloc[0]

    assert row["old_locus_tag"] == "MA_1234"
    assert list(dataframe.columns)[1:3] == ["description", "old_locus_tag"]


def test_records_to_dataframe_appends_positive_source_fields() -> None:
    """Positive source summary fields should be appended after existing columns."""
    record = ProteinRecord(
        protein_id="protein_1",
        positive_source_count=2,
        positive_sources_hit=["A", "B"],
        positive_sources_missing=["C"],
    )

    dataframe = records_to_dataframe({"protein_1": record})
    row = dataframe.iloc[0]

    assert list(dataframe.columns)[-3:] == [
        "positive_source_count",
        "positive_sources_hit",
        "positive_sources_missing",
    ]
    assert row["positive_source_count"] == 2
    assert row["positive_sources_hit"] == "A; B"
    assert row["positive_sources_missing"] == "C"


def test_records_to_dataframe_includes_negative_evidence_fields() -> None:
    """Negative hit evidence fields should be exported for classification sheets."""
    record = make_record()
    record.negative_best_identity = 45.0
    record.negative_best_query_coverage = 80.0
    record.negative_best_evalue = 1e-20
    record.negative_best_source = "Negative_source"
    record.negative_hit_strength = "strong"
    record.negative_strong_hit_count = 1
    record.negative_exclusion_reason = "excluded: strong negative hit"

    dataframe = records_to_dataframe({"protein_1": record})
    row = dataframe.iloc[0]

    assert row["negative_best_identity"] == 45.0
    assert row["negative_best_query_coverage"] == 80.0
    assert row["negative_best_evalue"] == 1e-20
    assert row["negative_best_source"] == "Negative_source"
    assert row["negative_hit_strength"] == "strong"
    assert row["negative_strong_hit_count"] == 1
    assert row["negative_exclusion_reason"] == "excluded: strong negative hit"


def test_records_to_dataframe_includes_domain_fields() -> None:
    """Domain details should be joined into stable Excel columns."""
    dataframe = records_to_dataframe({"protein_1": make_record()})
    row = dataframe.iloc[0]

    assert row["domain_sources"] == "CDD; Pfam"
    assert row["domain_names"] == "Thioredoxin_like; Domain; Second_domain"
    assert row["domain_accessions"] == "cd12345; PF00001; cd67890"
    assert row["domain_descriptions"] == "redox domain; second domain"
    assert row["domain_count"] == 3
    assert row["unique_domain_count"] == 3
    assert row["unique_domain_accessions"] == "cd12345; PF00001; cd67890"
    assert row["unique_domain_names"] == "Thioredoxin_like; Domain; Second_domain"


def test_records_to_dataframe_includes_unique_domain_summary() -> None:
    """Duplicate domain accessions and names should be summarized once."""
    record = ProteinRecord(
        protein_id="protein_1",
        domains=[
            DomainHit(source="Pfam", accession="PF00001", name="ABC"),
            DomainHit(source="Pfam", accession="PF00001", name="ABC"),
            DomainHit(source="CDD", accession="cd12345", name="CDD domain"),
        ],
    )

    dataframe = records_to_dataframe({"protein_1": record})
    row = dataframe.iloc[0]

    assert row["domain_count"] == 3
    assert row["unique_domain_count"] == 2
    assert row["unique_domain_accessions"] == "PF00001; cd12345"
    assert row["unique_domain_names"] == "ABC; CDD domain"


def test_records_to_dataframe_unique_names_skip_numeric_internal_ids() -> None:
    """Unique domain names should not include numeric-only internal ids."""
    record = ProteinRecord(
        protein_id="protein_1",
        domains=[
            DomainHit(source="Pfam", accession="PF01637.24", name="000001295"),
            DomainHit(source="Pfam", accession="PF03008.20", name="ABC_transporter"),
        ],
    )

    dataframe = records_to_dataframe({"protein_1": record})
    row = dataframe.iloc[0]

    assert row["unique_domain_names"] == "ABC_transporter"


def test_records_to_dataframe_empty_domains_are_blank() -> None:
    """Records without domains should use blank domain fields and count zero."""
    record = ProteinRecord(protein_id="protein_1", description="no domains")

    dataframe = records_to_dataframe({"protein_1": record})
    row = dataframe.iloc[0]

    assert row["domain_sources"] == ""
    assert row["domain_names"] == ""
    assert row["domain_accessions"] == ""
    assert row["domain_descriptions"] == ""
    assert row["domain_count"] == 0
    assert row["unique_domain_count"] == 0
    assert row["unique_domain_accessions"] == ""
    assert row["unique_domain_names"] == ""


def test_records_to_dataframe_selects_best_positive_hit() -> None:
    """Best positive hit should use highest bitscore, then lower e-value."""
    dataframe = records_to_dataframe({"protein_1": make_record()})
    row = dataframe.iloc[0]

    assert row["best_positive_hit"] == "positive_tie_best"
    assert row["best_positive_bitscore"] == 50.0
    assert row["best_positive_evalue"] == 1e-20


def test_records_to_dataframe_joins_motifs_and_notes() -> None:
    """Motifs and notes should be joined with a semicolon separator."""
    dataframe = records_to_dataframe({"protein_1": make_record()})
    row = dataframe.iloc[0]

    assert row["motifs"] == "CXXC; HXXH"
    assert row["notes"] == "reviewed; export ready"




def test_interaction_pair_columns_include_distance_independent_ranking() -> None:
    """Interaction sheets should expose distance-independent ranking fields."""
    columns = interaction_pair_columns(False)

    assert "distance_independent_score" in columns
    assert "distance_independent_rank" in columns
    assert "priority_group" in columns

def test_write_records_to_excel_creates_xlsx_file(tmp_path: Path) -> None:
    """Excel writer should create a readable xlsx file."""
    output_path = tmp_path / "reports" / "candidates.xlsx"

    result = write_records_to_excel({"protein_1": make_record()}, output_path)

    assert result == output_path.resolve()
    assert result.exists()

    dataframe = pd.read_excel(result, sheet_name="Candidates")
    assert list(dataframe.columns) == list(EXCEL_COLUMNS)
    assert dataframe.loc[0, "protein_id"] == "protein_1"


def test_write_classification_workbook_creates_expected_sheets(
    tmp_path: Path,
) -> None:
    """Classification workbook should include all BLAST category sheets."""
    candidate = ProteinRecord(
        protein_id="A_positive_only",
        positive_hits=[make_hit("positive", 50.0, 1e-20)],
        positive_source_count=1,
        positive_sources_hit=["A"],
    )
    no_hit = ProteinRecord(protein_id="B_no_hits")
    negative_only = ProteinRecord(
        protein_id="C_negative_only",
        negative_hits=[make_hit("negative", 40.0, 1e-10, source="negative")],
    )
    both = ProteinRecord(
        protein_id="D_both",
        positive_hits=[make_hit("positive", 50.0, 1e-20)],
        negative_hits=[make_hit("negative", 40.0, 1e-10, source="negative")],
    )
    output_path = tmp_path / "reports" / "classification.xlsx"

    result = write_classification_workbook(
        candidates={"A_positive_only": candidate},
        output_path=output_path,
        positive_all_sources={"A_positive_only": candidate},
        positive_source_summary={
            "A_positive_only": candidate,
            "B_no_hits": no_hit,
            "C_negative_only": negative_only,
            "D_both": both,
        },
        negative_unmatched={
            "A_positive_only": candidate,
            "B_no_hits": no_hit,
        },
        no_hit={"B_no_hits": no_hit},
        negative_hit={
            "C_negative_only": negative_only,
            "D_both": both,
        },
    )

    workbook = load_workbook(result)
    assert workbook.sheetnames == [
        "Index",
        "Candidates",
        "Candidates_relaxed",
        "Positive_all_sources",
        "Positive_source_summary",
        "Negative_unmatched",
        "No_hit",
        "Negative_hit",
        "Negative_strong_hit",
        "Negative_medium_hit",
        "Negative_weak_hit",
    ]
    assert workbook["Candidates"].max_row == 3
    assert workbook["Positive_all_sources"].max_row == 3
    assert workbook["Positive_source_summary"].max_row == 6
    assert workbook["Negative_unmatched"].max_row == 4
    assert workbook["No_hit"].max_row == 3
    assert workbook["Negative_hit"].max_row == 4

    negative_hit = pd.read_excel(result, sheet_name="Negative_hit", header=1)
    assert set(negative_hit["protein_id"]) == {"C_negative_only", "D_both"}




def test_index_sheet_explains_interaction_scoring_columns(tmp_path: Path) -> None:
    """Index should include short explanations for key interaction scoring columns."""
    output_path = tmp_path / "reports" / "index_scoring_explanations.xlsx"

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
    )

    workbook = load_workbook(result)
    index_values = [
        str(cell.value)
        for row in workbook["Index"].iter_rows()
        for cell in row
        if cell.value is not None
    ]
    index_text = "\n".join(index_values)

    assert "interaction_priority_score" in index_text
    assert "distance_independent_score" in index_text
    assert "priority_group" in index_text
    assert "protein_hunter_score" in index_text
    assert "alphafold_readiness_score" in index_text
    assert "not a direct protein-protein interaction probability" in index_text
    assert "string_ppi_score" in index_text
    assert "string-db.org" in index_text
    assert "CC BY 4.0" in index_text

def test_classification_workbook_index_links_all_sheets(tmp_path: Path) -> None:
    """Index should be first and link to every classification sheet."""
    output_path = tmp_path / "reports" / "classification_links.xlsx"

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
    )

    workbook = load_workbook(result)
    assert workbook.sheetnames[0] == "Index"
    index = workbook["Index"]
    expected_sheets = [row[0] for row in INDEX_ROWS]
    linked_sheets = [index.cell(row=row_index, column=1).value for row_index in range(2, 12)]
    assert linked_sheets == expected_sheets
    for row_index, sheet_name in enumerate(expected_sheets, start=2):
        assert index.cell(row=row_index, column=1).hyperlink.target == (
            f"#'{sheet_name}'!A1"
        )

    for sheet_name in expected_sheets:
        worksheet = workbook[sheet_name]
        assert worksheet["A1"].value == "Back to Index"
        assert worksheet["A1"].hyperlink.target == "#'Index'!A1"


def test_classification_workbook_adds_only_created_interaction_sheets(
    tmp_path: Path,
) -> None:
    """Interaction sheets should appear in Index only when actually created."""
    output_path = tmp_path / "reports" / "interaction.xlsx"
    interaction_result = InteractionScoringResult(
        query_rows=[
            {
                "query_id": "query_1",
                "input_protein_id": "query_1",
                "input_old_locus_tag": "",
                "resolved_protein_id": "query_1",
                "resolved_old_locus_tag": "",
                "sequence_length": 10,
                "resolution_status": "resolved",
                "description": "query",
                "notes": "",
            }
        ],
        source_rows={
            "Interaction_Positive_all": [
                {
                    "query_id": "query_1",
                    "query_protein_id": "query_1",
                    "query_old_locus_tag": "",
                    "candidate_rank": 1,
                    "candidate_protein_id": "candidate_1",
                    "candidate_old_locus_tag": "",
                    "candidate_source": "Positive_all_sources",
                    "candidate_description": "candidate",
                    "interaction_priority_score": 42.0,
                    "interaction_score_reasons": "candidate source: Positive_all_sources",
                    "candidate_priority_score": 30.0,
                    "same_gene_neighborhood_score": 0.0,
                    "distance_bp": None,
                    "co_occurrence_score": 0.0,
                    "domain_complementarity_score": 0.0,
                    "alphafold_readiness_score": 10.0,
                    "pair_total_length": 20,
                    "alphafold_recommended": True,
                }
            ]
        },
        neighborhood_rows=[],
        warnings=[],
    )

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
        interaction_result=interaction_result,
    )

    workbook = load_workbook(result)
    assert "Interaction_query" in workbook.sheetnames
    assert "Interaction_Positive_all" in workbook.sheetnames
    assert "Interaction_Positive_all_sources" not in workbook.sheetnames
    assert "Interaction_No_hit" not in workbook.sheetnames
    index_values = [cell.value for cell in workbook["Index"]["A"]]
    assert "Interaction_query" in index_values
    assert "Interaction_Positive_all" in index_values
    assert "Interaction_Positive_all_sources" not in index_values
    assert "Interaction_No_hit" not in index_values
    assert workbook["Interaction_query"]["A1"].hyperlink.target == "#'Index'!A1"
    assert workbook["Interaction_Positive_all"]["A1"].hyperlink.target == "#'Index'!A1"


def test_write_records_to_excel_applies_simple_formatting(tmp_path: Path) -> None:
    """Excel output should include basic readability formatting."""
    output_path = tmp_path / "reports" / "formatted_candidates.xlsx"

    result = write_records_to_excel({"protein_1": make_record()}, output_path)

    workbook = load_workbook(result)
    worksheet = workbook["Candidates"]

    assert worksheet.freeze_panes == "A2"
    assert worksheet.auto_filter.ref is not None
    assert worksheet["A1"].font.bold is True
    assert all(cell.font.bold for cell in worksheet[1])
    assert [cell.value for cell in worksheet[1]] == list(EXCEL_COLUMNS)
    assert worksheet.column_dimensions["B"].width >= 35


def test_write_records_to_excel_raises_excel_output_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Writer failures should be wrapped in ExcelOutputError."""

    class BrokenExcelWriter:
        """Context manager that fails while opening the workbook."""

        def __init__(self, *args: object, **kwargs: object) -> None:
            pass

        def __enter__(self) -> "BrokenExcelWriter":
            raise OSError("cannot write workbook")

        def __exit__(self, *args: object) -> None:
            pass

    monkeypatch.setattr("output.excel.pd.ExcelWriter", BrokenExcelWriter)

    with pytest.raises(ExcelOutputError, match="could not write the Excel file"):
        write_records_to_excel({"protein_1": make_record()}, tmp_path / "bad.xlsx")



def test_index_hyperlinks_point_to_existing_sheets(tmp_path: Path) -> None:
    """Every Index hyperlink should point to a sheet that exists."""
    output_path = tmp_path / "reports" / "interaction_links.xlsx"
    interaction_result = InteractionScoringResult(
        query_rows=[
            {
                "query_id": "query_1",
                "input_protein_id": "query_1",
                "input_old_locus_tag": "",
                "resolved_protein_id": "query_1",
                "resolved_old_locus_tag": "",
                "sequence_length": 10,
                "resolution_status": "resolved",
                "description": "query",
                "notes": "",
            }
        ],
        source_rows={
            "Interaction_Positive_all": [
                {
                    "query_id": "query_1",
                    "query_protein_id": "query_1",
                    "query_old_locus_tag": "",
                    "candidate_rank": 1,
                    "candidate_protein_id": "candidate_1",
                    "candidate_old_locus_tag": "",
                    "candidate_source": "Positive_all_sources",
                    "candidate_description": "candidate",
                    "same_contig": True,
                    "query_start": 1,
                    "query_end": 100,
                    "query_strand": "+",
                    "candidate_start": 200,
                    "candidate_end": 300,
                    "candidate_strand": "+",
                    "distance_bp": 100,
                    "strand_relation": "same_strand",
                    "same_gene_neighborhood_score": 25.0,
                    "interaction_priority_score": 42.0,
                    "interaction_score_reasons": "candidate source: Positive_all_sources",
                    "candidate_priority_score": 30.0,
                    "co_occurrence_score": 0.0,
                    "domain_complementarity_score": 0.0,
                    "alphafold_readiness_score": 10.0,
                    "pair_total_length": 20,
                    "alphafold_recommended": True,
                }
            ]
        },
        neighborhood_rows=[
            {
                "query_id": "query_1",
                "query_protein_id": "query_1",
                "query_old_locus_tag": "",
                "query_description": "query",
                "query_contig": "contig1",
                "query_start": 1,
                "query_end": 100,
                "query_strand": "+",
                "candidate_rank_by_distance": 1,
                "candidate_protein_id": "candidate_1",
                "candidate_old_locus_tag": "",
                "candidate_description": "candidate",
                "candidate_source": "Positive_all_sources",
                "candidate_contig": "contig1",
                "candidate_start": 200,
                "candidate_end": 300,
                "candidate_strand": "+",
                "distance_bp": 100,
                "strand_relation": "same_strand",
                "neighborhood_band": "<=5kb",
                "same_gene_neighborhood_score": 25.0,
                "interaction_priority_score": 42.0,
                "domain_complementarity_score": 0.0,
                "candidate_priority_score": 30.0,
                "co_occurrence_score": 0.0,
                "alphafold_recommended": True,
                "interaction_score_reasons": "candidate source: Positive_all_sources",
            }
        ],
        warnings=[],
    )

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
        interaction_result=interaction_result,
    )

    workbook = load_workbook(result)
    sheet_names = set(workbook.sheetnames)
    for cell in workbook["Index"]["A"]:
        if cell.hyperlink is None:
            continue
        target = cell.hyperlink.target
        linked_sheet = target.split("'")[1]
        assert linked_sheet in sheet_names
    assert "Interaction_Positive_all" in sheet_names
    assert "Interaction_Neighborhood" in sheet_names
    assert workbook["Interaction_Neighborhood"]["A1"].hyperlink.target == "#'Index'!A1"
    assert "Interaction_Positive_all_sources" not in sheet_names


def _minimal_pair_row(candidate_protein_id: str = "candidate_1") -> dict:
    """A minimal but complete Interaction_Candidates-shaped row for Excel tests."""
    return {
        "query_id": "query_1",
        "query_protein_id": "query_1",
        "query_old_locus_tag": "",
        "candidate_rank": 1,
        "candidate_protein_id": candidate_protein_id,
        "candidate_old_locus_tag": "",
        "candidate_source": "Candidates",
        "candidate_description": "candidate",
        "interaction_priority_score": 42.0,
        "interaction_score_reasons": "candidate source: Candidates",
        "candidate_priority_score": 30.0,
        "same_gene_neighborhood_score": 0.0,
        "distance_bp": None,
        "co_occurrence_score": 0.0,
        "domain_complementarity_score": 0.0,
        "alphafold_readiness_score": 10.0,
        "pair_total_length": 20,
        "alphafold_recommended": True,
        "protein_hunter_score": 14.0,
        "protein_hunter_score_components": "no_negative_hit=5.0; domain_hit=4.0",
        "protein_hunter_score_reasons": "This protein has no negative BLAST hits.",
    }


def test_interaction_candidates_sheet_shows_protein_hunter_score_column(
    tmp_path: Path,
) -> None:
    """protein_hunter_score reference columns should render in Interaction_* sheets."""
    output_path = tmp_path / "reports" / "protein_hunter_score.xlsx"
    interaction_result = InteractionScoringResult(
        query_rows=[],
        source_rows={"Interaction_Candidates": [_minimal_pair_row()]},
        neighborhood_rows=[],
        warnings=[],
    )

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
        interaction_result=interaction_result,
    )

    dataframe = pd.read_excel(result, sheet_name="Interaction_Candidates", header=1)
    assert dataframe.loc[0, "protein_hunter_score"] == 14.0
    assert dataframe.loc[0, "protein_hunter_score_components"] == (
        "no_negative_hit=5.0; domain_hit=4.0"
    )
    # interaction_priority_score is unaffected by the new reference columns.
    assert dataframe.loc[0, "interaction_priority_score"] == 42.0


def test_classification_workbook_writes_v2_evidence_detail_sheet(tmp_path: Path) -> None:
    """Interaction_Evidence_Detail should hold one row per component for v2 runs."""
    output_path = tmp_path / "reports" / "evidence_detail_v2.xlsx"
    interaction_result = InteractionScoringResult(
        query_rows=[],
        source_rows={"Interaction_Candidates": [_minimal_pair_row()]},
        neighborhood_rows=[],
        warnings=[],
        evidence_detail_rows=[
            {
                "query_id": "query_1",
                "query_protein_id": "query_1",
                "query_old_locus_tag": "",
                "candidate_protein_id": "candidate_1",
                "candidate_old_locus_tag": "",
                "candidate_source": "Candidates",
                "candidate_rank": 1,
                "category": "source_classification",
                "component_name": "source_classification",
                "status": "AVAILABLE",
                "raw_value": "Candidates",
                "normalized_value": 1.0,
                "weight": 1.0,
                "category_cap": 30.0,
                "is_negative": False,
                "explanation": "candidate source: Candidates",
            },
            {
                "query_id": "query_1",
                "query_protein_id": "query_1",
                "query_old_locus_tag": "",
                "candidate_protein_id": "candidate_1",
                "candidate_old_locus_tag": "",
                "candidate_source": "Candidates",
                "candidate_rank": 1,
                "category": "genomic_context",
                "component_name": "genomic_context",
                "status": "MISSING",
                "raw_value": None,
                "normalized_value": None,
                "weight": 0.0,
                "category_cap": None,
                "is_negative": False,
                "explanation": "no GFF coordinates",
            },
        ],
        evidence_detail_scoring_model="v2_evidence_based",
    )

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
        interaction_result=interaction_result,
    )

    workbook = load_workbook(result)
    assert "Interaction_Evidence_Detail" in workbook.sheetnames
    # Existing sheets are unaffected by adding the detail sheet.
    assert workbook["Interaction_Candidates"].max_row == 3  # back-link + header + 1 row

    detail = pd.read_excel(result, sheet_name="Interaction_Evidence_Detail", header=1)
    assert list(detail.columns) == list(INTERACTION_EVIDENCE_DETAIL_V2_COLUMNS)
    assert len(detail) == 2
    assert set(detail["component_name"]) == {"source_classification", "genomic_context"}
    source_row = detail[detail["component_name"] == "source_classification"].iloc[0]
    assert source_row["status"] == "AVAILABLE"
    assert source_row["category_cap"] == 30.0

    index_values = [cell.value for cell in workbook["Index"]["A"]]
    assert "Interaction_Evidence_Detail" in index_values
    assert workbook["Interaction_Evidence_Detail"]["A1"].hyperlink.target == "#'Index'!A1"


def test_classification_workbook_writes_legacy_evidence_detail_sheet(tmp_path: Path) -> None:
    """Interaction_Evidence_Detail should be a wide, one-row-per-pair projection for legacy runs."""
    output_path = tmp_path / "reports" / "evidence_detail_legacy.xlsx"
    interaction_result = InteractionScoringResult(
        query_rows=[],
        source_rows={"Interaction_Candidates": [_minimal_pair_row()]},
        neighborhood_rows=[],
        warnings=[],
        evidence_detail_rows=[
            {
                "query_id": "query_1",
                "query_protein_id": "query_1",
                "query_old_locus_tag": "",
                "candidate_protein_id": "candidate_1",
                "candidate_old_locus_tag": "",
                "candidate_source": "Candidates",
                "candidate_rank": 1,
                "candidate_priority_score": 30.0,
                "same_gene_neighborhood_score": 0.0,
                "co_occurrence_score": 0.0,
                "domain_complementarity_score": 0.0,
                "alphafold_readiness_score": 10.0,
                "interaction_score_reasons": "candidate source: Candidates",
            }
        ],
        evidence_detail_scoring_model="legacy_additive",
    )

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
        interaction_result=interaction_result,
    )

    workbook = load_workbook(result)
    assert "Interaction_Evidence_Detail" in workbook.sheetnames
    detail = pd.read_excel(result, sheet_name="Interaction_Evidence_Detail", header=1)
    assert list(detail.columns) == list(INTERACTION_EVIDENCE_DETAIL_LEGACY_COLUMNS)
    assert len(detail) == 1
    assert detail.iloc[0]["candidate_priority_score"] == 30.0


def test_classification_workbook_omits_evidence_detail_sheet_when_empty(tmp_path: Path) -> None:
    """No Interaction_Evidence_Detail sheet should be written when there are no detail rows."""
    output_path = tmp_path / "reports" / "evidence_detail_empty.xlsx"
    interaction_result = InteractionScoringResult(
        query_rows=[],
        source_rows={"Interaction_Candidates": [_minimal_pair_row()]},
        neighborhood_rows=[],
        warnings=[],
    )

    result = write_classification_workbook(
        candidates={},
        output_path=output_path,
        interaction_result=interaction_result,
    )

    workbook = load_workbook(result)
    assert "Interaction_Evidence_Detail" not in workbook.sheetnames
    assert "Interaction_Evidence_Detail" not in [cell.value for cell in workbook["Index"]["A"]]
    # Existing sheets are still written normally.
    assert "Interaction_Candidates" in workbook.sheetnames
