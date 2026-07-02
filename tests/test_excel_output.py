"""Tests for minimal Excel output helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from openpyxl import load_workbook

from core.exceptions import ExcelOutputError
from core.models import BlastHit, CandidateScore, DomainHit, ProteinRecord
from output.excel import EXCEL_COLUMNS, records_to_dataframe, write_records_to_excel


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


def test_write_records_to_excel_creates_xlsx_file(tmp_path: Path) -> None:
    """Excel writer should create a readable xlsx file."""
    output_path = tmp_path / "reports" / "candidates.xlsx"

    result = write_records_to_excel({"protein_1": make_record()}, output_path)

    assert result == output_path.resolve()
    assert result.exists()

    dataframe = pd.read_excel(result, sheet_name="Candidates")
    assert list(dataframe.columns) == list(EXCEL_COLUMNS)
    assert dataframe.loc[0, "protein_id"] == "protein_1"


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
