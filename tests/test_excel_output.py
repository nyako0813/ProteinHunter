"""Tests for minimal Excel output helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from core.exceptions import ExcelOutputError
from core.models import BlastHit, DomainHit, ProteinRecord
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
    return ProteinRecord(
        protein_id="protein_1",
        description="candidate protein",
        sequence="MSTNPKPQR",
        positive_hits=[
            make_hit("positive_low", 20.0, 1e-50),
            make_hit("positive_best", 50.0, 1e-5),
            make_hit("positive_tie_best", 50.0, 1e-20),
        ],
        negative_hits=[make_hit("negative_best", 10.0, 1e-3, source="negative")],
        domains=[DomainHit(source="pfam", accession="PF00001", name="Domain")],
        motifs=["CXXC", "HXXH"],
        uniprot_accession="P12345",
        alphafold_url="https://example.test/model",
        notes=["reviewed", "export ready"],
    )


def test_records_to_dataframe_column_order() -> None:
    """DataFrame columns should match the required Excel column order."""
    dataframe = records_to_dataframe({"protein_1": make_record()})

    assert list(dataframe.columns) == list(EXCEL_COLUMNS)


def test_empty_records_returns_expected_columns() -> None:
    """Empty records should return an empty DataFrame with stable columns."""
    dataframe = records_to_dataframe({})

    assert dataframe.empty
    assert list(dataframe.columns) == list(EXCEL_COLUMNS)


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
