"""Tests for directory-based FASTA input sources."""

from __future__ import annotations

from pathlib import Path

import pytest

from core.exceptions import FileValidationError
from core.fasta_sources import (
    combine_fasta_sources,
    discover_fasta_sources,
    prepare_directory_fasta,
)


def write_protein_faa(folder: Path, text: str = ">protein_1\nMSTNPKPQR\n") -> Path:
    """Create an NCBI-like protein.faa file under a source folder."""
    folder.mkdir(parents=True)
    fasta = folder / "protein.faa"
    fasta.write_text(text, encoding="utf-8")
    return fasta


def test_directory_mode_finds_immediate_child_protein_faa(tmp_path: Path) -> None:
    """Immediate child folders with protein.faa should become labeled sources."""
    root = tmp_path / "positive"
    write_protein_faa(root / "organism_A")
    write_protein_faa(root / "organism_B")

    sources, skipped, multiple_file_labels = discover_fasta_sources(root, "positive")

    assert [source.label for source in sources] == ["organism_A", "organism_B"]
    assert [source.path.name for source in sources] == ["protein.faa", "protein.faa"]
    assert skipped == []
    assert multiple_file_labels == []


def test_directory_mode_finds_nested_ncbi_protein_faa(tmp_path: Path) -> None:
    """Nested NCBI protein.faa files should use the source-label folder name."""
    root = tmp_path / "negative"
    nested = (
        root
        / "Sulfolobus_solfataricus"
        / "ncbi_dataset"
        / "data"
        / "GCF_002945325.1"
    )
    write_protein_faa(nested, ">nested_1\nMSTNPKPQR\n")

    sources, skipped, multiple_file_labels = discover_fasta_sources(root, "negative")

    assert len(sources) == 1
    assert sources[0].label == "Sulfolobus_solfataricus"
    assert sources[0].path == nested / "protein.faa"
    assert skipped == []
    assert multiple_file_labels == []


def test_directory_mode_combines_multiple_sources(tmp_path: Path) -> None:
    """Multiple protein.faa files should be written into one combined FASTA."""
    root = tmp_path / "target"
    write_protein_faa(root / "organism_A", ">target_1 first\nMSTNPKPQR\n")
    write_protein_faa(root / "organism_B", ">target_2 second\nMSTNPKPQR\n")

    result = prepare_directory_fasta(
        root,
        "target",
        combined_dir=tmp_path / "combined",
    )

    combined_text = result.combined_fasta.read_text(encoding="utf-8")
    assert result.record_count == 2
    assert result.combined_fasta == tmp_path / "combined" / "target.combined.faa"
    assert result.source_labels == ("organism_A", "organism_B")
    assert ">target_1 first [source=organism_A]" in combined_text
    assert ">target_2 second [source=organism_B]" in combined_text


def test_negative_directory_combines_all_immediate_child_sources(
    tmp_path: Path,
) -> None:
    """All valid negative source-label folders should be scanned and combined."""
    root = tmp_path / "negative"
    source_paths = {
        "Sulfolobus_solfataricus": "GCF_002945325.1",
        "Thermoproteus_tenax": "GCF_xxxxx",
        "Archaeoglobus_fulgidus": "GCF_yyyyy",
    }

    for label, assembly in source_paths.items():
        write_protein_faa(
            root / label / "ncbi_dataset" / "data" / assembly,
            f">{label}_protein\nMSTNPKPQR\n",
        )

    result = prepare_directory_fasta(
        root,
        "negative",
        combined_dir=tmp_path / "combined",
    )

    combined_text = result.combined_fasta.read_text(encoding="utf-8")
    assert result.source_labels == (
        "Archaeoglobus_fulgidus",
        "Sulfolobus_solfataricus",
        "Thermoproteus_tenax",
    )
    assert result.record_count == 3
    for label in source_paths:
        assert f">{label}_protein [source={label}]" in combined_text


def test_missing_protein_faa_child_is_skipped_when_valid_source_exists(
    tmp_path: Path,
) -> None:
    """Folders without protein.faa should be reported without failing the category."""
    root = tmp_path / "negative"
    write_protein_faa(root / "organism_A")
    (root / "missing_protein").mkdir(parents=True)

    sources, skipped, multiple_file_labels = discover_fasta_sources(root, "negative")

    assert [source.label for source in sources] == ["organism_A"]
    assert skipped == ["missing_protein"]
    assert multiple_file_labels == []


def test_no_valid_protein_faa_files_fails_clearly(tmp_path: Path) -> None:
    """A required category with no protein.faa sources should fail clearly."""
    root = tmp_path / "positive"
    (root / "missing_protein").mkdir(parents=True)

    with pytest.raises(FileValidationError, match="No valid protein.faa files"):
        discover_fasta_sources(root, "positive")


def test_duplicate_fasta_ids_are_reported_without_crashing(tmp_path: Path) -> None:
    """Duplicate IDs should be returned to the caller for warning logs."""
    root = tmp_path / "positive"
    write_protein_faa(root / "organism_A", ">shared_id one\nMSTNPKPQR\n")
    write_protein_faa(root / "organism_B", ">shared_id two\nMSTNPKPQR\n")
    sources, _skipped, _multiple_file_labels = discover_fasta_sources(root, "positive")

    _combined, duplicate_ids, count = combine_fasta_sources(
        sources,
        tmp_path / "combined" / "positive.combined.faa",
        "positive",
    )

    assert count == 2
    assert duplicate_ids == ("shared_id",)


def test_multiple_protein_faa_files_under_one_source_label_are_included(
    tmp_path: Path,
) -> None:
    """One source label may contain multiple nested NCBI protein.faa files."""
    root = tmp_path / "negative"
    source = root / "Sulfolobus_solfataricus"
    write_protein_faa(
        source / "ncbi_dataset" / "data" / "GCF_002945325.1",
        ">protein_a\nMSTNPKPQR\n",
    )
    write_protein_faa(
        source / "ncbi_dataset" / "data" / "GCF_999999999.1",
        ">protein_b\nMSTNPKPQR\n",
    )

    result = prepare_directory_fasta(
        root,
        "negative",
        combined_dir=tmp_path / "combined",
    )

    combined_text = result.combined_fasta.read_text(encoding="utf-8")
    assert [source.label for source in result.sources] == [
        "Sulfolobus_solfataricus",
        "Sulfolobus_solfataricus",
    ]
    assert result.source_labels == ("Sulfolobus_solfataricus",)
    assert result.multiple_file_labels == ("Sulfolobus_solfataricus",)
    assert result.record_count == 2
    assert "[source=Sulfolobus_solfataricus]" in combined_text
    assert "GCF_002945325.1]" not in combined_text
    assert "protein.faa]" not in combined_text
