"""Tests for applying UniProt and AlphaFold annotations to records."""

from __future__ import annotations

from unittest.mock import Mock

import pytest

from annotation.record_annotator import (
    annotate_alphafold,
    annotate_records_alphafold,
    annotate_records_uniprot,
    annotate_records_uniprot_and_alphafold,
    annotate_uniprot,
    annotate_uniprot_and_alphafold,
)
from core.exceptions import AlphaFoldAnnotationError, UniProtAnnotationError
from core.models import ProteinRecord


def make_record(protein_id: str = "protein_1") -> ProteinRecord:
    """Create a small ProteinRecord for annotation tests."""
    return ProteinRecord(protein_id=protein_id, sequence="MSTNPKPQR")


def test_successful_uniprot_accession_assignment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """UniProt metadata should set the record accession and annotation result."""
    metadata = {
        "query": "protein_1",
        "accession": "P12345",
        "id": "TEST_PROTEIN",
        "protein_name": "Test protein",
        "organism": "Test organism",
        "reviewed": True,
    }
    monkeypatch.setattr(
        "annotation.record_annotator.search_uniprot_by_protein_id",
        Mock(return_value=metadata),
    )
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        Mock(return_value=None),
    )
    record = make_record()

    result = annotate_uniprot(record, timeout=5)

    assert result is record
    assert record.uniprot_accession == "P12345"
    assert record.annotations["uniprot"].success is True
    assert record.annotations["uniprot"].metadata == metadata


def test_successful_alphafold_url_assignment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """AlphaFold URL should be stored when available."""
    alphafold_mock = Mock(return_value="https://alphafold.ebi.ac.uk/entry/P12345")
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        alphafold_mock,
    )
    record = make_record()
    record.uniprot_accession = "P12345"

    annotate_alphafold(record, timeout=7)

    assert record.alphafold_url == "https://alphafold.ebi.ac.uk/entry/P12345"
    assert record.annotations["alphafold"].success is True
    assert record.annotations["alphafold"].metadata["exists"] is True
    alphafold_mock.assert_called_once_with("P12345", cache=None, timeout=7)


def test_no_uniprot_accession_skips_alphafold(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """AlphaFold should not be checked when UniProt has no accession."""
    monkeypatch.setattr(
        "annotation.record_annotator.search_uniprot_by_protein_id",
        Mock(return_value={"query": "protein_1", "accession": None}),
    )
    alphafold_mock = Mock()
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        alphafold_mock,
    )
    record = make_record()

    annotate_uniprot_and_alphafold(record)

    assert record.uniprot_accession is None
    assert record.alphafold_url is None
    assert record.annotations["alphafold"].success is True
    assert record.annotations["alphafold"].error is None
    assert record.annotations["alphafold"].metadata["exists"] is False
    assert any("AlphaFold annotation skipped" in note for note in record.notes)
    alphafold_mock.assert_not_called()


def test_uniprot_failure_adds_note_and_failed_annotation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """UniProt failures should be recorded without raising."""
    monkeypatch.setattr(
        "annotation.record_annotator.search_uniprot_by_protein_id",
        Mock(side_effect=UniProtAnnotationError("network failed")),
    )
    alphafold_mock = Mock()
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        alphafold_mock,
    )
    record = make_record()

    result = annotate_uniprot(record)

    assert result is record
    assert record.annotations["uniprot"].success is False
    assert "UniProt annotation failed" in str(record.annotations["uniprot"].error)
    assert any("UniProt annotation failed" in note for note in record.notes)
    alphafold_mock.assert_not_called()


def test_alphafold_failure_adds_note_and_failed_annotation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """AlphaFold failures should be recorded without raising."""
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        Mock(side_effect=AlphaFoldAnnotationError("server error")),
    )
    record = make_record()
    record.uniprot_accession = "P12345"

    annotate_alphafold(record)

    assert record.uniprot_accession == "P12345"
    assert record.annotations["alphafold"].success is False
    assert "AlphaFold annotation failed" in str(record.annotations["alphafold"].error)
    assert any("AlphaFold annotation failed" in note for note in record.notes)


def test_annotate_uniprot_no_accession_stores_successful_result(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """UniProt with no accession should still store a successful result."""
    metadata = {"query": "protein_1", "accession": None}
    monkeypatch.setattr(
        "annotation.record_annotator.search_uniprot_by_protein_id",
        Mock(return_value=metadata),
    )
    record = make_record()

    annotate_uniprot(record)

    assert record.uniprot_accession is None
    assert record.annotations["uniprot"].success is True
    assert record.annotations["uniprot"].metadata == metadata


def test_annotate_alphafold_skips_when_accession_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """AlphaFold should be skipped when no UniProt accession is available."""
    alphafold_mock = Mock()
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        alphafold_mock,
    )
    record = make_record()

    annotate_alphafold(record)

    assert record.alphafold_url is None
    assert record.annotations["alphafold"].success is True
    assert record.annotations["alphafold"].metadata == {
        "accession": None,
        "url": None,
        "exists": False,
    }
    assert any("AlphaFold annotation skipped" in note for note in record.notes)
    alphafold_mock.assert_not_called()


def test_annotate_records_uniprot_and_alphafold_annotates_multiple_records(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Batch annotation should annotate each record and return the same dict."""
    metadata_by_id = {
        "protein_1": {"query": "protein_1", "accession": "P11111"},
        "protein_2": {"query": "protein_2", "accession": "P22222"},
    }

    def fake_uniprot(
        protein_id: str,
        *args: object,
        **kwargs: object,
    ) -> dict[str, str]:
        return metadata_by_id[protein_id]

    monkeypatch.setattr(
        "annotation.record_annotator.search_uniprot_by_protein_id",
        fake_uniprot,
    )
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        Mock(return_value=None),
    )
    records = {
        "protein_1": make_record("protein_1"),
        "protein_2": make_record("protein_2"),
    }

    result = annotate_records_uniprot_and_alphafold(records, timeout=9)

    assert result is records
    assert records["protein_1"].uniprot_accession == "P11111"
    assert records["protein_2"].uniprot_accession == "P22222"
    assert records["protein_1"].annotations["uniprot"].success is True
    assert records["protein_2"].annotations["uniprot"].success is True


def test_separate_batch_annotators_return_same_dict(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Separate batch annotators should be usable independently."""
    monkeypatch.setattr(
        "annotation.record_annotator.search_uniprot_by_protein_id",
        Mock(return_value={"query": "protein_1", "accession": "P12345"}),
    )
    monkeypatch.setattr(
        "annotation.record_annotator.get_alphafold_url_if_exists",
        Mock(return_value="https://alphafold.ebi.ac.uk/entry/P12345"),
    )
    records = {"protein_1": make_record("protein_1")}

    uniprot_result = annotate_records_uniprot(records)
    alphafold_result = annotate_records_alphafold(records)

    assert uniprot_result is records
    assert alphafold_result is records
    assert records["protein_1"].uniprot_accession == "P12345"
    assert records["protein_1"].alphafold_url == (
        "https://alphafold.ebi.ac.uk/entry/P12345"
    )
