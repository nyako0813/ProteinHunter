"""Tests for applying CDD domain annotations to records."""

from __future__ import annotations

from unittest.mock import Mock

import pytest

from annotation.domain_annotator import (
    annotate_cdd_domains,
    annotate_pfam_domains,
    annotate_records_cdd,
    annotate_records_pfam,
    filter_pfam_domains,
)
from core.exceptions import CDDAnnotationError, PfamAnnotationError
from core.models import DomainHit, ProteinRecord


def make_record(protein_id: str = "protein_1") -> ProteinRecord:
    """Create a small protein record for CDD annotation tests."""
    return ProteinRecord(protein_id=protein_id, sequence="MSTNPKPQR")


def make_domain(accession: str = "cd12345") -> DomainHit:
    """Create a fake CDD domain hit."""
    return DomainHit(
        source="CDD",
        accession=accession,
        name="Thioredoxin_like",
        description="redox domain",
        evalue=1e-20,
        bitscore=55.5,
        start=10,
        end=80,
    )


def make_pfam_domain(accession: str = "PF00001") -> DomainHit:
    """Create a fake Pfam domain hit."""
    return DomainHit(
        source="Pfam",
        accession=accession,
        name="Pfam_domain",
        description="pfam domain",
        evalue=1e-10,
        bitscore=42.0,
        start=5,
        end=60,
    )


def test_successful_cdd_domain_assignment(monkeypatch: pytest.MonkeyPatch) -> None:
    """CDD hits should be appended to the record and stored as annotation."""
    domain = make_domain()
    search_mock = Mock(return_value=[domain])
    monkeypatch.setattr(
        "annotation.domain_annotator.search_cdd_by_sequence",
        search_mock,
    )
    record = make_record()

    result = annotate_cdd_domains(record, timeout=5)

    assert result is record
    assert record.domains == [domain]
    assert record.annotations["cdd"].success is True
    assert record.annotations["cdd"].domains == [domain]
    assert record.annotations["cdd"].metadata == {"domain_count": 1}
    search_mock.assert_called_once_with(
        protein_id="protein_1",
        sequence="MSTNPKPQR",
        cache=None,
        timeout=5,
    )


def test_no_cdd_domains_stores_successful_empty_annotation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No CDD hits should still produce a successful empty annotation."""
    monkeypatch.setattr(
        "annotation.domain_annotator.search_cdd_by_sequence",
        Mock(return_value=[]),
    )
    record = make_record()

    annotate_cdd_domains(record)

    assert record.domains == []
    assert record.annotations["cdd"].success is True
    assert record.annotations["cdd"].domains == []
    assert record.annotations["cdd"].metadata == {"domain_count": 0}


def test_cdd_failure_adds_note_and_failed_annotation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """CDD failures should be recorded without raising."""
    monkeypatch.setattr(
        "annotation.domain_annotator.search_cdd_by_sequence",
        Mock(side_effect=CDDAnnotationError("CDD service failed")),
    )
    record = make_record()

    result = annotate_cdd_domains(record)

    assert result is record
    assert record.domains == []
    assert record.annotations["cdd"].success is False
    assert "CDD annotation failed" in str(record.annotations["cdd"].error)
    assert any("CDD annotation failed" in note for note in record.notes)


def test_annotate_records_cdd_annotates_multiple_records(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Batch CDD annotation should annotate each record and return the same dict."""
    domains_by_id = {
        "protein_1": [make_domain("cd11111")],
        "protein_2": [make_domain("cd22222")],
    }

    def fake_search(
        protein_id: str,
        *args: object,
        **kwargs: object,
    ) -> list[DomainHit]:
        return domains_by_id[protein_id]

    monkeypatch.setattr(
        "annotation.domain_annotator.search_cdd_by_sequence",
        fake_search,
    )
    records = {
        "protein_1": make_record("protein_1"),
        "protein_2": make_record("protein_2"),
    }

    result = annotate_records_cdd(records, timeout=7)

    assert result is records
    assert records["protein_1"].domains == domains_by_id["protein_1"]
    assert records["protein_2"].domains == domains_by_id["protein_2"]
    assert records["protein_1"].annotations["cdd"].success is True
    assert records["protein_2"].annotations["cdd"].success is True


def test_existing_domains_are_preserved_and_cdd_domains_are_appended(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Existing domains should remain when new CDD domains are added."""
    existing = DomainHit(source="Pfam", accession="PF00001", name="Existing")
    new_domain = make_domain("cd33333")
    monkeypatch.setattr(
        "annotation.domain_annotator.search_cdd_by_sequence",
        Mock(return_value=[new_domain]),
    )
    record = make_record()
    record.domains.append(existing)

    annotate_cdd_domains(record)

    assert record.domains == [existing, new_domain]
    assert record.annotations["cdd"].domains == [new_domain]


def test_successful_pfam_domain_assignment(monkeypatch: pytest.MonkeyPatch) -> None:
    """Pfam hits should be appended to the record and stored as annotation."""
    domain = make_pfam_domain()
    search_mock = Mock(return_value=[domain])
    monkeypatch.setattr(
        "annotation.domain_annotator.search_pfam_by_sequence",
        search_mock,
    )
    record = make_record()

    result = annotate_pfam_domains(record, timeout=5)

    assert result is record
    assert record.domains == [domain]
    assert record.annotations["pfam"].success is True
    assert record.annotations["pfam"].domains == [domain]
    assert record.annotations["pfam"].metadata == {"domain_count": 1}
    search_mock.assert_called_once_with(
        protein_id="protein_1",
        sequence="MSTNPKPQR",
        cache=None,
        timeout=5,
    )


def test_pfam_evalue_filter_keeps_strong_hit(monkeypatch: pytest.MonkeyPatch) -> None:
    """Pfam hits at or below the e-value threshold should be kept."""
    strong = make_pfam_domain("PF00001")
    strong.evalue = 1e-10
    monkeypatch.setattr(
        "annotation.domain_annotator.search_pfam_by_sequence",
        Mock(return_value=[strong]),
    )
    record = make_record()

    annotate_pfam_domains(record, evalue_threshold=1e-5)

    assert record.domains == [strong]
    assert record.annotations["pfam"].metadata == {"domain_count": 1}


def test_pfam_evalue_filter_removes_weak_hit(monkeypatch: pytest.MonkeyPatch) -> None:
    """Pfam hits above the e-value threshold should be filtered out."""
    weak = make_pfam_domain("PF00002")
    weak.evalue = 0.01
    monkeypatch.setattr(
        "annotation.domain_annotator.search_pfam_by_sequence",
        Mock(return_value=[weak]),
    )
    record = make_record()

    annotate_pfam_domains(record, evalue_threshold=1e-5)

    assert record.domains == []
    assert record.annotations["pfam"].domains == []
    assert record.annotations["pfam"].metadata == {"domain_count": 0}


def test_pfam_evalue_filter_keeps_missing_evalue() -> None:
    """Pfam hits without e-values are kept for now."""
    missing = make_pfam_domain("PF00003")
    missing.evalue = None

    assert filter_pfam_domains([missing], evalue_threshold=1e-5) == [missing]


def test_no_pfam_domains_stores_successful_empty_annotation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No Pfam hits should still produce a successful empty annotation."""
    monkeypatch.setattr(
        "annotation.domain_annotator.search_pfam_by_sequence",
        Mock(return_value=[]),
    )
    record = make_record()

    annotate_pfam_domains(record)

    assert record.domains == []
    assert record.annotations["pfam"].success is True
    assert record.annotations["pfam"].domains == []
    assert record.annotations["pfam"].metadata == {"domain_count": 0}


def test_pfam_failure_adds_note_and_failed_annotation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Pfam failures should be recorded without raising."""
    monkeypatch.setattr(
        "annotation.domain_annotator.search_pfam_by_sequence",
        Mock(side_effect=PfamAnnotationError("Pfam service failed")),
    )
    record = make_record()

    result = annotate_pfam_domains(record)

    assert result is record
    assert record.domains == []
    assert record.annotations["pfam"].success is False
    assert "Pfam annotation failed" in str(record.annotations["pfam"].error)
    assert any("Pfam annotation failed" in note for note in record.notes)


def test_annotate_records_pfam_annotates_multiple_records(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Batch Pfam annotation should annotate each record and return the same dict."""
    domains_by_id = {
        "protein_1": [make_pfam_domain("PF00001")],
        "protein_2": [make_pfam_domain("PF00002")],
    }

    def fake_search(
        protein_id: str,
        *args: object,
        **kwargs: object,
    ) -> list[DomainHit]:
        return domains_by_id[protein_id]

    monkeypatch.setattr(
        "annotation.domain_annotator.search_pfam_by_sequence",
        fake_search,
    )
    records = {
        "protein_1": make_record("protein_1"),
        "protein_2": make_record("protein_2"),
    }

    result = annotate_records_pfam(records, timeout=7)

    assert result is records
    assert records["protein_1"].domains == domains_by_id["protein_1"]
    assert records["protein_2"].domains == domains_by_id["protein_2"]
    assert records["protein_1"].annotations["pfam"].success is True
    assert records["protein_2"].annotations["pfam"].success is True


def test_existing_domains_are_preserved_and_pfam_domains_are_appended(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Existing domains should remain when new Pfam domains are added."""
    existing = DomainHit(source="CDD", accession="cd12345", name="Existing")
    new_domain = make_pfam_domain("PF00003")
    monkeypatch.setattr(
        "annotation.domain_annotator.search_pfam_by_sequence",
        Mock(return_value=[new_domain]),
    )
    record = make_record()
    record.domains.append(existing)

    annotate_pfam_domains(record)

    assert record.domains == [existing, new_domain]
    assert record.annotations["pfam"].domains == [new_domain]
