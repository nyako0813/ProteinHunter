"""Tests for the evidence model (core/evidence.py)."""

from __future__ import annotations

import pytest

from core.evidence import (
    EvidenceComponent,
    EvidenceStatus,
    NON_CONTRIBUTING_STATUSES,
    clamp01,
    linear_normalize,
)


def test_clamp01_bounds() -> None:
    assert clamp01(-0.5) == 0.0
    assert clamp01(1.5) == 1.0
    assert clamp01(0.3) == 0.3


def test_linear_normalize_basic() -> None:
    assert linear_normalize(0, 0, 100) == 0.0
    assert linear_normalize(100, 0, 100) == 1.0
    assert linear_normalize(50, 0, 100) == 0.5
    assert linear_normalize(-10, 0, 100) == 0.0
    assert linear_normalize(200, 0, 100) == 1.0


def test_linear_normalize_inverted_scale() -> None:
    # smaller distance -> higher support
    assert linear_normalize(0, 100, 0) == 1.0
    assert linear_normalize(100, 100, 0) == 0.0
    assert linear_normalize(25, 100, 0) == 0.75


def test_linear_normalize_rejects_equal_bounds() -> None:
    with pytest.raises(ValueError):
        linear_normalize(1, 5, 5)


def test_available_component_contribution() -> None:
    component = EvidenceComponent.available(
        "domain_pair", "functional_annotation", normalized_value=0.8, weight=20.0
    )
    assert component.status is EvidenceStatus.AVAILABLE
    assert component.effective_weight == 20.0
    assert component.contribution == pytest.approx(16.0)


def test_available_component_clamps_out_of_range_value() -> None:
    component = EvidenceComponent.available(
        "domain_pair", "functional_annotation", normalized_value=1.4, weight=10.0
    )
    assert component.normalized_value == 1.0


@pytest.mark.parametrize(
    "status",
    [
        EvidenceStatus.MISSING,
        EvidenceStatus.NOT_RUN,
        EvidenceStatus.NOT_APPLICABLE,
        EvidenceStatus.FAILED,
        EvidenceStatus.MALFORMED,
        EvidenceStatus.EXCLUDED,
    ],
)
def test_unavailable_component_never_contributes(status: EvidenceStatus) -> None:
    component = EvidenceComponent.unavailable("gene_neighborhood", "genomic_context", status)
    assert component.effective_weight == 0.0
    assert component.contribution == 0.0
    assert component.normalized_value is None
    assert status in NON_CONTRIBUTING_STATUSES


def test_unavailable_rejects_available_status() -> None:
    with pytest.raises(ValueError):
        EvidenceComponent.unavailable("x", "y", EvidenceStatus.AVAILABLE)


def test_available_construction_requires_normalized_value() -> None:
    with pytest.raises(ValueError):
        EvidenceComponent(
            name="x",
            category="y",
            status=EvidenceStatus.AVAILABLE,
            normalized_value=None,
        )


def test_available_construction_rejects_out_of_range_direct_init() -> None:
    with pytest.raises(ValueError):
        EvidenceComponent(
            name="x",
            category="y",
            status=EvidenceStatus.AVAILABLE,
            normalized_value=1.5,
        )


def test_available_construction_rejects_negative_weight() -> None:
    with pytest.raises(ValueError):
        EvidenceComponent(
            name="x",
            category="y",
            status=EvidenceStatus.AVAILABLE,
            normalized_value=0.5,
            weight=-1.0,
        )


def test_negative_flagged_component_still_reports_contribution() -> None:
    # The scoring engine, not the component, decides how a negative-flagged
    # contribution is applied (as a penalty). The component itself only
    # reports its own weighted value.
    component = EvidenceComponent.available(
        "incompatible_localization",
        "cellular_compatibility",
        normalized_value=1.0,
        weight=15.0,
        is_negative=True,
    )
    assert component.is_negative is True
    assert component.contribution == pytest.approx(15.0)
