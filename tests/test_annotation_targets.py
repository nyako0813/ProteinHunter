"""Tests for per-sheet annotation target selection."""

from __future__ import annotations

from types import SimpleNamespace

from config import AnnotationTargetConfig, AnnotationTargetsConfig
from core.models import ProteinRecord
from main import _records_for_annotation_step


def targets(
    candidates: AnnotationTargetConfig,
    positive_all_sources: AnnotationTargetConfig,
    no_hit: AnnotationTargetConfig,
    negative_unmatched: AnnotationTargetConfig,
    negative_hit: AnnotationTargetConfig,
) -> AnnotationTargetsConfig:
    """Build annotation target settings for tests."""
    return AnnotationTargetsConfig(
        candidates=candidates,
        positive_all_sources=positive_all_sources,
        no_hit=no_hit,
        negative_unmatched=negative_unmatched,
        negative_hit=negative_hit,
    )


def target(
    gff: bool = True,
    pfam: bool = False,
    uniprot: bool = False,
    alphafold: bool = False,
) -> AnnotationTargetConfig:
    """Build one annotation target setting."""
    return AnnotationTargetConfig(
        gff=gff,
        pfam=pfam,
        uniprot=uniprot,
        alphafold=alphafold,
    )


def test_records_for_annotation_step_selects_no_hit_when_enabled() -> None:
    """No_hit records should be selected when no_hit.pfam is true."""
    candidate = ProteinRecord(protein_id="candidate")
    no_hit = ProteinRecord(protein_id="no_hit")
    classification = SimpleNamespace(
        positive_only_records={"candidate": candidate},
        positive_all_sources_records={},
        no_hit_records={"no_hit": no_hit},
        negative_unmatched_records={"candidate": candidate, "no_hit": no_hit},
        negative_hit_records={},
    )
    annotation_targets = targets(
        candidates=target(pfam=True),
        positive_all_sources=target(pfam=True),
        no_hit=target(pfam=True),
        negative_unmatched=target(),
        negative_hit=target(),
    )

    selected = _records_for_annotation_step(
        classification,
        annotation_targets,
        "pfam",
    )

    assert selected == {
        "candidate": candidate,
        "no_hit": no_hit,
    }


def test_records_for_annotation_step_avoids_duplicate_shared_records() -> None:
    """Records shared by Candidates and Positive_all_sources should appear once."""
    shared = ProteinRecord(protein_id="shared")
    classification = SimpleNamespace(
        positive_only_records={"shared": shared},
        positive_all_sources_records={"shared": shared},
        no_hit_records={},
        negative_unmatched_records={"shared": shared},
        negative_hit_records={},
    )
    annotation_targets = targets(
        candidates=target(pfam=True),
        positive_all_sources=target(pfam=True),
        no_hit=target(),
        negative_unmatched=target(),
        negative_hit=target(),
    )

    selected = _records_for_annotation_step(
        classification,
        annotation_targets,
        "pfam",
    )

    assert list(selected) == ["shared"]
    assert selected["shared"] is shared
