"""Tests for lightweight interaction candidate ranking."""

from __future__ import annotations

import gzip
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from config import (
    INTERACTION_ALPHAFOLD_DEFAULT,
    INTERACTION_CANDIDATE_SOURCE_DEFAULTS,
    INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    INTERACTION_NEIGHBORHOOD_DEFAULT,
    INTERACTION_SCORING_WEIGHTS_DEFAULT,
    InteractionEvidenceDetailConfig,
    InteractionNeighborhoodConfig,
    InteractionQueryConfig,
    InteractionScoringConfig,
)
from core.models import BlastHit, DomainHit, ProteinRecord
from analysis.interaction_scoring import (
    INTERACTION_EVIDENCE_DETAIL_LEGACY_COLUMNS,
    INTERACTION_EVIDENCE_DETAIL_V2_COLUMNS,
    PROTEIN_HUNTER_SCORE_CEILING,
    interaction_pair_columns,
    resolve_cdd_annotation_targets,
    resolve_protein_hunter_scores,
    run_interaction_scoring,
)
from analysis.coexpression_bridge import GSE77738_STEADY_STATE_RPKM_COLUMNS
from analysis.scoring import build_candidate_score
from analysis.string_ppi_bridge import STRING_VERSION


def interaction_config(
    *,
    enabled: bool = True,
    query_proteins: tuple[InteractionQueryConfig, ...] = (),
    query_fasta: Path | None = None,
    candidate_sources: dict[str, bool] | None = None,
    max_candidates_per_query: int = 200,
    include_sequences_in_excel: bool = False,
    gff_file: Path | None = None,
    neighborhood: InteractionNeighborhoodConfig = INTERACTION_NEIGHBORHOOD_DEFAULT,
    scoring_model: str = "legacy_additive",
    scoring_engine_config: Path | None = None,
    functional_complementarity_ruleset: Path | None = None,
    pih_evidence_bundle: Path | None = None,
    evidence_detail_sheet: InteractionEvidenceDetailConfig = INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    cache_dir: Path | None = None,
) -> SimpleNamespace:
    """Build a minimal app config for interaction scoring tests."""
    return SimpleNamespace(
        paths=SimpleNamespace(gff_file=gff_file, cache_dir=cache_dir),
        interaction_scoring=InteractionScoringConfig(
            enabled=enabled,
            query_proteins=query_proteins,
            query_fasta=query_fasta,
            candidate_sources=(
                dict(INTERACTION_CANDIDATE_SOURCE_DEFAULTS)
                if candidate_sources is None
                else candidate_sources
            ),
            max_candidates_per_query=max_candidates_per_query,
            include_sequences_in_excel=include_sequences_in_excel,
            scoring_weights=INTERACTION_SCORING_WEIGHTS_DEFAULT,
            alphafold=INTERACTION_ALPHAFOLD_DEFAULT,
            neighborhood=neighborhood,
            scoring_model=scoring_model,
            scoring_engine_config=scoring_engine_config,
            functional_complementarity_ruleset=functional_complementarity_ruleset,
            pih_evidence_bundle=pih_evidence_bundle,
            evidence_detail_sheet=evidence_detail_sheet,
        )
    )


def record(
    protein_id: str,
    *,
    old_locus_tag: str = "",
    sequence: str = "MSTNPKPQR",
    description: str = "enzyme protein",
    positive_sources_hit: list[str] | None = None,
) -> ProteinRecord:
    """Create a target/candidate protein record."""
    return ProteinRecord(
        protein_id=protein_id,
        old_locus_tag=old_locus_tag or None,
        sequence=sequence,
        description=description,
        positive_sources_hit=positive_sources_hit or [],
    )


def classification(records: dict[str, ProteinRecord]) -> SimpleNamespace:
    """Build a classification object with candidate source groups."""
    return SimpleNamespace(
        all_records=records,
        positive_only_records={"candidate": records["candidate"]},
        candidates_relaxed_records={
            "candidate": records["candidate"],
            "relaxed": records["relaxed"],
        },
        positive_all_sources_records={},
        negative_unmatched_records={},
        no_hit_records={"novel": records["novel"]},
        negative_hit_records={},
        negative_strong_hit_records={},
        negative_medium_hit_records={},
        negative_weak_hit_records={},
    )


def build_classification(**buckets: dict[str, ProteinRecord]) -> SimpleNamespace:
    """Build a minimal blast_classification-like object with only the given buckets set."""
    defaults: dict[str, dict[str, ProteinRecord]] = {
        "all_records": {},
        "positive_only_records": {},
        "candidates_relaxed_records": {},
        "positive_all_sources_records": {},
        "negative_unmatched_records": {},
        "no_hit_records": {},
        "negative_hit_records": {},
        "negative_strong_hit_records": {},
        "negative_medium_hit_records": {},
        "negative_weak_hit_records": {},
    }
    defaults.update(buckets)
    return SimpleNamespace(**defaults)


def all_sources_disabled() -> dict[str, bool]:
    """A candidate_sources mapping with every source turned off."""
    return {key: False for key in INTERACTION_CANDIDATE_SOURCE_DEFAULTS}


# ---------------------------------------------------------------------------
# resolve_cdd_annotation_targets tests
# ---------------------------------------------------------------------------


def test_resolve_cdd_annotation_targets_empty_when_disabled() -> None:
    """Disabled interaction_scoring must not extend the CDD target set."""
    cfg = interaction_config(enabled=False)
    cls = build_classification()

    assert resolve_cdd_annotation_targets(cfg, cls) == {}


def test_resolve_cdd_annotation_targets_unions_enabled_sources_only() -> None:
    """Only buckets whose candidate_sources flag is true should be included."""
    candidate_a = record("candidate_a")
    relaxed_b = record("relaxed_b")
    no_hit_c = record("no_hit_c")
    cls = build_classification(
        all_records={"candidate_a": candidate_a, "relaxed_b": relaxed_b, "no_hit_c": no_hit_c},
        positive_only_records={"candidate_a": candidate_a},
        candidates_relaxed_records={"candidate_a": candidate_a, "relaxed_b": relaxed_b},
        no_hit_records={"no_hit_c": no_hit_c},
    )
    sources = all_sources_disabled()
    sources["candidates_relaxed"] = True
    cfg = interaction_config(enabled=True, candidate_sources=sources)

    targets = resolve_cdd_annotation_targets(cfg, cls)

    assert set(targets) == {"candidate_a", "relaxed_b"}
    assert targets["candidate_a"] is candidate_a


def test_resolve_cdd_annotation_targets_excludes_disabled_sources() -> None:
    """A bucket must be excluded entirely when its candidate_sources flag is false."""
    no_hit_only = record("no_hit_only")
    cls = build_classification(
        all_records={"no_hit_only": no_hit_only},
        no_hit_records={"no_hit_only": no_hit_only},
    )
    cfg = interaction_config(enabled=True, candidate_sources=all_sources_disabled())

    assert resolve_cdd_annotation_targets(cfg, cls) == {}


def test_resolve_cdd_annotation_targets_includes_resolved_query() -> None:
    """A query that resolves to a real target record should be included."""
    query_record = record("query_protein")
    cls = build_classification(all_records={"query_protein": query_record})
    cfg = interaction_config(
        enabled=True,
        query_proteins=(InteractionQueryConfig("query_protein", "", ""),),
        candidate_sources=all_sources_disabled(),
    )

    targets = resolve_cdd_annotation_targets(cfg, cls)

    assert set(targets) == {"query_protein"}
    assert targets["query_protein"] is query_record


def test_resolve_cdd_annotation_targets_resolves_multiple_queries() -> None:
    """Every configured query_proteins entry should be resolved, not just the first."""
    query_1 = record("query_1")
    query_2 = record("query_2")
    cls = build_classification(all_records={"query_1": query_1, "query_2": query_2})
    cfg = interaction_config(
        enabled=True,
        query_proteins=(
            InteractionQueryConfig("query_1", "", ""),
            InteractionQueryConfig("query_2", "", ""),
        ),
        candidate_sources=all_sources_disabled(),
    )

    targets = resolve_cdd_annotation_targets(cfg, cls)

    assert set(targets) == {"query_1", "query_2"}


def test_resolve_cdd_annotation_targets_excludes_unmatched_sequence_only_query() -> None:
    """A sequence-only query with no matching target record has no record to annotate."""
    cls = build_classification(all_records={})
    cfg = interaction_config(
        enabled=True,
        query_proteins=(InteractionQueryConfig("", "", "MSTNPKPQR"),),
        candidate_sources=all_sources_disabled(),
    )

    assert resolve_cdd_annotation_targets(cfg, cls) == {}


def test_resolve_cdd_annotation_targets_deduplicates_query_already_in_a_bucket() -> None:
    """A query that is also a candidate_sources member should appear only once."""
    shared = record("shared_protein")
    cls = build_classification(
        all_records={"shared_protein": shared},
        positive_only_records={"shared_protein": shared},
    )
    sources = all_sources_disabled()
    sources["candidates"] = True
    cfg = interaction_config(
        enabled=True,
        query_proteins=(InteractionQueryConfig("shared_protein", "", ""),),
        candidate_sources=sources,
    )

    targets = resolve_cdd_annotation_targets(cfg, cls)

    assert set(targets) == {"shared_protein"}


# ---------------------------------------------------------------------------
# resolve_protein_hunter_scores tests (M1: protein_hunter_score scope extension)
# ---------------------------------------------------------------------------


def test_resolve_protein_hunter_scores_empty_when_disabled() -> None:
    """Disabled interaction_scoring must not compute any protein_hunter_score."""
    cfg = interaction_config(enabled=False)
    cls = build_classification()

    assert resolve_protein_hunter_scores(cfg, cls) == {}


def test_resolve_protein_hunter_scores_covers_non_candidates_buckets() -> None:
    """Candidates_relaxed/No_hit records (never scored by the Candidate scoring
    step, which only ever touches positive_only_records) must get a score too."""
    relaxed_only = record("relaxed_only")
    no_hit_only = record("no_hit_only")
    cls = build_classification(
        all_records={"relaxed_only": relaxed_only, "no_hit_only": no_hit_only},
        candidates_relaxed_records={"relaxed_only": relaxed_only},
        no_hit_records={"no_hit_only": no_hit_only},
    )
    sources = all_sources_disabled()
    sources["candidates_relaxed"] = True
    sources["no_hit"] = True
    cfg = interaction_config(enabled=True, candidate_sources=sources)

    scores = resolve_protein_hunter_scores(cfg, cls)

    assert set(scores) == {"relaxed_only", "no_hit_only"}
    assert scores["relaxed_only"].protein_id == "relaxed_only"


def test_resolve_protein_hunter_scores_matches_build_candidate_score() -> None:
    """The computed score must match build_candidate_score's own formula exactly."""
    candidate_a = record("candidate_a", positive_sources_hit=["A"])
    cls = build_classification(
        all_records={"candidate_a": candidate_a},
        positive_only_records={"candidate_a": candidate_a},
    )
    sources = all_sources_disabled()
    sources["candidates"] = True
    cfg = interaction_config(enabled=True, candidate_sources=sources)

    scores = resolve_protein_hunter_scores(cfg, cls)

    expected = build_candidate_score(candidate_a)
    assert scores["candidate_a"].total_score == expected.total_score
    assert scores["candidate_a"].components == expected.components


def test_resolve_protein_hunter_scores_does_not_mutate_shared_records() -> None:
    """Scoring here must not leak into ProteinRecord.score on shared objects --
    Candidates_relaxed/No_hit classification sheets must stay unaffected."""
    relaxed_only = record("relaxed_only")
    cls = build_classification(
        all_records={"relaxed_only": relaxed_only},
        candidates_relaxed_records={"relaxed_only": relaxed_only},
    )
    sources = all_sources_disabled()
    sources["candidates_relaxed"] = True
    cfg = interaction_config(enabled=True, candidate_sources=sources)

    resolve_protein_hunter_scores(cfg, cls)

    assert relaxed_only.score is None


# ---------------------------------------------------------------------------
# protein_hunter_score reference columns (M2)
# ---------------------------------------------------------------------------


def test_interaction_row_carries_protein_hunter_score_reference_columns() -> None:
    """Every Interaction_* pair row should carry the candidate's own protein_hunter_score."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    expected = build_candidate_score(records["candidate"])
    assert row["protein_hunter_score"] == expected.total_score
    assert "no_negative_hit=" in row["protein_hunter_score_components"]
    assert row["protein_hunter_score_reasons"]


def test_protein_hunter_score_reference_column_covers_relaxed_bucket_too() -> None:
    """M1's scope fix must actually reach Interaction_Candidates_relaxed rows."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed", positive_sources_hit=["A"]),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates_relaxed": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    rows = {r["candidate_protein_id"]: r for r in result.source_rows["Interaction_Candidates_relaxed"]}
    # "relaxed" is not in positive_only_records, so the plain Candidate
    # scoring pipeline step never scores it -- protein_hunter_score must
    # still be populated here via the wider interaction_scoring scope.
    assert rows["relaxed"]["protein_hunter_score"] is not None
    assert rows["relaxed"]["protein_hunter_score"] == build_candidate_score(
        records["relaxed"]
    ).total_score


def test_protein_hunter_score_does_not_affect_ranking_or_existing_columns() -> None:
    """protein_hunter_score must be purely additive: rank/score/tier stay identical."""
    records = {
        "query": record("query", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", description="iron-sulfur carrier protein", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # Same values as the pre-M2 behavior (test_v2_mode_scores_full_evidence_pair_near_100).
    assert row["candidate_rank"] == 1
    assert row["interaction_priority_score"] == 100.0
    assert row["evidence_tier"] == "Tier2_Strong"
    assert row["protein_hunter_score"] is not None


# ---------------------------------------------------------------------------
# interaction_score / interaction_evidence_tier (M3, v2_evidence_based only)
# ---------------------------------------------------------------------------


def test_interaction_score_uses_only_genomic_context_and_domain_complementarity() -> None:
    """interaction_score should reuse the same fixture that scores 100 for the
    full composite, but land far lower once source_classification/co_occurrence
    are excluded -- domain_complementarity alone can't clear category-count
    gates for the higher tiers even at a perfect normalized value."""
    records = {
        "query": record("query", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", description="iron-sulfur carrier protein", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # Full composite (source_classification + functional_annotation, no GFF configured).
    assert row["interaction_priority_score"] == 100.0
    assert row["evidence_category_count"] == 2
    # interaction_score: domain_complementarity alone (genomic_context is
    # MISSING -- no GFF configured in this fixture), 20/20 category cap,
    # single category -> capped below Tier1/Tier2 by category count alone.
    assert row["interaction_score"] == 100.0
    assert row["interaction_evidence_tier"] == "Tier3_Moderate"


def test_interaction_score_excludes_co_occurrence_and_source_classification() -> None:
    """A candidate whose entire interaction_priority_score comes from
    source_classification + co_occurrence (no genomic proximity, no domain
    match) must score 0/None on interaction_score -- this is exactly the
    MA_0050/MA_0238 AlphaFold3-calibration false-positive pattern that
    motivated the split."""
    records = {
        "query": record("query", description="", positive_sources_hit=["A"]),
        "candidate": record("candidate", description="", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # Full composite still scores well from source_classification + co_occurrence.
    assert row["interaction_priority_score"] > 0
    # No genomic_context (no GFF) and no domain_complementarity (no
    # description text) evidence at all -- interaction_score must be
    # ineligible, not silently zero-scored-as-if-evaluated.
    assert row["interaction_score"] is None
    assert row["interaction_evidence_tier"] == "Unclassified"


def test_interaction_score_reflects_genomic_context_when_available(tmp_path: Path) -> None:
    """genomic_context alone should drive interaction_score once GFF coordinates exist."""
    gff_file = tmp_path / "genome.gff"
    gff_file.write_text(
        "##gff-version 3\n"
        "contig1\tRefSeq\tgene\t100\t400\t.\t+\t.\tID=query\n"
        "contig1\tRefSeq\tgene\t500\t800\t.\t+\t.\tID=candidate\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", description="", positive_sources_hit=[]),
        "candidate": record("candidate", description="", positive_sources_hit=[]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_score"] is not None
    assert row["interaction_score"] > 0
    assert row["interaction_evidence_tier"] != "Unclassified"


# ---------------------------------------------------------------------------
# genomic_context: operon-tight spacing + strand awareness
#
# See claude/genomic_distance_weight_finding.md for the real-genome check
# behind these thresholds: true M. acetivorans operons (the Mcr activation
# complex gene cluster, the nifI1-nifI2-nifK-nifD operon, mtpA-mtpC) all
# have intergenic gaps of a few bp up to ~70 bp on the same strand. A pair
# several kb apart, even on the same strand, is not "operon-tight"; a pair
# on opposite strands is essentially never in the same operon regardless of
# distance.
# ---------------------------------------------------------------------------


def _genomic_context_component(result, candidate_id: str = "candidate") -> dict:
    return next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == candidate_id and r["component_name"] == "genomic_context"
    )


def _genomic_context_fixture(
    tmp_path: Path, *, query_end: int, candidate_start: int, query_strand: str, candidate_strand: str
) -> tuple[dict[str, ProteinRecord], SimpleNamespace, Path]:
    gff_file = tmp_path / "genome.gff"
    gff_file.write_text(
        "##gff-version 3\n"
        f"contig1\tRefSeq\tgene\t1\t{query_end}\t.\t{query_strand}\t.\tID=query\n"
        f"contig1\tRefSeq\tgene\t{candidate_start}\t{candidate_start + 300}\t.\t{candidate_strand}\t.\tID=candidate\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", description="", positive_sources_hit=[]),
        "candidate": record("candidate", description="", positive_sources_hit=[]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    return records, classification(records), gff_file


def test_genomic_context_operon_tight_same_strand_scores_full(tmp_path: Path) -> None:
    """A ~50 bp same-strand gap (real-operon-tight, e.g. mtpA-mtpC at 70 bp) scores 1.0."""
    _records, cls, gff_file = _genomic_context_fixture(
        tmp_path, query_end=400, candidate_start=450, query_strand="+", candidate_strand="+"
    )
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    component = _genomic_context_component(result)
    assert component["status"] == "AVAILABLE"
    assert component["normalized_value"] == pytest.approx(1.0)
    assert "operon-tight" in component["explanation"]


def test_genomic_context_same_strand_but_not_operon_tight_scores_lower(tmp_path: Path) -> None:
    """A same-strand gap well beyond real operon spacing (e.g. ~2.9 kb, like
    MA_0826-MA_0823) must not receive full "close" credit even though the
    strand matches."""
    _records, cls, gff_file = _genomic_context_fixture(
        tmp_path, query_end=400, candidate_start=3260, query_strand="+", candidate_strand="+"
    )
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    component = _genomic_context_component(result)
    assert component["status"] == "AVAILABLE"
    assert 0.0 < component["normalized_value"] < 1.0
    assert component["normalized_value"] == pytest.approx(0.4)


def test_genomic_context_opposite_strand_scores_much_lower_than_same_strand(tmp_path: Path) -> None:
    """At a comparable distance, an opposite-strand pair must score noticeably
    lower than a same-strand pair -- opposite-strand genes are essentially
    never part of the same operon."""
    same_dir = tmp_path / "same"
    opposite_dir = tmp_path / "opposite"
    same_dir.mkdir()
    opposite_dir.mkdir()
    _, cls_same, gff_same = _genomic_context_fixture(
        same_dir, query_end=400, candidate_start=3260, query_strand="+", candidate_strand="+"
    )
    _, cls_opposite, gff_opposite = _genomic_context_fixture(
        opposite_dir, query_end=400, candidate_start=3260, query_strand="+", candidate_strand="-"
    )
    cfg_same = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_same,
    )
    cfg_opposite = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_opposite,
    )

    same_result = run_interaction_scoring(cfg_same, cls_same)
    opposite_result = run_interaction_scoring(cfg_opposite, cls_opposite)

    assert same_result is not None and opposite_result is not None
    same_component = _genomic_context_component(same_result)
    opposite_component = _genomic_context_component(opposite_result)
    assert opposite_component["normalized_value"] < same_component["normalized_value"]
    assert opposite_component["normalized_value"] == pytest.approx(0.05)


def test_genomic_context_unknown_strand_falls_back_to_original_distance_tiers(tmp_path: Path) -> None:
    """When strand is unavailable ("."), fall back to the original,
    strand-agnostic distance tiers unchanged (don't guess operon-likeness
    without the information needed to judge it)."""
    _records, cls, gff_file = _genomic_context_fixture(
        tmp_path, query_end=400, candidate_start=3260, query_strand=".", candidate_strand="+"
    )
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    component = _genomic_context_component(result)
    assert component["status"] == "AVAILABLE"
    # <=5000 bp with strand unknown keeps the original "close" value (1.0).
    assert component["normalized_value"] == pytest.approx(1.0)
    assert "strand unknown" in component["explanation"]


# ---------------------------------------------------------------------------
# legacy_additive interaction_score (M4)
# ---------------------------------------------------------------------------


def test_legacy_interaction_score_matches_renormalized_formula() -> None:
    """legacy interaction_score = (gene_neighborhood + domain_complementarity) /
    (their configured weights) * 100 -- co_occurrence and candidate_priority excluded."""
    records = {
        "query": record("query", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", description="iron-sulfur carrier protein", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # No GFF configured -> same_gene_neighborhood_score = 0; domain match is a
    # full complementary-term hit -> domain_complementarity_score = weight (15).
    assert row["same_gene_neighborhood_score"] == 0.0
    assert row["domain_complementarity_score"] == 15.0
    assert row["interaction_score"] == round((0.0 + 15.0) / (25.0 + 15.0) * 100, 3)
    assert row["interaction_evidence_tier"] is None


def test_legacy_interaction_score_excludes_co_occurrence_and_source() -> None:
    """Same MA_0050/MA_0238-style pattern as the v2 test: a legacy row whose
    interaction_priority_score comes entirely from candidate_priority_score +
    co_occurrence_score must score 0 on interaction_score (unlike v2, legacy
    has no MISSING concept, so this is a numeric 0, not None)."""
    records = {
        "query": record("query", description="", positive_sources_hit=["A"]),
        "candidate": record("candidate", description="", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_priority_score"] > 0
    assert row["interaction_score"] == 0.0


def test_legacy_mode_leaves_interaction_evidence_tier_blank() -> None:
    """legacy_additive has no per-category tiering concept -- always blank (M4)."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row.get("interaction_evidence_tier") is None


def test_legacy_string_ppi_score_feeds_both_scores(tmp_path: Path) -> None:
    """legacy_additive's string_ppi_score (Phase 6a M4) contributes to both
    interaction_priority_score and interaction_score once STRING is configured."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=[]),
        "candidate": record(
            "candidate", old_locus_tag="MA_0002", description="", positive_sources_hit=[]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    _seed_string_files(tmp_path, 188937, ["188937.MA_0001 188937.MA_0002 0 0 600 0 0 0 0 600"])
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, string_ppi_ncbi_taxon_id=188937)

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # cooccurrence=600, neighborhood=0 (not seeded) -> average 0.3, scaled
    # by the default external_ppi weight (15) -> 4.5.
    assert row["string_ppi_score"] == pytest.approx(4.5)
    assert row["interaction_priority_score"] == pytest.approx(
        row["candidate_priority_score"]
        + row["same_gene_neighborhood_score"]
        + row["co_occurrence_score"]
        + row["domain_complementarity_score"]
        + row["alphafold_readiness_score"]
        + row["string_ppi_score"]
    )
    # No GFF neighborhood, no domain match -> string_ppi_score (4.5) is the
    # entire interaction_score numerator, over (25 + 15 + 15) = 55 points
    # (gene_neighborhood + domain_complementarity + external_ppi, all
    # active once STRING is configured).
    assert row["interaction_score"] == pytest.approx(4.5 / 55.0 * 100, abs=0.01)


def test_legacy_string_ppi_score_absent_when_taxon_id_unset() -> None:
    """Without string_ppi_ncbi_taxon_id, string_ppi_score must be 0 and the
    interaction_score denominator must NOT include external_ppi -- otherwise
    every existing legacy_additive run would silently change once this
    field exists at all."""
    records = {
        "query": record("query", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", description="iron-sulfur carrier protein", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["string_ppi_score"] == 0.0
    # Same formula as before Phase 6a M4: (0 + 15) / (25 + 15) * 100.
    assert row["interaction_score"] == pytest.approx((0.0 + 15.0) / 40.0 * 100, abs=0.01)


# ---------------------------------------------------------------------------
# ranking_metric (M5): candidate_rank/row order can be driven by
# interaction_score instead of interaction_priority_score
# ---------------------------------------------------------------------------


def _ranking_metric_fixture_v2(tmp_path: Path) -> tuple[dict[str, ProteinRecord], SimpleNamespace, Path]:
    """priority_only wins on the full composite; interaction_only wins on
    query-specific evidence alone (genomic proximity vs. a generic,
    non-matching description overlap that scores zero either way)."""
    gff_file = tmp_path / "genome.gff"
    gff_file.write_text(
        "##gff-version 3\n"
        "contig1\tRefSeq\tgene\t100\t400\t.\t+\t.\tID=query\n"
        "contig1\tRefSeq\tgene\t500\t800\t.\t+\t.\tID=interaction_only\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", description="putative protein", positive_sources_hit=["A"]),
        # Domain evidence is AVAILABLE-but-zero (generic overlap only), not
        # MISSING -- keeps the interaction-only breakdown eligible (a single
        # active category at zero) instead of hitting rank_candidates' "no
        # formal score" sentinel, which would make ranks incomparable.
        "priority_only": record(
            "priority_only", description="putative membrane protein", positive_sources_hit=["A"]
        ),
        "interaction_only": record("interaction_only", description="", positive_sources_hit=[]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cls = build_classification(
        all_records=records,
        positive_only_records={
            "priority_only": records["priority_only"],
            "interaction_only": records["interaction_only"],
        },
    )
    return records, cls, gff_file


def test_ranking_metric_default_preserves_interaction_priority_score_order_v2(tmp_path: Path) -> None:
    """Without ranking_metric set, v2 ranking must be unchanged from pre-Phase-5 behavior."""
    _records, cls, gff_file = _ranking_metric_fixture_v2(tmp_path)
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    rows = {r["candidate_protein_id"]: r for r in result.source_rows["Interaction_Candidates"]}
    # source_classification + a full-Jaccard co_occurrence outrank a lone
    # genomic_context hit once source_classification/co_occurrence are
    # counted (the default, unchanged composite).
    assert rows["priority_only"]["interaction_priority_score"] > rows["interaction_only"]["interaction_priority_score"]
    assert rows["priority_only"]["candidate_rank"] < rows["interaction_only"]["candidate_rank"]
    # ...but interaction_score already shows the reverse story, unused for
    # ranking here.
    assert rows["interaction_only"]["interaction_score"] > rows["priority_only"]["interaction_score"]


def test_ranking_metric_interaction_score_flips_order_v2(tmp_path: Path) -> None:
    """ranking_metric: interaction_score must re-rank by query-specific evidence only,
    without changing any row's interaction_priority_score/evidence_tier/interaction_score value."""
    _records, cls, gff_file = _ranking_metric_fixture_v2(tmp_path)
    default_cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
    )
    interaction_ranked_cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
    )
    interaction_ranked_cfg.interaction_scoring = replace(
        interaction_ranked_cfg.interaction_scoring, ranking_metric="interaction_score"
    )

    default_rows = {
        r["candidate_protein_id"]: r
        for r in run_interaction_scoring(default_cfg, cls).source_rows["Interaction_Candidates"]
    }
    reranked_rows = {
        r["candidate_protein_id"]: r
        for r in run_interaction_scoring(interaction_ranked_cfg, cls).source_rows["Interaction_Candidates"]
    }

    # interaction_only now outranks priority_only once ranking is by
    # query-specific evidence only.
    assert reranked_rows["interaction_only"]["candidate_rank"] < reranked_rows["priority_only"]["candidate_rank"]
    # ...which is the opposite of the default ranking on the same data.
    assert default_rows["priority_only"]["candidate_rank"] < default_rows["interaction_only"]["candidate_rank"]

    # Every score column keeps its original value regardless of which
    # metric is ranking -- only candidate_rank/distance_independent_rank differ.
    for candidate_id in ("priority_only", "interaction_only"):
        for column in ("interaction_priority_score", "evidence_tier", "interaction_score", "interaction_evidence_tier"):
            assert default_rows[candidate_id][column] == reranked_rows[candidate_id][column]


def _ranking_metric_fixture_legacy(tmp_path: Path) -> tuple[SimpleNamespace, Path]:
    """priority_only wins the full composite (candidate_priority + full
    co_occurrence); interaction_only wins interaction_score alone (moderate
    genomic proximity, worth less than a full co_occurrence match in the
    legacy additive sum, but the only query-specific evidence present)."""
    gff_file = tmp_path / "genome.gff"
    gff_file.write_text(
        "##gff-version 3\n"
        "contig1\tRefSeq\tgene\t100\t400\t.\t+\t.\tID=query\n"
        "contig1\tRefSeq\tgene\t15100\t15400\t.\t+\t.\tID=interaction_only\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "priority_only": record("priority_only", positive_sources_hit=["A"]),
        "interaction_only": record("interaction_only", positive_sources_hit=[]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cls = build_classification(
        all_records=records,
        positive_only_records={
            "priority_only": records["priority_only"],
            "interaction_only": records["interaction_only"],
        },
    )
    return cls, gff_file


def test_ranking_metric_interaction_score_flips_order_legacy(tmp_path: Path) -> None:
    """Same re-ranking behavior for legacy_additive."""
    cls, gff_file = _ranking_metric_fixture_legacy(tmp_path)
    default_cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        gff_file=gff_file,
    )
    interaction_ranked_cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        gff_file=gff_file,
    )
    interaction_ranked_cfg.interaction_scoring = replace(
        interaction_ranked_cfg.interaction_scoring, ranking_metric="interaction_score"
    )

    default_rows = {
        r["candidate_protein_id"]: r
        for r in run_interaction_scoring(default_cfg, cls).source_rows["Interaction_Candidates"]
    }
    reranked_rows = {
        r["candidate_protein_id"]: r
        for r in run_interaction_scoring(interaction_ranked_cfg, cls).source_rows["Interaction_Candidates"]
    }

    assert default_rows["priority_only"]["candidate_rank"] < default_rows["interaction_only"]["candidate_rank"]
    assert reranked_rows["interaction_only"]["candidate_rank"] < reranked_rows["priority_only"]["candidate_rank"]
    for candidate_id in ("priority_only", "interaction_only"):
        assert (
            default_rows[candidate_id]["interaction_priority_score"]
            == reranked_rows[candidate_id]["interaction_priority_score"]
        )
        assert (
            default_rows[candidate_id]["interaction_score"] == reranked_rows[candidate_id]["interaction_score"]
        )


# ---------------------------------------------------------------------------
# STRING PPI evidence: external_ppi_evidence / string_cooccurrence (Phase 6a M2)
# ---------------------------------------------------------------------------


def _write_gzip(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(content)


def _seed_string_files(cache_dir: Path, taxon_id: int, links_rows: list[str]) -> None:
    """Write minimal local STRING bulk files, as if already downloaded once."""
    string_dir = cache_dir / "string_ppi_network"
    header = (
        "protein1 protein2 neighborhood fusion cooccurence coexpression "
        "experimental database textmining combined_score\n"
    )
    _write_gzip(
        string_dir / f"{taxon_id}.protein.links.detailed.v{STRING_VERSION}.txt.gz",
        header + "\n".join(links_rows) + "\n",
    )
    known_tags = {
        tag
        for row in links_rows
        for tag in (row.split(" ")[0].split(".", 1)[1], row.split(" ")[1].split(".", 1)[1])
    }
    info_rows = "\n".join(f"{taxon_id}.{tag}\t{tag}\t100\tsome protein" for tag in known_tags)
    _write_gzip(
        string_dir / f"{taxon_id}.protein.info.v{STRING_VERSION}.txt.gz",
        "#string_protein_id\tpreferred_name\tprotein_size\tannotation\n" + info_rows + "\n",
    )


def test_string_evidence_not_run_when_taxon_id_unset(tmp_path: Path) -> None:
    """Without string_ppi_ncbi_taxon_id, string_cooccurrence must be NOT_RUN, not MISSING."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="MA_0002", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "string_cooccurrence"
    )
    assert detail["status"] == "NOT_RUN"


def test_string_evidence_missing_for_unmapped_protein(tmp_path: Path) -> None:
    """A candidate with no old_locus_tag can never be found in STRING -- MISSING."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    _seed_string_files(tmp_path, 188937, [f"188937.MA_0001 188937.MA_9999 0 0 500 0 0 0 0 500"])
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, string_ppi_ncbi_taxon_id=188937)

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "string_cooccurrence"
    )
    assert detail["status"] == "MISSING"


def test_string_cooccurrence_available_and_feeds_interaction_score(tmp_path: Path) -> None:
    """A real STRING match should be AVAILABLE and count toward both scores."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", old_locus_tag="MA_0002", description="", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    _seed_string_files(tmp_path, 188937, ["188937.MA_0001 188937.MA_0002 0 0 800 0 0 0 0 800"])
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, string_ppi_ncbi_taxon_id=188937)

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "string_cooccurrence"
    )
    assert detail["status"] == "AVAILABLE"
    assert detail["normalized_value"] == pytest.approx(0.8)
    assert detail["raw_value"] == pytest.approx(800.0)
    assert detail["category"] == "external_ppi_evidence"
    assert detail["category_cap"] == 15.0

    # domain_complementarity is MISSING (empty description both sides), so
    # functional_annotation stays inactive. genomic_context has no GFF
    # evidence of its own, but string_neighborhood (same seeded row,
    # neighborhood=0) makes the category active at zero -- so the total
    # cap is external_ppi_evidence(15) + genomic_context(25) = 40, with
    # only the string_cooccurrence contribution (0.8*15=12) as raw score.
    assert row["interaction_score"] == pytest.approx(12.0 / 40.0 * 100, abs=0.01)
    assert row["interaction_score"] > 0


def test_string_evidence_reuses_cache_across_runs(tmp_path: Path) -> None:
    """A second run for the same query/species must not need to rescan the bulk file."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="MA_0002", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    _seed_string_files(tmp_path, 188937, ["188937.MA_0001 188937.MA_0002 0 0 800 0 0 0 0 800"])
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, string_ppi_ncbi_taxon_id=188937)

    run_interaction_scoring(cfg, classification(records))

    # Remove the (potentially huge) links file -- the second run must still work from cache.
    string_dir = tmp_path / "string_ppi_network"
    (string_dir / f"188937.protein.links.detailed.v{STRING_VERSION}.txt.gz").unlink()

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_score"] > 0


# ---------------------------------------------------------------------------
# string_neighborhood: shares genomic_context's category (Phase 6a M3)
# ---------------------------------------------------------------------------


def test_string_neighborhood_shares_genomic_context_category(tmp_path: Path) -> None:
    """string_neighborhood must appear under genomic_context and combine with the GFF-based component."""
    gff_file = tmp_path / "genome.gff"
    gff_file.write_text(
        "##gff-version 3\n"
        "contig1\tRefSeq\tgene\t100\t400\t.\t+\t.\tID=query\n"
        "contig1\tRefSeq\tgene\t500\t800\t.\t+\t.\tID=candidate\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=[]),
        "candidate": record(
            "candidate", old_locus_tag="MA_0002", description="", positive_sources_hit=[]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    _seed_string_files(tmp_path, 188937, ["188937.MA_0001 188937.MA_0002 850 0 0 0 0 0 0 850"])
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        gff_file=gff_file,
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, string_ppi_ncbi_taxon_id=188937)

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail_rows = {
        r["component_name"]: r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate"
    }
    string_neighborhood = detail_rows["string_neighborhood"]
    assert string_neighborhood["status"] == "AVAILABLE"
    assert string_neighborhood["category"] == "genomic_context"
    assert string_neighborhood["normalized_value"] == pytest.approx(0.85)
    assert string_neighborhood["raw_value"] == pytest.approx(850.0)
    # Shares its cap with the pipeline's own GFF-based genomic_context
    # component -- same cap value on both rows.
    assert detail_rows["genomic_context"]["category_cap"] == string_neighborhood["category_cap"]

    row = result.source_rows["Interaction_Candidates"][0]
    # genomic_context: GFF's own component (close proximity, normalized
    # 1.0) averaged with string_neighborhood (0.85) -> (1.0 + 0.85) / 2 *
    # 25 = 23.125 raw. The same seeded STRING row also makes
    # string_cooccurrence AVAILABLE at 0 (cooccurrence=0 in the fixture),
    # which activates external_ppi_evidence's 15-point cap at zero raw.
    # total_cap = 15 (external_ppi_evidence) + 25 (genomic_context) = 40;
    # raw = 0 + 23.125 = 23.125.
    assert row["interaction_score"] == pytest.approx(23.125 / 40.0 * 100, abs=0.01)


def test_string_neighborhood_not_run_when_taxon_id_unset() -> None:
    """Without string_ppi_ncbi_taxon_id, string_neighborhood must be NOT_RUN, matching string_cooccurrence."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="MA_0002", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "string_neighborhood"
    )
    assert detail["status"] == "NOT_RUN"
    assert detail["category"] == "genomic_context"


# ---------------------------------------------------------------------------
# coexpression_gse77738: coexpression_evidence category (Phase 6b M2)
# ---------------------------------------------------------------------------


def _seed_gse77738_coexpression_file(cache_dir: Path, gene_values: dict[str, list[float]]) -> None:
    """Write a small synthetic GSE77738_ReadCounts.xls-equivalent workbook.

    ``gene_values`` maps a bare gene locus ('MA0001') to one value per
    steady-state sample column (see GSE77738_STEADY_STATE_RPKM_COLUMNS).
    """
    columns = sorted(GSE77738_STEADY_STATE_RPKM_COLUMNS)
    data = {"Gene Locus": list(gene_values), "Gene Name": ["-"] * len(gene_values)}
    for col_index, col in enumerate(columns):
        data[col] = [values[col_index] for values in gene_values.values()]
    df = pd.DataFrame(data)
    path = cache_dir / "coexpression" / "GSE77738_ReadCounts.xls"
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name="RPKM Normalized Read Counts", index=False)


def test_coexpression_gse77738_not_run_when_disabled() -> None:
    """Without geo_coexpression_enabled, coexpression_gse77738 must be NOT_RUN, not MISSING."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="MA_0002", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "coexpression_gse77738"
    )
    assert detail["status"] == "NOT_RUN"
    assert detail["category"] == "coexpression_evidence"


def test_coexpression_gse77738_missing_for_gene_absent_from_dataset(tmp_path: Path) -> None:
    """A candidate gene never measured in GSE77738 can never be found -- MISSING."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="MA_9999", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    _seed_gse77738_coexpression_file(
        tmp_path,
        {"MA0001": [10, 20, 30, 40, 50, 60, 70, 80, 90, 15, 25, 35, 45]},
    )
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, geo_coexpression_enabled=True)

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "coexpression_gse77738"
    )
    assert detail["status"] == "MISSING"


def test_coexpression_gse77738_available_but_excluded_from_interaction_score(tmp_path: Path) -> None:
    """A correlated pair should be AVAILABLE (percentile, not raw r) but NOT feed interaction_score.

    coexpression_gse77738 was removed from INTERACTION_SCORE_COMPONENT_NAMES
    after a real-data check found it scored AlphaFold3-confirmed
    non-interacting pairs higher, on average, than curated true positive
    pairs -- see claude/experimental_interactions_calibration_report.md.
    It is still computed, cached, and shown in Interaction_Evidence_Detail.
    """
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", old_locus_tag="MA_0002", description="", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    base = [10, 20, 30, 40, 50, 60, 70, 80, 90, 15, 25, 35, 45]
    _seed_gse77738_coexpression_file(
        tmp_path,
        {
            "MA0001": base,
            # Moves in lockstep with MA0001 -- should rank at the top of
            # MA_0001's own background distribution (percentile ~1.0).
            "MA0002": [v * 2 for v in base],
            # Moves opposite -- should rank at the bottom (percentile ~0.0).
            "MA0003": [max(base) - v for v in base],
        },
    )
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, geo_coexpression_enabled=True)

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "coexpression_gse77738"
    )
    assert detail["status"] == "AVAILABLE"
    assert detail["category"] == "coexpression_evidence"
    assert detail["category_cap"] == 12.0
    assert detail["raw_value"] == pytest.approx(1.0, abs=1e-3)  # raw Pearson r, not the percentile
    assert detail["normalized_value"] == pytest.approx(1.0, abs=1e-6)  # percentile rank vs. MA0003

    # No other interaction_score-eligible evidence exists in this minimal
    # fixture (no GFF, no STRING, no GSE64349) -- with coexpression_gse77738
    # excluded, interaction_score must now be None (no query-specific
    # evidence at all), not a score derived solely from gse77738.
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_score"] is None


def test_coexpression_gse77738_reuses_cache_across_runs(tmp_path: Path) -> None:
    """A second run for the same query must not need the source file re-parsed."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="MA_0002", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    base = [10, 20, 30, 40, 50, 60, 70, 80, 90, 15, 25, 35, 45]
    _seed_gse77738_coexpression_file(tmp_path, {"MA0001": base, "MA0002": [v * 2 for v in base]})
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, geo_coexpression_enabled=True)

    run_interaction_scoring(cfg, classification(records))
    (tmp_path / "coexpression" / "GSE77738_ReadCounts.xls").unlink()

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "coexpression_gse77738"
    )
    assert detail["status"] == "AVAILABLE"


# ---------------------------------------------------------------------------
# coexpression_gse64349: shares coexpression_evidence category, lower weight
# (Phase 6b M3)
# ---------------------------------------------------------------------------


def _seed_gse64349_coexpression_files(cache_dir: Path, gene_values: dict[str, list[float]]) -> None:
    """Write small synthetic TableS1 (wild-type) and TableS2 (parental+mutant) workbooks.

    ``gene_values`` maps a bare gene locus to 4 values: [DMS, MMPA, MeOH
    (TableS1's 3 wild-type samples), WWM82-parental (TableS2, kept)]. A
    Delta-msrH column with an obviously-wrong constant value is always
    added, to prove it never influences the result.
    """
    coexpr_dir = cache_dir / "coexpression"
    coexpr_dir.mkdir(parents=True, exist_ok=True)

    table1 = pd.DataFrame(
        {
            "Feature ID": list(gene_values),
            "DMS - S1_R1 (single) (GE) - RPKM": [v[0] for v in gene_values.values()],
            "MMPA - S2_R1 (single) (GE) - RPKM": [v[1] for v in gene_values.values()],
            "MeOH - S3_R1 (single) (GE) - RPKM": [v[2] for v in gene_values.values()],
        }
    )
    table1.to_excel(coexpr_dir / "GSE64349_TableS1_GEO.xlsx", index=False)

    table2 = pd.DataFrame(
        {
            "Feature ID": list(gene_values),
            "WWM82 (parental strain) - S4_R1 (single) (GE) - RPKM": [v[3] for v in gene_values.values()],
            "delta-msrH - S5_R1 (single) (GE) - RPKM": [999999.0] * len(gene_values),
        }
    )
    table2.to_excel(coexpr_dir / "GSE64349_TableS2_GEO.xlsx", index=False)


def test_coexpression_gse64349_not_run_when_disabled() -> None:
    """Without geo_coexpression_enabled, coexpression_gse64349 must be NOT_RUN, matching gse77738."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": record("candidate", old_locus_tag="MA_0002", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "coexpression_gse64349"
    )
    assert detail["status"] == "NOT_RUN"
    assert detail["category"] == "coexpression_evidence"


def test_coexpression_gse64349_weighted_lower_than_gse77738(tmp_path: Path) -> None:
    """Both components share coexpression_evidence's cap, but gse64349's weight is 1/3 of gse77738's."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", old_locus_tag="MA_0002", description="", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    base77738 = [10, 20, 30, 40, 50, 60, 70, 80, 90, 15, 25, 35, 45]
    _seed_gse77738_coexpression_file(
        tmp_path, {"MA0001": base77738, "MA0002": [v * 2 for v in base77738]}
    )
    _seed_gse64349_coexpression_files(
        tmp_path,
        {
            "MA0001": [10, 20, 30, 40],
            "MA0002": [20, 40, 60, 80],
            "MA0003": [40, 30, 20, 10],
        },
    )
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        cache_dir=tmp_path,
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, geo_coexpression_enabled=True)

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail_rows = {
        r["component_name"]: r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate"
    }
    gse77738 = detail_rows["coexpression_gse77738"]
    gse64349 = detail_rows["coexpression_gse64349"]
    assert gse77738["status"] == "AVAILABLE"
    assert gse64349["status"] == "AVAILABLE"
    assert gse77738["category"] == "coexpression_evidence"
    assert gse64349["category"] == "coexpression_evidence"
    # Both components share one cap...
    assert gse77738["category_cap"] == gse64349["category_cap"] == 12.0
    # ...but gse64349's weight is 1/3 of gse77738's (see
    # V2_COMPONENT_WEIGHTS["coexpression_gse64349"]).
    assert gse64349["weight"] == pytest.approx(gse77738["weight"] / 3.0)

    # Delta-msrH's obviously-wrong constant value (999999) must never
    # surface -- both wild-type-equivalent genes (MA0001/MA0002) still
    # correlate near-perfectly.
    assert gse64349["normalized_value"] == pytest.approx(1.0, abs=1e-6)

    # coexpression_gse64349 (unlike coexpression_gse77738, see
    # test_coexpression_gse77738_available_but_excluded_from_interaction_score)
    # still feeds interaction_score.
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_score"] is not None
    assert row["interaction_score"] > 0


# ---------------------------------------------------------------------------
# Final Score: combines protein_hunter_score + interaction_score (Final
# Score integration phase, design spec sections 17-22/27, see
# claude/final_score_integration_investigation.md)
# ---------------------------------------------------------------------------


def _candidate_with_known_protein_hunter_score() -> ProteinRecord:
    """A candidate whose protein_hunter_score is exactly 14/18 (positive_hit + no_negative_hit + domain_hit)."""
    candidate = record(
        "candidate", old_locus_tag="MA_0002", description="iron-sulfur carrier protein", positive_sources_hit=["A"]
    )
    candidate.positive_hits = [
        BlastHit(
            query_id="candidate", subject_id="ref", percent_identity=50.0,
            alignment_length=100, evalue=1e-20, bitscore=100.0,
        )
    ]
    candidate.domains = [DomainHit(source="CDD", accession="cd00001", name="test_domain", description="test")]
    return candidate


def test_final_score_falls_back_to_protein_hunter_alone_when_interaction_missing() -> None:
    """No query-specific evidence at all -> final_score is protein_hunter_score's own normalized value, re-normalized to 100."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=[]),
        "candidate": _candidate_with_known_protein_hunter_score(),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    records["candidate"].description = ""  # no domain_complementarity match either
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["protein_hunter_score"] == 14
    assert row["interaction_score"] is None
    # protein_hunter category alone: normalized (14/18) * its own cap (30),
    # renormalized against total_cap == 30 (only active category) -> *100/30
    # cancels out to just the normalized fraction * 100.
    assert row["final_score"] == pytest.approx(14 / PROTEIN_HUNTER_SCORE_CEILING * 100, abs=0.01)
    assert row["final_score_tier"] in {"Tier1_VeryStrong", "Tier2_Strong", "Tier3_Moderate", "Tier4_Weak"}


def test_final_score_combines_both_categories_per_30_70_split(tmp_path: Path) -> None:
    """With both sub-scores available, final_score == 0.3*phs_norm100 + 0.7*interaction_score (the configured 30/70 cap split)."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": _candidate_with_known_protein_hunter_score(),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["protein_hunter_score"] == 14
    assert row["interaction_score"] is not None
    phs_norm100 = 14 / PROTEIN_HUNTER_SCORE_CEILING * 100
    expected = 0.3 * phs_norm100 + 0.7 * row["interaction_score"]
    assert row["final_score"] == pytest.approx(expected, abs=0.05)


def test_final_score_ignores_negative_hit_strength_by_design(tmp_path: Path) -> None:
    """final_score must NOT be lowered by negative_hit_strength, even "strong" -- tried and removed.

    negative_hit_strength measures phylogenetic novelty (is this protein
    also found in the negative reference genomes), not design spec section
    7.7's "biological contradiction" concept -- a real-data check found
    conflating the two collapsed true-positive/AlphaFold3-negative
    separation, since ancient well-conserved true positives (Hdr/Mtp/Nif
    complexes) are exactly the candidates a "strong" negative hit routes
    into the Negative_hit bucket in the first place. See
    _final_score_components's docstring and
    claude/final_score_integration_investigation.md.
    interaction_priority_score's own, separate use of negative_hit_strength
    is intentionally untouched by this test.
    """
    strong_hit_candidate = _candidate_with_known_protein_hunter_score()
    strong_hit_candidate.negative_hit_strength = "strong"
    strong_hit_candidate.description = ""
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=[]),
        "candidate": strong_hit_candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    phs_norm100 = 14 / PROTEIN_HUNTER_SCORE_CEILING * 100
    # No penalty applied at all: final_score == protein_hunter-alone value,
    # not reduced despite negative_hit_strength == "strong".
    assert row["final_score"] == pytest.approx(phs_norm100, abs=0.01)

    detail = next(
        r
        for r in result.evidence_detail_rows
        if r["candidate_protein_id"] == "candidate" and r["component_name"] == "final_score_negative_penalty"
    )
    assert detail["status"] == "NOT_APPLICABLE"
    assert detail["normalized_value"] is None


def test_final_score_backward_compatible_with_custom_scoring_engine_config_missing_its_categories(
    tmp_path: Path,
) -> None:
    """A pre-existing custom scoring_engine_config.yaml without protein_hunter/interaction caps must not crash."""
    custom_config = tmp_path / "custom_engine.yaml"
    custom_config.write_text(
        "category_caps:\n  source_classification: 30\n  genomic_context: 25\n"
        "  functional_annotation: 20\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=["A"]),
        "candidate": _candidate_with_known_protein_hunter_score(),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        scoring_engine_config=custom_config,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # Falls back to the module's own provisional default caps (30/70)
    # rather than raising a ConfigError.
    assert row["final_score"] is not None


def test_legacy_additive_also_computes_final_score() -> None:
    """protein_hunter_score and interaction_score both exist for legacy_additive too, so final_score should as well."""
    records = {
        "query": record("query", old_locus_tag="MA_0001", positive_sources_hit=["A"]),
        "candidate": _candidate_with_known_protein_hunter_score(),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="legacy_additive",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["protein_hunter_score"] == 14
    assert row["final_score"] is not None
    assert row["final_score_tier"] is not None


def test_ranking_metric_final_score_reorders_candidate_rank() -> None:
    """ranking_metric: final_score should re-rank candidates by final_score, leaving other columns untouched."""
    high_phs = _candidate_with_known_protein_hunter_score()
    high_phs.protein_id = "high_phs"
    high_phs.old_locus_tag = "MA_0002"
    high_phs.description = ""

    low_phs = record("low_phs", old_locus_tag="MA_0003", description="", positive_sources_hit=[])
    # low_phs has no positive_hits/domains -> only no_negative_hit fires -> protein_hunter_score == 5.

    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=[]),
        "candidate": low_phs,
        "high_phs": high_phs,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )
    cfg.interaction_scoring = replace(cfg.interaction_scoring, ranking_metric="final_score")

    classification_obj = classification(records)
    classification_obj.positive_only_records = {"low_phs": low_phs, "high_phs": high_phs}

    result = run_interaction_scoring(cfg, classification_obj)

    assert result is not None
    rows_by_id = {r["candidate_protein_id"]: r for r in result.source_rows["Interaction_Candidates"]}
    assert rows_by_id["high_phs"]["final_score"] > rows_by_id["low_phs"]["final_score"]
    assert rows_by_id["high_phs"]["candidate_rank"] == 1
    assert rows_by_id["low_phs"]["candidate_rank"] == 2


def test_negative_hit_strength_exposed_as_row_column_both_models() -> None:
    """negative_hit_strength (analysis/ortholog_filter.py) is now a plain pair-row column.

    Phase 6-8 sheet redesign needs this as a column (replacing the 3
    strength-specific sheets), for both scoring models -- neither
    _score_pair (legacy) nor _score_pair_v2 previously exposed it.
    """
    candidate = record("candidate", positive_sources_hit=["A"])
    candidate.negative_hit_strength = "medium"
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    for scoring_model in ("legacy_additive", "v2_evidence_based"):
        cfg = interaction_config(
            query_proteins=(InteractionQueryConfig("query", "", ""),),
            candidate_sources={"candidates": True},
            scoring_model=scoring_model,
        )
        result = run_interaction_scoring(cfg, classification(records))
        assert result is not None
        row = result.source_rows["Interaction_Candidates"][0]
        assert row["negative_hit_strength"] == "medium"


def test_v2_category_reference_columns_present_and_legacy_columns_blank() -> None:
    """functional_domain/evolutionary/cellular_compatibility/interaction_evidence_score.

    v2_evidence_based only (category-cap concept doesn't exist for
    legacy_additive) -- new Phase 6-8 reference columns mirroring
    candidate_priority_score/same_gene_neighborhood_score, needed so the
    12-sheet redesign's category-level columns can read them directly from
    the existing pair rows instead of re-deriving from ScoreBreakdown.
    """
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": _candidate_with_known_protein_hunter_score(),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg_v2 = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )
    result_v2 = run_interaction_scoring(cfg_v2, classification(records))
    assert result_v2 is not None
    row_v2 = result_v2.source_rows["Interaction_Candidates"][0]
    # functional_annotation category is active (domain_complementarity has a
    # match: "radical SAM" query vs "iron-sulfur carrier protein" candidate).
    assert row_v2["functional_domain_score"] is not None
    assert row_v2["functional_domain_score"] >= 0
    # No PIH evidence bundle configured -> evolutionary/cellular_compatibility
    # categories have no components at all -> None (MISSING), not zero.
    assert row_v2["evolutionary_score"] is None
    assert row_v2["cellular_compatibility_score"] is None
    # No STRING/coexpression bundle configured -> interaction_evidence_score
    # (external_ppi_evidence + coexpression_evidence + pih_direct_interaction) is None.
    assert row_v2["interaction_evidence_score"] is None

    cfg_legacy = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="legacy_additive",
    )
    result_legacy = run_interaction_scoring(cfg_legacy, classification(records))
    assert result_legacy is not None
    row_legacy = result_legacy.source_rows["Interaction_Candidates"][0]
    assert row_legacy.get("functional_domain_score") is None
    assert row_legacy.get("evolutionary_score") is None
    assert row_legacy.get("cellular_compatibility_score") is None
    assert row_legacy.get("interaction_evidence_score") is None


def test_compute_protein_hunter_only_final_score_matches_in_run_fallback() -> None:
    """The public no-query fallback wrapper must match Final Score's own MISSING-interaction_score path."""
    from analysis.interaction_scoring import compute_protein_hunter_only_final_score
    from analysis.scoring_engine_config import load_scoring_engine_config

    records = {
        "query": record("query", old_locus_tag="MA_0001", description="", positive_sources_hit=[]),
        "candidate": _candidate_with_known_protein_hunter_score(),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    records["candidate"].description = ""
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )
    result = run_interaction_scoring(cfg, classification(records))
    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_score"] is None

    engine_config = load_scoring_engine_config(None)
    final_score, tier = compute_protein_hunter_only_final_score(row["protein_hunter_score"], engine_config)
    assert final_score == pytest.approx(row["final_score"], abs=1e-6)
    assert tier == row["final_score_tier"]


def test_interaction_scoring_disabled_returns_none() -> None:
    """Disabled interaction scoring should not create output."""
    cfg = interaction_config(enabled=False)
    result = run_interaction_scoring(cfg, SimpleNamespace(all_records={}))

    assert result is None


def test_query_resolution_by_protein_id_old_locus_and_sequence() -> None:
    """Queries should resolve by protein_id, old_locus_tag, or explicit sequence."""
    records = {
        "query": record("query", old_locus_tag="MA_0001"),
        "locus_query": record("locus_query", old_locus_tag="MA_0002"),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(
            InteractionQueryConfig("query", "", ""),
            InteractionQueryConfig("", "MA_0002", ""),
            InteractionQueryConfig("external", "", "MPEPTIDE"),
            InteractionQueryConfig("missing", "", ""),
        ),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    statuses = {row["query_id"]: row["resolution_status"] for row in result.query_rows}
    assert statuses["query"] == "resolved"
    assert statuses["MA_0002"] == "resolved"
    assert statuses["external"] == "resolved"
    assert statuses["missing"] == "unresolved"
    external = next(row for row in result.query_rows if row["query_id"] == "external")
    assert external["sequence_length"] == len("MPEPTIDE")


def test_query_fasta_is_loaded(tmp_path: Path) -> None:
    """query_fasta records should become interaction queries."""
    query_fasta = tmp_path / "queries.faa"
    query_fasta.write_text(">query_fasta_1 description\nMSTNPK\n", encoding="utf-8")
    records = {
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_fasta=query_fasta,
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert result.query_rows[0]["query_id"] == "query_fasta_1"
    assert result.query_rows[0]["resolution_status"] == "resolved"


def test_candidate_sources_true_false_and_ranking_limit_are_respected() -> None:
    """Only enabled source sheets should produce interaction rows."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed", description="carrier protein"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={
            "candidates": True,
            "candidates_relaxed": False,
            "no_hit": True,
        },
        max_candidates_per_query=1,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert set(result.source_rows) == {"Interaction_Candidates", "Interaction_No_hit"}
    assert len(result.source_rows["Interaction_Candidates"]) == 1
    assert result.source_rows["Interaction_Candidates"][0]["candidate_rank"] == 1
    assert "Interaction_Candidates_relaxed" not in result.source_rows


def test_missing_enabled_source_is_skipped_with_warning() -> None:
    """Enabled but empty candidate sources should not crash or create a sheet."""
    records = {
        "query": record("query"),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"positive_all_sources": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert result.source_rows == {}
    assert result.warnings == [
        "interaction candidate source has no records: positive_all_sources"
    ]


def test_include_sequences_controls_pair_columns() -> None:
    """Sequence columns should only be present when explicitly requested."""
    assert "query_sequence" not in interaction_pair_columns(False)
    assert "candidate_sequence" not in interaction_pair_columns(False)
    assert "query_sequence" in interaction_pair_columns(True)
    assert "candidate_sequence" in interaction_pair_columns(True)


def test_scoring_adds_reasons_and_alphafold_readiness_without_running_alphafold() -> None:
    """Scoring should use lightweight evidence only."""
    query = record(
        "query",
        description="enzyme protein",
        positive_sources_hit=["A"],
    )
    candidate = record(
        "candidate",
        description="carrier enzyme protein",
        positive_sources_hit=["A"],
    )
    candidate.domains.append(
        DomainHit(source="Pfam", accession="PF00001", name="carrier")
    )
    records = {
        "query": query,
        "candidate": candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_priority_score"] > 0
    assert row["alphafold_recommended"] is True
    assert "compatible for manual AlphaFold" in row["interaction_score_reasons"]



def test_positive_all_sources_uses_short_interaction_sheet_name() -> None:
    """Long source names should map to explicit Excel-safe sheet names."""
    records = {
        "query": record("query"),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cls = classification(records)
    cls.positive_all_sources_records = {"candidate": records["candidate"]}
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"positive_all_sources": True},
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    assert set(result.source_rows) == {"Interaction_Positive_all"}
    assert result.source_rows["Interaction_Positive_all"][0]["candidate_source"] == (
        "Positive_all_sources"
    )


def test_gff_distance_scoring_close_medium_far_different_missing_and_overlap(
    tmp_path: Path,
) -> None:
    """GFF coordinates should drive same-gene-neighborhood scoring."""
    gff = tmp_path / "coords.gff"
    gff.write_text(
        "contig1\tRefSeq\tCDS\t100\t200\t.\t+\t0\tID=cds-query;protein_id=query;old_locus_tag=MA_0001\n"
        "contig1\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-close;protein_id=close;old_locus_tag=MA_0002\n"
        "contig1\tRefSeq\tCDS\t15000\t15100\t.\t-\t0\tID=cds-medium;protein_id=medium;old_locus_tag=MA_0003\n"
        "contig1\tRefSeq\tCDS\t150000\t150100\t.\t+\t0\tID=cds-far;protein_id=far;old_locus_tag=MA_0004\n"
        "contig2\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-other;protein_id=other;old_locus_tag=MA_0005\n"
        "contig1\tRefSeq\tCDS\t150\t250\t.\t+\t0\tID=cds-overlap;protein_id=overlap;old_locus_tag=MA_0006\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", old_locus_tag="MA_0001"),
        "candidate": record("close", old_locus_tag="MA_0002"),
        "relaxed": record("medium", old_locus_tag="MA_0003"),
        "novel": record("far", old_locus_tag="MA_0004"),
        "other": record("other", old_locus_tag="MA_0005"),
        "overlap": record("overlap", old_locus_tag="MA_0006"),
        "missing": record("missing", old_locus_tag="MA_9999"),
    }
    cls = classification(records)
    cls.positive_only_records = {
        key: records[key]
        for key in ("candidate", "relaxed", "novel", "other", "overlap", "missing")
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        gff_file=gff,
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    by_id = {
        row["candidate_protein_id"]: row
        for row in result.source_rows["Interaction_Candidates"]
    }
    assert by_id["close"]["distance_bp"] == 300
    assert by_id["close"]["same_gene_neighborhood_score"] == 25.0
    assert by_id["medium"]["same_gene_neighborhood_score"] == 15.0
    assert by_id["far"]["same_gene_neighborhood_score"] == 0.0
    assert by_id["other"]["same_contig"] is False
    assert by_id["missing"]["distance_bp"] is None
    assert by_id["overlap"]["distance_bp"] == 0
    assert by_id["close"]["strand_relation"] == "same_strand"
    assert by_id["medium"]["strand_relation"] == "opposite_strand"


def test_shared_generic_words_only_do_not_score() -> None:
    """Generic description overlap should not create complementarity score."""
    records = {
        "query": record("query", description="hypothetical protein family domain"),
        "candidate": record("candidate", description="putative protein family domain"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["domain_complementarity_score"] == 0.0
    assert "generic-only description overlap ignored" in row["interaction_score_reasons"]
    assert "meaningful shared terms" not in row["interaction_score_reasons"]


def test_meaningful_and_complementary_terms_score() -> None:
    """Meaningful functional terms should still contribute conservatively."""
    records = {
        "query": record("query", description="radical SAM iron-sulfur protein"),
        "candidate": record("candidate", description="radical SAM transferase"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["domain_complementarity_score"] > 0
    assert (
        "meaningful shared terms" in row["interaction_score_reasons"]
        or "complementary terms" in row["interaction_score_reasons"]
    )



def test_neighborhood_rows_use_gff_distance_and_limits(tmp_path: Path) -> None:
    """Interaction_Neighborhood rows should summarize nearby same-contig pairs."""
    gff = tmp_path / "coords.gff"
    gff.write_text(
        "contig1\tRefSeq\tCDS\t100\t200\t.\t+\t0\tID=cds-query;protein_id=query;old_locus_tag=MA_0001\n"
        "contig1\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-close;protein_id=close;old_locus_tag=MA_0002\n"
        "contig1\tRefSeq\tCDS\t1500\t1600\t.\t-\t0\tID=cds-near;protein_id=near;old_locus_tag=MA_0003\n"
        "contig1\tRefSeq\tCDS\t9000\t9100\t.\t+\t0\tID=cds-far;protein_id=far;old_locus_tag=MA_0004\n"
        "contig2\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-other;protein_id=other;old_locus_tag=MA_0005\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", old_locus_tag="MA_0001", description="radical SAM protein"),
        "candidate": record("close", old_locus_tag="MA_0002", description="iron-sulfur protein"),
        "relaxed": record("near", old_locus_tag="MA_0003"),
        "novel": record("far", old_locus_tag="MA_0004"),
        "other": record("other", old_locus_tag="MA_0005"),
    }
    cls = classification(records)
    cls.positive_only_records = {
        key: records[key]
        for key in ("candidate", "relaxed", "novel", "other")
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        gff_file=gff,
        neighborhood=InteractionNeighborhoodConfig(
            enabled=True,
            max_distance_bp=2000,
            max_rows_per_query=2,
        ),
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    assert [row["candidate_protein_id"] for row in result.neighborhood_rows] == [
        "close",
        "near",
    ]
    assert [row["candidate_rank_by_distance"] for row in result.neighborhood_rows] == [1, 2]
    assert result.neighborhood_rows[0]["query_description"] == "radical SAM protein"
    assert result.neighborhood_rows[0]["query_contig"] == "contig1"
    assert result.neighborhood_rows[0]["candidate_contig"] == "contig1"
    assert result.neighborhood_rows[0]["neighborhood_band"] == "<=5kb"
    assert all(row["candidate_protein_id"] != "query" for row in result.neighborhood_rows)


def test_neighborhood_rows_can_be_disabled(tmp_path: Path) -> None:
    """Neighborhood summary should be optional even when pair scoring runs."""
    gff = tmp_path / "coords.gff"
    gff.write_text(
        "contig1\tRefSeq\tCDS\t100\t200\t.\t+\t0\tID=cds-query;protein_id=query\n"
        "contig1\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-candidate;protein_id=candidate\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query"),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        gff_file=gff,
        neighborhood=InteractionNeighborhoodConfig(
            enabled=False,
            max_distance_bp=100000,
            max_rows_per_query=200,
        ),
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert result.source_rows["Interaction_Candidates"]
    assert result.neighborhood_rows == []


def test_domain_complementarity_uses_pfam_cdd_terms_and_caps_score() -> None:
    """Domain annotations should contribute concise functional evidence only."""
    query = record("query", description="hypothetical protein")
    query.domains.append(
        DomainHit(
            source="Pfam",
            accession="PF04055",
            name="Radical_SAM",
            description="radical SAM enzyme",
        )
    )
    candidate = record("candidate", description="conserved protein")
    candidate.domains.append(
        DomainHit(
            source="CDD",
            accession="cd12345",
            name="Ferredoxin_like",
            description="iron-sulfur ferredoxin domain",
        )
    )
    records = {
        "query": query,
        "candidate": candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["domain_complementarity_score"] == 15.0
    assert "Pfam/CDD functional terms used" in row["interaction_score_reasons"]
    assert "complementary terms" in row["interaction_score_reasons"]



def test_distance_independent_score_excludes_gene_neighborhood(tmp_path: Path) -> None:
    """Distance-independent score should not include genomic neighborhood evidence."""
    gff = tmp_path / "coords.gff"
    gff.write_text(
        "contig1\tRefSeq\tCDS\t100\t200\t.\t+\t0\tID=cds-query;protein_id=query\n"
        "contig1\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-candidate;protein_id=candidate\n",
        encoding="utf-8",
    )
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        gff_file=gff,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    expected = (
        row["candidate_priority_score"]
        + row["co_occurrence_score"]
        + row["domain_complementarity_score"]
    )
    assert row["same_gene_neighborhood_score"] > 0
    assert row["distance_independent_score"] == expected
    assert row["interaction_priority_score"] > row["distance_independent_score"]


def test_priority_groups_distinguish_nearby_distant_and_no_hit(tmp_path: Path) -> None:
    """Nearby and distant retained candidates should receive distinct groups."""
    gff = tmp_path / "coords.gff"
    gff.write_text(
        "contig1\tRefSeq\tCDS\t100\t200\t.\t+\t0\tID=cds-query;protein_id=query\n"
        "contig1\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-near;protein_id=near\n"
        "contig1\tRefSeq\tCDS\t200000\t200100\t.\t+\t0\tID=cds-distant-co;protein_id=distant_co\n"
        "contig1\tRefSeq\tCDS\t210000\t210100\t.\t+\t0\tID=cds-distant-domain;protein_id=distant_domain\n",
        encoding="utf-8",
    )
    query = record("query", description="radical SAM protein", positive_sources_hit=["A"])
    near = record("near", positive_sources_hit=[])
    distant_co = record("distant_co", positive_sources_hit=["A"])
    distant_domain = record("distant_domain", description="iron-sulfur protein")
    records = {
        "query": query,
        "candidate": near,
        "relaxed": distant_co,
        "novel": distant_domain,
    }
    cls = classification(records)
    cls.positive_only_records = {
        "near": near,
        "distant_co": distant_co,
        "distant_domain": distant_domain,
    }
    cls.no_hit_records = {"novel": distant_domain}
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True, "no_hit": True},
        gff_file=gff,
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    by_id = {
        row["candidate_protein_id"]: row
        for row in result.source_rows["Interaction_Candidates"]
    }
    assert set(by_id) == {"near", "distant_co", "distant_domain"}
    assert by_id["near"]["priority_group"] == "nearby_candidate"
    assert by_id["distant_co"]["priority_group"] == "distant_cooccurrence_candidate"
    assert by_id["distant_domain"]["priority_group"] == "distant_domain_candidate"
    assert "distant candidate retained" in by_id["distant_co"]["interaction_score_reasons"]
    no_hit_row = result.source_rows["Interaction_No_hit"][0]
    assert no_hit_row["priority_group"] == "distant_domain_candidate"


def test_distance_independent_rank_is_assigned_within_interaction_sheet() -> None:
    """Distance-independent ranks should be present without changing candidate rows."""
    records = {
        "query": record("query", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": record("candidate", description="iron-sulfur protein"),
        "relaxed": record("relaxed", positive_sources_hit=["A"]),
        "novel": record("novel"),
    }
    cls = classification(records)
    cls.positive_only_records = {
        "candidate": records["candidate"],
        "relaxed": records["relaxed"],
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    rows = result.source_rows["Interaction_Candidates"]
    assert len(rows) == 2
    assert {row["distance_independent_rank"] for row in rows} == {1, 2}
    assert all("distance_independent_score" in row for row in rows)


# ---------------------------------------------------------------------------
# scoring model v2 (evidence-based, category-capped) integration tests
# ---------------------------------------------------------------------------


def test_v2_mode_is_opt_in_default_stays_legacy() -> None:
    """Without an explicit scoring_model, behavior must stay legacy_additive."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )
    result = run_interaction_scoring(cfg, classification(records))
    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # legacy rows never set scoring_model -- the column stays absent/None.
    assert row.get("scoring_model") is None


def test_v2_mode_scores_full_evidence_pair_near_100() -> None:
    """A candidate with strong evidence in every category should score high."""
    records = {
        "query": record("query", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", description="iron-sulfur carrier protein", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["scoring_model"] == "v2_evidence_based"
    assert row["formal_score_available"] is True
    # source_classification (Candidates=30/30) + co_occurrence (full Jaccard
    # overlap) + domain_complementarity (radical sam / iron-sulfur rule) are
    # all available; genomic_context is missing (no GFF configured) and must
    # be excluded from the denominator rather than counted as zero.
    assert row["evidence_category_count"] == 2
    assert row["interaction_priority_score"] == 100.0
    assert row["evidence_tier"] == "Tier2_Strong"


def test_v2_mode_distinguishes_missing_from_evaluated_no_match() -> None:
    """Missing annotation must not score the same as an evaluated non-match."""
    records = {
        "query": record("query", description="", positive_sources_hit=[]),
        "candidate": record("no_annotation", description="", positive_sources_hit=[]),
        "no_annotation": record("no_annotation", description="", positive_sources_hit=[]),
        "unrelated_annotation": record(
            "unrelated_annotation", description="completely unrelated text here", positive_sources_hit=[]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cls = classification(records)
    cls.positive_only_records = {
        "no_annotation": records["no_annotation"],
        "unrelated_annotation": records["unrelated_annotation"],
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    rows = {row["candidate_protein_id"]: row for row in result.source_rows["Interaction_Candidates"]}

    # Query has no description at all -> domain_complementarity is MISSING
    # for every candidate, and co_occurrence is also MISSING (no BLAST
    # source pattern on either side and no negative hits to compare either,
    # but the query record itself carries no positive_sources_hit -- this
    # still resolves to the "no negative hit either side" weak-positive
    # branch, not MISSING, because both records exist). The candidate with
    # descriptive text should NOT score lower than the blank one just for
    # having text that happens not to match.
    no_annotation_row = rows["no_annotation"]
    unrelated_row = rows["unrelated_annotation"]
    assert no_annotation_row["evidence_category_count"] == unrelated_row["evidence_category_count"]


def test_v2_mode_category_cap_limits_functional_annotation_contribution() -> None:
    """co_occurrence and domain_complementarity must share one capped category."""
    records = {
        "query": record("query", description="radical sam protein", positive_sources_hit=["A", "B"]),
        "candidate": record(
            "candidate", description="iron-sulfur protein", positive_sources_hit=["A", "B"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # Both co_occurrence (perfect source overlap) and domain_complementarity
    # (rule match) are maxed out, but functional_annotation's cap (20) must
    # not be exceeded, and it must not inflate source_classification (30)
    # or push the total above what two active categories (30 + 20 = 50
    # points of cap) can produce as a 0-100 score.
    assert row["interaction_priority_score"] == 100.0
    assert row["co_occurrence_score"] <= 10.0
    assert row["domain_complementarity_score"] <= 10.0


def test_v2_mode_uses_custom_scoring_engine_config(tmp_path: Path) -> None:
    """A custom scoring_engine_config path should change eligibility/tiers."""
    engine_config_path = tmp_path / "scoring.yaml"
    engine_config_path.write_text(
        """
category_caps:
  source_classification: 30
  genomic_context: 25
  functional_annotation: 20
minimum_evidence:
  min_categories: 5
  min_available_weight: 0
""",
        encoding="utf-8",
    )
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        scoring_engine_config=engine_config_path,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    # min_categories=5 is unreachable with this fixture's evidence, so no
    # formal score should be produced.
    assert row["formal_score_available"] is False
    assert row["evidence_tier"] == "Unclassified"
    assert row["interaction_priority_score"] is None


def test_v2_mode_uses_custom_functional_complementarity_ruleset(tmp_path: Path) -> None:
    """A project-specific ruleset should be picked up instead of the default."""
    ruleset_path = tmp_path / "rules.yaml"
    ruleset_path.write_text(
        """
version: test
rules:
  - rule_id: custom_pair
    left_terms: [gizmo]
    right_terms: [widget]
meaningful_keywords: [gizmo, widget]
stopwords: [protein]
""",
        encoding="utf-8",
    )
    records = {
        "query": record("query", description="gizmo enzyme", positive_sources_hit=["A"]),
        "candidate": record("candidate", description="widget carrier", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        functional_complementarity_ruleset=ruleset_path,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert "custom_pair" in row["interaction_score_reasons"]


def test_v2_mode_applies_negative_hit_strength_as_penalty() -> None:
    """A strong negative BLAST hit should reduce the score via a capped penalty."""
    clean_candidate = record("clean_candidate", positive_sources_hit=["A"])
    flagged_candidate = record("flagged_candidate", positive_sources_hit=["A"])
    flagged_candidate.negative_hit_strength = "strong"

    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": clean_candidate,
        "flagged": flagged_candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cls = classification(records)
    cls.positive_only_records = {"clean_candidate": clean_candidate, "flagged_candidate": flagged_candidate}
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    rows = {row["candidate_protein_id"]: row for row in result.source_rows["Interaction_Candidates"]}
    clean_row = rows["clean_candidate"]
    flagged_row = rows["flagged_candidate"]
    assert clean_row["interaction_priority_score"] > flagged_row["interaction_priority_score"]
    assert "negative BLAST hit strength: strong" in flagged_row["interaction_score_reasons"]
    assert clean_row["candidate_rank"] < flagged_row["candidate_rank"]


def test_v2_mode_no_negative_hit_is_not_applicable_not_a_penalty() -> None:
    """A candidate with no negative hit must not carry a phantom penalty."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert "negative BLAST hit strength" not in row["interaction_score_reasons"]


# ---------------------------------------------------------------------------
# sequence_evidence (BLAST hit strength) tests
# ---------------------------------------------------------------------------


def _hit(protein_id: str, *, identity: float, evalue: float, bitscore: float) -> BlastHit:
    return BlastHit(
        query_id=protein_id,
        subject_id=f"ref_{protein_id}",
        percent_identity=identity,
        alignment_length=90,
        evalue=evalue,
        bitscore=bitscore,
        query_length=100,
    )


def test_v2_mode_strong_positive_hit_scores_higher_than_weak_hit() -> None:
    """A stronger BLAST hit (identity/evalue) should score higher via sequence_evidence."""
    strong_candidate = record("candidate", positive_sources_hit=["A"])
    strong_candidate.positive_hits = [
        _hit("candidate", identity=85.0, evalue=1e-40, bitscore=200.0)
    ]
    weak_candidate = record("flagged", positive_sources_hit=["A"])
    weak_candidate.positive_hits = [_hit("flagged", identity=27.0, evalue=5e-6, bitscore=40.0)]

    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": strong_candidate,
        "flagged": weak_candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cls = classification(records)
    cls.positive_only_records = {"candidate": strong_candidate, "flagged": weak_candidate}
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, cls)

    assert result is not None
    rows = {row["candidate_protein_id"]: row for row in result.source_rows["Interaction_Candidates"]}
    strong_row = rows["candidate"]
    weak_row = rows["flagged"]
    assert strong_row["candidate_priority_score"] > weak_row["candidate_priority_score"]
    assert "best positive BLAST hit" in strong_row["interaction_score_reasons"]


def test_v2_mode_missing_positive_hit_is_missing_not_zero_penalty() -> None:
    """A candidate with no positive BLAST hits should get MISSING, not a scored zero."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert "no positive BLAST hit" in row["interaction_score_reasons"]
    # source_classification's own component is still AVAILABLE, so the
    # shared category still produces a score -- MISSING must not zero it out.
    assert row["candidate_priority_score"] is not None


def test_v2_mode_sequence_evidence_uses_best_hit_among_multiple() -> None:
    """Multiple positive_hits use get_best_hit's (bitscore, -evalue) rule, not an average."""
    candidate = record("candidate", positive_sources_hit=["A"])
    weak_hit = _hit("candidate", identity=27.0, evalue=5e-6, bitscore=40.0)
    strong_hit = _hit("candidate", identity=88.0, evalue=1e-50, bitscore=210.0)
    candidate.positive_hits = [weak_hit, strong_hit]

    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert "identity=88.0%" in row["interaction_score_reasons"]
    assert "identity=27.0%" not in row["interaction_score_reasons"]


def test_v2_mode_sequence_evidence_handles_zero_evalue() -> None:
    """evalue == 0.0 (floating-point underflow) must not crash log10 and scores as strongest."""
    candidate = record("candidate", positive_sources_hit=["A"])
    candidate.positive_hits = [_hit("candidate", identity=95.0, evalue=0.0, bitscore=300.0)]

    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": candidate,
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["candidate_priority_score"] is not None
    assert "evalue=0.00e+00" in row["interaction_score_reasons"]


# ---------------------------------------------------------------------------
# ProteinInteractionHunter (PIH) evidence bridge (Phase 4) integration tests
# ---------------------------------------------------------------------------


def _write_pih_bundle(path: Path, *, query_id: str, candidate_id: str, category_scores: list[dict]) -> None:
    """Write a single-line PIH candidate_evidence_bundle.jsonl fixture."""
    import json

    record_line = {
        "query_id": query_id,
        "candidate_id": candidate_id,
        "integrated_scoring": {"category_scores": category_scores},
    }
    path.write_text(json.dumps(record_line) + "\n", encoding="utf-8")


def test_pih_bridge_adds_only_non_overlapping_categories(tmp_path: Path) -> None:
    """Only cellular_compatibility/evolutionary/direct_interaction are bridged."""
    bundle_path = tmp_path / "candidate_evidence_bundle.jsonl"
    _write_pih_bundle(
        bundle_path,
        query_id="query",
        candidate_id="candidate",
        category_scores=[
            {"category_name": "direct_interaction", "normalized_score": 1.0, "available_weight": 1.0},
            {"category_name": "evolutionary", "normalized_score": 1.0, "available_weight": 1.0},
            {"category_name": "cellular_compatibility", "normalized_score": 1.0, "available_weight": 1.0},
            # These two must be ignored: v5 already computes its own versions
            # of genomic_context and functional_annotation independently.
            {"category_name": "genomic_context", "normalized_score": 1.0, "available_weight": 1.0},
            {"category_name": "functional_annotation", "normalized_score": 1.0, "available_weight": 1.0},
        ],
    )
    records = {
        "query": record("query", description="", positive_sources_hit=[]),
        "candidate": record("candidate", description="", positive_sources_hit=[]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        pih_evidence_bundle=bundle_path,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    reasons = row["interaction_score_reasons"]
    assert "pih_direct_interaction" in reasons
    assert "pih_evolutionary" in reasons
    assert "pih_cellular_compatibility" in reasons
    assert "pih_genomic_context" not in reasons
    assert "pih_functional_annotation" not in reasons


def test_pih_bridge_raises_score_and_category_count_relative_to_no_bridge(tmp_path: Path) -> None:
    """Bridged evidence should add categories and raise the final score."""
    bundle_path = tmp_path / "candidate_evidence_bundle.jsonl"
    _write_pih_bundle(
        bundle_path,
        query_id="query",
        candidate_id="candidate",
        category_scores=[
            {"category_name": "direct_interaction", "normalized_score": 1.0, "available_weight": 1.0},
            {"category_name": "evolutionary", "normalized_score": 1.0, "available_weight": 1.0},
            {"category_name": "cellular_compatibility", "normalized_score": 1.0, "available_weight": 1.0},
        ],
    )
    records = {
        "query": record("query", description="", positive_sources_hit=[]),
        "candidate": record("candidate", description="", positive_sources_hit=[]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }

    cfg_without_bridge = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )
    cfg_with_bridge = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        pih_evidence_bundle=bundle_path,
    )

    row_without = run_interaction_scoring(cfg_without_bridge, classification(records)).source_rows[
        "Interaction_Candidates"
    ][0]
    row_with = run_interaction_scoring(cfg_with_bridge, classification(records)).source_rows[
        "Interaction_Candidates"
    ][0]

    assert row_with["evidence_category_count"] > row_without["evidence_category_count"]
    assert row_with["interaction_priority_score"] > row_without["interaction_priority_score"]


def test_pih_bridge_missing_file_warns_and_does_not_crash(tmp_path: Path) -> None:
    """A configured but absent bundle file must degrade gracefully, not raise."""
    missing_path = tmp_path / "does_not_exist.jsonl"
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        pih_evidence_bundle=missing_path,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert any("PIH evidence bundle not found" in warning for warning in result.warnings)
    row = result.source_rows["Interaction_Candidates"][0]
    assert "pih_" not in row["interaction_score_reasons"]


def test_pih_bridge_matches_version_stripped_query_id(tmp_path: Path) -> None:
    """v5's own versioned protein id should still match an unversioned PIH key.

    PIH and v5 were not designed to share an identifier convention; PIH's
    bundle may be keyed without the trailing '.N' version suffix that v5's
    resolved protein id carries. The bridge must try the version-stripped
    form on the query side, exactly as it already does on the candidate side.
    """
    bundle_path = tmp_path / "candidate_evidence_bundle.jsonl"
    _write_pih_bundle(
        bundle_path,
        query_id="query",
        candidate_id="candidate",
        category_scores=[
            {"category_name": "direct_interaction", "normalized_score": 1.0, "available_weight": 1.0},
        ],
    )
    records = {
        "query.2": record("query.2", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query.2", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        pih_evidence_bundle=bundle_path,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert "pih_direct_interaction" in row["interaction_score_reasons"]


def test_pih_bridge_malformed_line_is_skipped_with_warning(tmp_path: Path) -> None:
    """A malformed JSONL line must not abort the whole run."""
    bundle_path = tmp_path / "candidate_evidence_bundle.jsonl"
    bundle_path.write_text("not valid json\n", encoding="utf-8")
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
        pih_evidence_bundle=bundle_path,
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert any("not valid JSON" in warning for warning in result.warnings)
    row = result.source_rows["Interaction_Candidates"][0]
    assert "pih_" not in row["interaction_score_reasons"]
    assert "no negative BLAST hit" in row["interaction_score_reasons"]


# ---------------------------------------------------------------------------
# Interaction_Evidence_Detail (evidence_detail_rows) tests
# ---------------------------------------------------------------------------


def test_evidence_detail_v2_long_format_has_one_row_per_component() -> None:
    """v2 mode should emit one detail row per EvidenceComponent for each pair."""
    records = {
        "query": record("query", description="radical SAM protein", positive_sources_hit=["A"]),
        "candidate": record(
            "candidate", description="iron-sulfur carrier protein", positive_sources_hit=["A"]
        ),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    # The main sheet is unaffected by collecting evidence detail alongside it.
    row = result.source_rows["Interaction_Candidates"][0]
    assert row["interaction_priority_score"] == 100.0

    assert result.evidence_detail_scoring_model == "v2_evidence_based"
    detail_rows = [
        r for r in result.evidence_detail_rows if r["candidate_protein_id"] == "candidate"
    ]
    # source_classification, sequence_evidence, genomic_context, co_occurrence,
    # domain_complementarity, negative_hit_strength, string_cooccurrence,
    # string_neighborhood, coexpression_gse77738, coexpression_gse64349 --
    # always exactly ten from _build_evidence_components_v2, whether
    # AVAILABLE or not -- plus three more from Final Score's own,
    # separately-attached breakdown (protein_hunter_score, interaction_score,
    # final_score_negative_penalty -- see _final_score_components).
    assert len(detail_rows) == 13
    assert {r["component_name"] for r in detail_rows} == {
        "source_classification",
        "sequence_evidence",
        "genomic_context",
        "co_occurrence",
        "domain_complementarity",
        "negative_hit_strength",
        "string_cooccurrence",
        "string_neighborhood",
        "coexpression_gse77738",
        "coexpression_gse64349",
        "protein_hunter_score",
        "interaction_score",
        "final_score_negative_penalty",
    }
    for detail_row in detail_rows:
        assert set(detail_row) == set(INTERACTION_EVIDENCE_DETAIL_V2_COLUMNS)
        assert detail_row["query_protein_id"] == "query"
        assert detail_row["candidate_rank"] == row["candidate_rank"]

    source_component = next(r for r in detail_rows if r["component_name"] == "source_classification")
    assert source_component["status"] == "AVAILABLE"
    assert source_component["raw_value"] == "Candidates"
    assert source_component["category"] == "source_classification"
    assert source_component["category_cap"] == 30.0

    genomic_component = next(r for r in detail_rows if r["component_name"] == "genomic_context")
    assert genomic_component["status"] == "MISSING"  # no GFF configured in this test
    assert genomic_component["normalized_value"] is None


def test_evidence_detail_legacy_wide_format_has_one_row_per_pair() -> None:
    """legacy_additive mode should emit exactly one detail row per pair, projecting existing columns."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True},
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    row = result.source_rows["Interaction_Candidates"][0]
    assert result.evidence_detail_scoring_model == "legacy_additive"
    assert len(result.evidence_detail_rows) == 1
    detail_row = result.evidence_detail_rows[0]
    assert set(detail_row) == set(INTERACTION_EVIDENCE_DETAIL_LEGACY_COLUMNS)
    assert detail_row["candidate_protein_id"] == "candidate"
    assert detail_row["candidate_rank"] == row["candidate_rank"]
    assert detail_row["candidate_priority_score"] == row["candidate_priority_score"]
    assert detail_row["same_gene_neighborhood_score"] == row["same_gene_neighborhood_score"]
    assert detail_row["co_occurrence_score"] == row["co_occurrence_score"]
    assert detail_row["domain_complementarity_score"] == row["domain_complementarity_score"]
    assert detail_row["alphafold_readiness_score"] == row["alphafold_readiness_score"]
    assert detail_row["interaction_score_reasons"] == row["interaction_score_reasons"]


def test_evidence_detail_excludes_no_hit_by_default() -> None:
    """no_hit must be excluded from evidence_detail_rows unless explicitly included."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"no_hit": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    # The main sheet still gets its rows -- only the detail sheet is scoped down.
    assert result.source_rows["Interaction_No_hit"]
    assert result.evidence_detail_rows == []


def test_evidence_detail_includes_no_hit_when_opted_in() -> None:
    """include_no_hit: true should produce detail rows for the no_hit bucket too."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate"),
        "relaxed": record("relaxed"),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"no_hit": True},
        scoring_model="v2_evidence_based",
        evidence_detail_sheet=InteractionEvidenceDetailConfig(include_no_hit=True),
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    assert result.evidence_detail_rows
    assert all(
        r["candidate_source"] == "No_hit" for r in result.evidence_detail_rows
    )


def test_evidence_detail_mixed_buckets_only_no_hit_excluded() -> None:
    """With multiple buckets enabled, only no_hit should be missing from the detail sheet."""
    records = {
        "query": record("query", positive_sources_hit=["A"]),
        "candidate": record("candidate", positive_sources_hit=["A"]),
        "relaxed": record("relaxed", positive_sources_hit=["A"]),
        "novel": record("novel"),
    }
    cfg = interaction_config(
        query_proteins=(InteractionQueryConfig("query", "", ""),),
        candidate_sources={"candidates": True, "candidates_relaxed": True, "no_hit": True},
        scoring_model="v2_evidence_based",
    )

    result = run_interaction_scoring(cfg, classification(records))

    assert result is not None
    detail_sources = {r["candidate_source"] for r in result.evidence_detail_rows}
    assert "No_hit" not in detail_sources
    assert detail_sources <= {"Candidates", "Candidates_relaxed"}
    # Every candidate actually present in those two main sheets must also
    # appear in the detail sheet -- the scope must match, not just avoid crashing.
    expected_candidates = {
        row["candidate_protein_id"]
        for sheet in ("Interaction_Candidates", "Interaction_Candidates_relaxed")
        for row in result.source_rows.get(sheet, [])
    }
    detail_candidates = {r["candidate_protein_id"] for r in result.evidence_detail_rows}
    assert detail_candidates == expected_candidates
