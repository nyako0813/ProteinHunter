"""Tests for lightweight interaction candidate ranking."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from config import (
    INTERACTION_ALPHAFOLD_DEFAULT,
    INTERACTION_CANDIDATE_SOURCE_DEFAULTS,
    INTERACTION_NEIGHBORHOOD_DEFAULT,
    INTERACTION_SCORING_WEIGHTS_DEFAULT,
    InteractionNeighborhoodConfig,
    InteractionQueryConfig,
    InteractionScoringConfig,
)
from core.models import DomainHit, ProteinRecord
from analysis.interaction_scoring import (
    interaction_pair_columns,
    run_interaction_scoring,
)


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
) -> SimpleNamespace:
    """Build a minimal app config for interaction scoring tests."""
    return SimpleNamespace(
        paths=SimpleNamespace(gff_file=gff_file),
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
