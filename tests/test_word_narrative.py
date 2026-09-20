"""Tests for the deterministic Word-report narrative text (Phase 6-8 Stage 2)."""

from __future__ import annotations

from output.word_narrative import (
    CategoryRef,
    build_biological_interpretation,
    build_evolutionary_closer,
    build_why_ranks_highly,
)

V2_REFS: tuple[CategoryRef, ...] = (
    CategoryRef("candidate_priority_score", "sequence", "Sequence/Source Classification", 30.0),
    CategoryRef("same_gene_neighborhood_score", "genomic_context", "Genomic Context", 25.0),
    CategoryRef("functional_domain_score", "functional_domain", "Functional/Domain", 20.0),
    CategoryRef("interaction_evidence_score", "interaction", "Interaction", 47.0),
    CategoryRef("evolutionary_score", "evolutionary", "Evolutionary", 10.0),
    CategoryRef("cellular_compatibility_score", "cellular_compatibility", "Cellular Compatibility", 5.0),
)


def _row(**overrides: object) -> dict:
    row = {
        "query_id": "query_1",
        "candidate_protein_id": "candidate_1",
        "candidate_source": "Candidates",
        "negative_hit_strength": "none",
        "final_score": 42.884,
        "final_score_tier": "Tier3_Moderate",
        "evidence_category_count": 2,
        "candidate_priority_score": 30.0,
        "same_gene_neighborhood_score": 12.5,
        "functional_domain_score": None,
        "interaction_evidence_score": None,
        "evolutionary_score": None,
        "cellular_compatibility_score": None,
    }
    row.update(overrides)
    return row


def test_why_ranks_highly_reproducible_for_identical_input() -> None:
    """Same row + same refs must produce byte-identical text (design spec section 45)."""
    row = _row()
    first = build_why_ranks_highly(row, V2_REFS)
    second = build_why_ranks_highly(dict(row), V2_REFS)
    assert first == second


def test_why_ranks_highly_tier1_opening() -> None:
    row = _row(final_score_tier="Tier1_VeryStrong", final_score=88.4, evidence_category_count=4)
    text = build_why_ranks_highly(row, V2_REFS)
    assert "Tier 1" in text
    assert "88.4/100" in text
    assert "4 independent" in text


def test_why_ranks_highly_unclassified_opening_has_no_final_score_sentence() -> None:
    row = _row(final_score_tier="Unclassified", final_score=None)
    text = build_why_ranks_highly(row, V2_REFS)
    assert "did not receive a formal Final Score" in text
    assert "None" not in text


def test_why_ranks_highly_lists_only_contributing_categories() -> None:
    row = _row(candidate_priority_score=30.0, same_gene_neighborhood_score=0.0)
    text = build_why_ranks_highly(row, V2_REFS)
    assert "Sequence/Source Classification (30.0/30)" in text
    assert "Genomic Context (0.0/25)" not in text


def test_why_ranks_highly_no_contributing_categories() -> None:
    row = _row(candidate_priority_score=None, same_gene_neighborhood_score=None)
    text = build_why_ranks_highly(row, V2_REFS)
    assert "No individual evidence category scored above zero" in text


def test_why_ranks_highly_candidate_source_sentences_are_distinct() -> None:
    sources = (
        "Candidates",
        "Positive_all_sources",
        "Candidates_relaxed",
        "No_hit",
        "Negative_unmatched",
        "Negative_hit",
    )
    texts = {source: build_why_ranks_highly(_row(candidate_source=source), V2_REFS) for source in sources}
    assert len(set(texts.values())) == len(sources)
    assert "strict positive candidate" in texts["Candidates"]
    assert "lineage-specific or novel" in texts["No_hit"]
    assert "included in this ranking despite that hit" in texts["Negative_hit"]


def test_why_ranks_highly_negative_hit_strength_caveat_outside_negative_bucket() -> None:
    row = _row(candidate_source="Candidates_relaxed", negative_hit_strength="weak")
    text = build_why_ranks_highly(row, V2_REFS)
    assert "weak BLAST hit to a negative-reference sequence" in text


def test_why_ranks_highly_no_duplicate_caveat_for_negative_hit_bucket() -> None:
    row = _row(candidate_source="Negative_hit", negative_hit_strength="strong")
    text = build_why_ranks_highly(row, V2_REFS)
    assert text.count("negative-reference") == 1


def test_evolutionary_closer_legacy_model() -> None:
    text = build_evolutionary_closer("legacy_additive", pih_bundle_configured=False)
    assert "legacy_additive scoring model" in text
    assert "does not evaluate Evolutionary" in text


def test_evolutionary_closer_v2_without_bundle() -> None:
    text = build_evolutionary_closer("v2_evidence_based", pih_bundle_configured=False)
    assert "was not supplied" in text
    assert "neither category was evaluated" in text


def test_evolutionary_closer_v2_with_bundle() -> None:
    text = build_evolutionary_closer("v2_evidence_based", pih_bundle_configured=True)
    assert "were evaluated for this run" in text
    assert "not that the category was skipped" in text


def test_biological_interpretation_never_asserts_identity() -> None:
    """Design spec section 35: must hedge, never assert the candidate IS the target."""
    row = _row(candidate_priority_score=30.0, same_gene_neighborhood_score=12.5)
    closer = build_evolutionary_closer("v2_evidence_based", pih_bundle_configured=False)
    text = build_biological_interpretation(row, rank=1, n_candidates=42, category_refs=V2_REFS, evolutionary_closer=closer)
    assert "not a confirmed identification" in text
    assert "best-supported by currently available evidence" in text
    assert "ranked 1 of 42" in text


def test_biological_interpretation_includes_color_sentence_for_evaluated_category() -> None:
    row = _row(same_gene_neighborhood_score=12.5)
    closer = build_evolutionary_closer("v2_evidence_based", pih_bundle_configured=False)
    text = build_biological_interpretation(row, rank=1, n_candidates=10, category_refs=V2_REFS, evolutionary_closer=closer)
    assert "genomic proximity" in text


def test_biological_interpretation_omits_color_sentence_for_unevaluated_category() -> None:
    row = _row(evolutionary_score=None)
    closer = build_evolutionary_closer("v2_evidence_based", pih_bundle_configured=False)
    text = build_biological_interpretation(row, rank=1, n_candidates=10, category_refs=V2_REFS, evolutionary_closer=closer)
    assert "Evolutionary/phylogenetic profile evidence" not in text


def test_biological_interpretation_always_states_negative_evidence_is_reserved() -> None:
    row = _row()
    closer = build_evolutionary_closer("v2_evidence_based", pih_bundle_configured=True)
    text = build_biological_interpretation(row, rank=1, n_candidates=10, category_refs=V2_REFS, evolutionary_closer=closer)
    assert "reserved category with no implemented signal" in text


def test_biological_interpretation_reproducible_for_identical_input() -> None:
    row = _row()
    closer = build_evolutionary_closer("v2_evidence_based", pih_bundle_configured=False)
    first = build_biological_interpretation(row, rank=3, n_candidates=15, category_refs=V2_REFS, evolutionary_closer=closer)
    second = build_biological_interpretation(dict(row), rank=3, n_candidates=15, category_refs=V2_REFS, evolutionary_closer=closer)
    assert first == second


# ---------------------------------------------------------------------------
# Japanese wording (language="ja")
# ---------------------------------------------------------------------------


def test_default_language_is_english_and_ja_is_opt_in() -> None:
    row = _row()

    assert build_why_ranks_highly(row, V2_REFS) == build_why_ranks_highly(row, V2_REFS, "en")
    assert build_why_ranks_highly(row, V2_REFS, "ja") != build_why_ranks_highly(row, V2_REFS, "en")


def test_ja_why_ranks_highly_keeps_numbers_and_glosses_labels() -> None:
    row = _row(final_score_tier="Tier1_VeryStrong", final_score=88.4, evidence_category_count=4)

    text = build_why_ranks_highly(row, V2_REFS, "ja")

    assert "88.4/100" in text and "4 個" in text
    assert "Tier 1 — Very Strong(非常に強い)" in text
    assert "Sequence/Source Classification(配列・ソース分類) (30.0/30)" in text
    assert "Genomic Context(ゲノム文脈) (12.5/25)" in text


def test_ja_negative_hit_caveat_glosses_the_strength_but_keeps_its_value() -> None:
    row = _row(candidate_source="Candidates_relaxed", negative_hit_strength="medium")

    text = build_why_ranks_highly(row, V2_REFS, "ja")

    assert "medium(中)" in text
    assert "Candidates_relaxed(候補(緩和条件))" in text
    assert "candidate_source" in text  # column name stays English


def test_ja_unknown_candidate_source_adds_no_sentence_like_english() -> None:
    row = _row(candidate_source="Negative_strong_hit")

    en = build_why_ranks_highly(row, V2_REFS, "en")
    ja = build_why_ranks_highly(row, V2_REFS, "ja")

    assert "Negative_strong_hit" not in en and "Negative_strong_hit" not in ja


def test_ja_biological_interpretation_stays_hedged_and_names_the_candidate() -> None:
    row = _row(same_gene_neighborhood_score=12.5)
    closer = build_evolutionary_closer("v2_evidence_based", False, "ja")

    text = build_biological_interpretation(row, 3, 20, V2_REFS, closer, "ja")

    assert "candidate_1" in text and "query_1" in text and "20 個の候補のうち 3 位" in text
    assert "確定的に同定したものではない" in text
    assert "立証するものではない" in text  # genomic-context colour sentence
    assert "評価されていない" in text  # the PIH-not-supplied closer


def test_ja_evolutionary_closer_covers_all_three_branches() -> None:
    legacy = build_evolutionary_closer("legacy_additive", False, "ja")
    with_pih = build_evolutionary_closer("v2_evidence_based", True, "ja")
    without_pih = build_evolutionary_closer("v2_evidence_based", False, "ja")

    assert len({legacy, with_pih, without_pih}) == 3
    assert "legacy_additive" in legacy and "PIH" in with_pih and "提供されなかった" in without_pih


def test_ja_output_is_reproducible() -> None:
    row = _row()

    assert build_why_ranks_highly(row, V2_REFS, "ja") == build_why_ranks_highly(dict(row), V2_REFS, "ja")


def test_unsupported_language_raises() -> None:
    import pytest

    with pytest.raises(ValueError):
        build_why_ranks_highly(_row(), V2_REFS, "fr")

