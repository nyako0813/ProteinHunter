"""Tests for analysis/pih_evidence_bridge.py against PIH's real record shape.

Records come from tests/pih_fixtures.py, which follows PIH's own schema:
``integrated_scoring`` is a list of zero or one IntegratedScore. The bridge used
to expect a dict there, so it never imported a single category from a real PIH
bundle; ``test_real_schema_record_fills_the_bridged_categories`` is the direct
regression test for that.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from pih_fixtures import bundle_record, integrated_score, score_category, write_bundle

from analysis.pih_evidence_bridge import BRIDGED_PIH_CATEGORIES, load_pih_evidence_bundle

QUERY = "WP_011024006.1"
CANDIDATE = "WP_011024007.1"


def scored_record(categories: list[dict], *, query: str = QUERY, candidate: str = CANDIDATE) -> dict:
    return bundle_record(query, candidate, [integrated_score(query, candidate, categories)])


def load(tmp_path: Path, records: list[dict]):
    return load_pih_evidence_bundle(write_bundle(tmp_path / "candidate_evidence_bundle.jsonl", records))


def test_real_schema_record_fills_the_bridged_categories(tmp_path: Path) -> None:
    """A record in PIH's actual shape yields the three bridged categories -- not an empty result."""
    bundle = load(
        tmp_path,
        [
            scored_record(
                [
                    score_category("genomic_context", 0.8, 1.0),
                    score_category("functional_annotation", 0.4, 1.0),
                    score_category("cellular_compatibility", 1.0, 0.5),
                    score_category("evolutionary", 0.6, 1.75),
                    score_category("direct_interaction", 0.9, 1.5),
                ]
            )
        ],
    )

    found = bundle.lookup([QUERY], [CANDIDATE])

    assert bundle.warnings == ()
    assert set(found) == set(BRIDGED_PIH_CATEGORIES)
    assert (found["evolutionary"].normalized_score, found["evolutionary"].available_weight) == (0.6, 1.75)
    assert (found["cellular_compatibility"].normalized_score, found["cellular_compatibility"].available_weight) == (
        1.0,
        0.5,
    )
    assert found["direct_interaction"].normalized_score == 0.9


def test_only_the_three_non_overlapping_categories_are_bridged(tmp_path: Path) -> None:
    bundle = load(
        tmp_path,
        [
            scored_record(
                [
                    score_category("genomic_context", 1.0),
                    score_category("functional_annotation", 1.0),
                    score_category("evolutionary", 1.0),
                ]
            )
        ],
    )

    assert set(bundle.lookup([QUERY], [CANDIDATE])) == {"evolutionary"}


def test_categories_without_available_weight_are_not_active(tmp_path: Path) -> None:
    """PIH lists all five categories on every scored pair; weight 0 means no evidence for that one."""
    bundle = load(
        tmp_path,
        [
            scored_record(
                [
                    score_category("cellular_compatibility", 0.0, 0.0),
                    score_category("evolutionary", 0.5, 1.0),
                    score_category("direct_interaction", 0.0, 0.0),
                ]
            )
        ],
    )

    assert set(bundle.lookup([QUERY], [CANDIDATE])) == {"evolutionary"}


def test_empty_integrated_scoring_list_means_pih_did_not_score_the_pair(tmp_path: Path) -> None:
    """The graceful path: no entry, no evidence, no warning, and the rest of the bundle still loads."""
    bundle = load(
        tmp_path,
        [
            bundle_record(QUERY, "WP_011024008.1", None),
            scored_record([score_category("evolutionary", 0.5, 1.0)]),
        ],
    )

    assert bundle.warnings == ()
    assert bundle.lookup([QUERY], ["WP_011024008.1"]) == {}
    assert set(bundle.lookup([QUERY], [CANDIDATE])) == {"evolutionary"}


def test_missing_integrated_scoring_key_is_treated_like_an_empty_list(tmp_path: Path) -> None:
    record = bundle_record(QUERY, CANDIDATE, None)
    del record["integrated_scoring"]

    bundle = load(tmp_path, [record])

    assert bundle.pairs == {}
    assert bundle.warnings == ()


def test_the_separate_score_object_is_never_read_as_a_category_breakdown(tmp_path: Path) -> None:
    """PIH's `score` (CandidateScore) is a different structure; the bridge must not fall back to it."""
    record = bundle_record(QUERY, CANDIDATE, None)
    record["score"] = {"category_scores": [score_category("evolutionary", 1.0, 1.0)]}

    bundle = load(tmp_path, [record])

    assert bundle.lookup([QUERY], [CANDIDATE]) == {}


def test_old_dict_shape_is_not_accepted_and_is_reported(tmp_path: Path) -> None:
    """A dict `integrated_scoring` is not PIH output. It is skipped -- loudly, not silently."""
    record = bundle_record(QUERY, CANDIDATE, None)
    record["integrated_scoring"] = {"category_scores": [score_category("evolutionary", 1.0, 1.0)]}

    bundle = load(tmp_path, [record])

    assert bundle.pairs == {}
    assert len(bundle.warnings) == 1
    assert "not a list" in bundle.warnings[0]


def test_shape_problems_are_aggregated_into_one_warning_per_kind(tmp_path: Path) -> None:
    """A wholesale mismatch (every line of a large bundle) must not flood the log."""
    broken = [dict(bundle_record(QUERY, f"WP_0000{i}.1", None), integrated_scoring={"x": 1}) for i in range(3)]
    good = scored_record([score_category("evolutionary", 0.5, 1.0)])

    bundle = load(tmp_path, [good, *broken])

    assert len(bundle.warnings) == 1
    assert "3 record(s)" in bundle.warnings[0]
    assert "first at line 2" in bundle.warnings[0]
    assert set(bundle.lookup([QUERY], [CANDIDATE])) == {"evolutionary"}  # the good record still loads


def test_entry_without_category_scores_is_reported(tmp_path: Path) -> None:
    record = bundle_record(QUERY, CANDIDATE, [{"query_protein_id": QUERY, "candidate_protein_id": CANDIDATE}])

    bundle = load(tmp_path, [record])

    assert bundle.pairs == {}
    assert any("category_scores" in warning for warning in bundle.warnings)


def test_more_than_one_integrated_scoring_entry_warns_and_uses_the_first(tmp_path: Path) -> None:
    first = integrated_score(QUERY, CANDIDATE, [score_category("evolutionary", 0.25, 1.0)])
    second = integrated_score(QUERY, CANDIDATE, [score_category("evolutionary", 0.75, 1.0)])

    bundle = load(tmp_path, [bundle_record(QUERY, CANDIDATE, [first, second])])

    assert bundle.lookup([QUERY], [CANDIDATE])["evolutionary"].normalized_score == 0.25
    assert any("more than one" in warning for warning in bundle.warnings)


# --- negative normalized_score: sign information is currently discarded ---------


def test_negative_normalized_score_is_clamped_to_zero_but_stays_available(tmp_path: Path) -> None:
    """PIH's category score is in [-1, 1]. The bridge keeps the [0, 1] clamp (no penalty is modelled).

    A negative value (e.g. an incompatible localization) becomes 0.0 evidence:
    the category is still active -- it counts toward the score's denominator --
    but contributes nothing. The sign is discarded on purpose until real bundles
    show how negative values are distributed.
    """
    bundle = load(
        tmp_path,
        [
            scored_record(
                [
                    score_category("cellular_compatibility", -0.25, 0.5),
                    score_category("evolutionary", -1.0, 1.75),
                    score_category("direct_interaction", 0.5, 1.5),
                ]
            )
        ],
    )

    found = bundle.lookup([QUERY], [CANDIDATE])

    assert found["cellular_compatibility"].normalized_score == 0.0
    assert found["cellular_compatibility"].available_weight == 0.5
    assert found["evolutionary"].normalized_score == 0.0
    assert found["direct_interaction"].normalized_score == 0.5  # positives are untouched


def test_normalized_score_above_one_is_clamped_to_one(tmp_path: Path) -> None:
    bundle = load(tmp_path, [scored_record([score_category("evolutionary", 1.4, 1.0)])])

    assert bundle.lookup([QUERY], [CANDIDATE])["evolutionary"].normalized_score == 1.0


# --- graceful degradation is unchanged --------------------------------------------


def test_missing_file_and_bad_json_still_degrade_gracefully(tmp_path: Path) -> None:
    missing = load_pih_evidence_bundle(tmp_path / "nope.jsonl")
    assert missing.pairs == {} and "not found" in missing.warnings[0]

    path = tmp_path / "candidate_evidence_bundle.jsonl"
    good = json.dumps(scored_record([score_category("evolutionary", 0.5, 1.0)]))
    path.write_text(f"not json\n\n[1, 2]\n{good}\n", encoding="utf-8")
    bundle = load_pih_evidence_bundle(path)

    assert set(bundle.lookup([QUERY], [CANDIDATE])) == {"evolutionary"}
    assert any("not valid JSON" in warning for warning in bundle.warnings)
    assert any("not a JSON object" in warning for warning in bundle.warnings)


def test_non_numeric_category_score_is_skipped_with_a_warning(tmp_path: Path) -> None:
    record = scored_record([score_category("evolutionary", 0.5, 1.0), score_category("direct_interaction", 0.5, 1.0)])
    record["integrated_scoring"][0]["category_scores"][0]["normalized_score"] = "high"

    bundle = load(tmp_path, [record])

    assert set(bundle.lookup([QUERY], [CANDIDATE])) == {"direct_interaction"}
    assert any("non-numeric" in warning for warning in bundle.warnings)


@pytest.mark.parametrize("query_key,candidate_key", [(QUERY, CANDIDATE), ("WP_011024006", "WP_011024007")])
def test_lookup_tries_each_spelling_in_order(tmp_path: Path, query_key: str, candidate_key: str) -> None:
    bundle = load(tmp_path, [scored_record([score_category("evolutionary", 0.5, 1.0)], query=query_key, candidate=candidate_key)])

    assert bundle.lookup(["", "MA_4115", QUERY, "WP_011024006"], ["", "MA_4116", CANDIDATE, "WP_011024007"])
