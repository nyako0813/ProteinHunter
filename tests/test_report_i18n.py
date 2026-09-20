"""Tests for the report wording tables (output/report_i18n.py)."""

from __future__ import annotations

import re

import pytest

from analysis.interaction_scoring import CANDIDATE_SOURCE_MAP
from analysis.scoring_engine_config import load_scoring_engine_config
from config import INTERACTION_SCORING_WEIGHTS_DEFAULT, REPORT_LANGUAGES
from core.provenance import REFERENCE_GENOME_WARNING_TEXT
from output.report_i18n import LABEL_JA, STRINGS, SUPPORTED_LANGUAGES, bilingual, normalize_language, t
from output.word_report import category_refs_for_scoring_model

_PLACEHOLDER = re.compile(r"\{([a-z_]+)(?::[^}]*)?\}")
_JAPANESE = re.compile(r"[぀-ヿ㐀-鿿]")


def test_both_languages_define_exactly_the_same_keys() -> None:
    assert set(STRINGS["en"]) == set(STRINGS["ja"])
    assert STRINGS["en"], "empty string table"


def test_every_key_uses_the_same_placeholders_in_both_languages() -> None:
    """A placeholder present in one language only would raise (or silently drop data) when formatted."""
    for key, english in STRINGS["en"].items():
        assert set(_PLACEHOLDER.findall(english)) == set(_PLACEHOLDER.findall(STRINGS["ja"][key])), key


def test_every_japanese_string_actually_contains_japanese() -> None:
    """Catches a Japanese entry that was left as the English text (the two punctuation-only separators are exempt)."""
    for key, japanese in STRINGS["ja"].items():
        if key in {"sentence_joiner", "why.listing_separator"}:
            continue
        assert _JAPANESE.search(japanese), key


def test_supported_languages_match_the_config_validator() -> None:
    assert tuple(SUPPORTED_LANGUAGES) == tuple(REPORT_LANGUAGES)


def test_english_reference_genome_warning_is_the_one_the_excel_index_uses() -> None:
    assert STRINGS["en"]["report.reference_genome_warning"] == REFERENCE_GENOME_WARNING_TEXT


def test_t_formats_and_fails_loudly_on_a_missing_key() -> None:
    assert t("ranking.query_heading", "en", index=2, query_id="q1") == "7.2 Query: q1"
    assert t("ranking.query_heading", "ja", index=2, query_id="q1") == "7.2 クエリ: q1"
    with pytest.raises(KeyError):
        t("no.such.key", "en")


def test_unsupported_language_is_rejected_and_empty_means_default() -> None:
    assert normalize_language(None) == "en" and normalize_language("") == "en"
    with pytest.raises(ValueError, match="unsupported report language"):
        normalize_language("fr")


def test_bilingual_glosses_labels_in_japanese_only() -> None:
    assert bilingual("Strong", "en") == "Strong"
    assert bilingual("Strong", "ja") == "Strong(強)"
    assert bilingual("medium", "ja") == "medium(中)"  # data value keeps its own spelling
    assert bilingual("Candidates_relaxed", "ja") == "Candidates_relaxed(候補(緩和条件))"
    assert bilingual("Something New", "ja") == "Something New"  # unknown labels are left alone


def test_label_table_covers_every_candidate_source_the_pipeline_can_produce() -> None:
    for label, _records, _sheet in CANDIDATE_SOURCE_MAP.values():
        assert label.lower() in LABEL_JA, label


def test_label_table_covers_strengths_tiers_and_every_category_label() -> None:
    for word in ("strong", "medium", "weak", "none", "very strong", "moderate"):
        assert word in LABEL_JA, word
    engine = load_scoring_engine_config(None)
    for model in ("v2_evidence_based", "legacy_additive"):
        for ref in category_refs_for_scoring_model(model, engine, INTERACTION_SCORING_WEIGHTS_DEFAULT):
            assert ref.label.lower() in LABEL_JA, ref.label


def test_only_the_sub_bucket_sources_lack_a_narrative_sentence() -> None:
    """The narrative has a sentence per candidate_source except the three strength sub-buckets (as before this change)."""
    with_sentence = {key.split(".", 1)[1] for key in STRINGS["en"] if key.startswith("source.")}
    all_sources = {label for label, _records, _sheet in CANDIDATE_SOURCE_MAP.values()}

    assert all_sources - with_sentence == {"Negative_strong_hit", "Negative_medium_hit", "Negative_weak_hit"}
    assert with_sentence <= all_sources


def test_hedged_sentences_stay_hedged_in_japanese() -> None:
    """The design forbids asserting that a candidate IS the partner; the Japanese must keep the hedge."""
    ja = STRINGS["ja"]
    assert "確定的に同定したものではない" in ja["interpretation.opening"]
    assert "立証するものではない" in ja["color.genomic_context"]
    assert "立証できない" in ja["color.functional_domain"]
    assert "意味しない" in ja["negative_evidence_closer"]
