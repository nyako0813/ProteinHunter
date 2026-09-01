"""Tests for main.py command-line config handling."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import yaml

from config import (
    INTERACTION_ALPHAFOLD_DEFAULT,
    INTERACTION_CANDIDATE_SOURCE_DEFAULTS,
    INTERACTION_NEIGHBORHOOD_DEFAULT,
    INTERACTION_SCORING_WEIGHTS_DEFAULT,
    AnnotationTargetConfig,
    InteractionQueryConfig,
    InteractionScoringConfig,
)
from core.models import ProteinRecord
from main import build_arg_parser, _build_cdd_target_records, _records_enabled_for_annotation


def _interaction_config(
    *,
    enabled: bool,
    candidate_sources: dict[str, bool] | None = None,
    query_proteins: tuple[InteractionQueryConfig, ...] = (),
) -> SimpleNamespace:
    """Build a minimal app config exposing only interaction_scoring."""
    return SimpleNamespace(
        interaction_scoring=InteractionScoringConfig(
            enabled=enabled,
            query_proteins=query_proteins,
            query_fasta=None,
            candidate_sources=(
                dict(INTERACTION_CANDIDATE_SOURCE_DEFAULTS)
                if candidate_sources is None
                else candidate_sources
            ),
            max_candidates_per_query=200,
            include_sequences_in_excel=False,
            scoring_weights=INTERACTION_SCORING_WEIGHTS_DEFAULT,
            alphafold=INTERACTION_ALPHAFOLD_DEFAULT,
            neighborhood=INTERACTION_NEIGHBORHOOD_DEFAULT,
            scoring_model="legacy_additive",
            scoring_engine_config=None,
            functional_complementarity_ruleset=None,
            pih_evidence_bundle=None,
        )
    )


def _classification(**buckets: dict[str, ProteinRecord]) -> SimpleNamespace:
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


def test_build_cdd_target_records_extends_scope_without_mutating_candidates_dict() -> None:
    """Extending CDD's target set must not change the Candidates-sheet record count.

    This is the wiring guarantee Step 15's design relies on: main.py keeps
    using its own `records` (positive_only_records) for candidate scoring
    and the Excel "Candidates" sheet, while CDD itself is given a larger,
    separate dict when interaction_scoring needs more records annotated.
    """
    candidate = ProteinRecord(protein_id="candidate_1", sequence="MSTNPKPQR")
    no_hit = ProteinRecord(protein_id="no_hit_1", sequence="MSTNPKPQR")
    positive_only_records = {"candidate_1": candidate}
    cls = _classification(
        all_records={"candidate_1": candidate, "no_hit_1": no_hit},
        positive_only_records=positive_only_records,
        no_hit_records={"no_hit_1": no_hit},
    )
    sources = {key: False for key in INTERACTION_CANDIDATE_SOURCE_DEFAULTS}
    sources["no_hit"] = True
    config = _interaction_config(enabled=True, candidate_sources=sources)

    target_records = _build_cdd_target_records(positive_only_records, config, cls)

    assert set(target_records) == {"candidate_1", "no_hit_1"}
    # The caller's own positive_only_records dict (used later for
    # score_records / the Excel "Candidates" sheet) must stay exactly as
    # it was -- still just 1 record, not silently grown in place.
    assert positive_only_records == {"candidate_1": candidate}


def test_build_cdd_target_records_unaffected_when_interaction_scoring_disabled() -> None:
    """With interaction_scoring disabled, the CDD target set is exactly Candidates."""
    candidate = ProteinRecord(protein_id="candidate_1", sequence="MSTNPKPQR")
    positive_only_records = {"candidate_1": candidate}
    cls = _classification(
        all_records={"candidate_1": candidate},
        positive_only_records=positive_only_records,
        no_hit_records={"no_hit_1": ProteinRecord(protein_id="no_hit_1")},
    )
    config = _interaction_config(enabled=False)

    target_records = _build_cdd_target_records(positive_only_records, config, cls)

    assert set(target_records) == {"candidate_1"}


def test_main_config_arg_defaults_to_config_yaml() -> None:
    """No --config argument should keep the production config default."""
    args = build_arg_parser().parse_args([])

    assert args.config == "config.yaml"
    assert args.check_only is False


def test_main_config_arg_accepts_custom_config() -> None:
    """--config should allow a custom YAML config path."""
    args = build_arg_parser().parse_args(["--config", "config.demo.yaml"])

    assert args.config == "config.demo.yaml"


def test_main_check_only_arg_is_accepted() -> None:
    """--check-only should enable validation-only mode."""
    args = build_arg_parser().parse_args(["--check-only"])

    assert args.check_only is True


def test_main_check_only_arg_works_with_custom_config() -> None:
    """--check-only should work together with --config."""
    args = build_arg_parser().parse_args(
        ["--config", "config.demo.yaml", "--check-only"]
    )

    assert args.config == "config.demo.yaml"
    assert args.check_only is True


def test_main_check_only_works_with_directory_mode(
    tmp_path: Path,
    monkeypatch,
) -> None:
    """--check-only should validate directory mode and stop before BLAST."""
    from main import main

    class PassingStartupChecker:
        def run(self) -> bool:
            return True

    def write_source(root: Path, label: str, record_id: str) -> None:
        source_dir = root / label
        source_dir.mkdir(parents=True)
        (source_dir / "protein.faa").write_text(
            f">{record_id}\nMSTNPKPQR\n",
            encoding="utf-8",
        )

    target_dir = tmp_path / "target"
    positive_dir = tmp_path / "positive"
    negative_dir = tmp_path / "negative"
    write_source(target_dir, "target_source", "target_1")
    write_source(positive_dir, "positive_source", "positive_1")
    write_source(negative_dir, "negative_source", "negative_1")

    config_data = {
        "project": {"name": "ProteinHunter", "version": "5.0"},
        "input_mode": "directory",
        "paths": {
            "target_dir": str(target_dir),
            "positive_dir": str(positive_dir),
            "negative_dir": str(negative_dir),
            "output_excel": str(tmp_path / "out" / "results.xlsx"),
            "cache_dir": str(tmp_path / ".cache"),
            "log_dir": str(tmp_path / "logs"),
        },
        "blast": {"evalue": 1e-5, "max_target_seqs": 10, "threads": 1},
        "annotation": {
            "enable_cdd": False,
            "enable_pfam": False,
            "enable_alphafold": False,
            "enable_uniprot": False,
            "enable_gene_context": False,
            "cdd_threads": 1,
            "pfam_threads": 1,
            "alphafold_threads": 1,
        },
        "cache": {"enabled": True, "overwrite": False},
        "score": {
            "blast_weight": 5,
            "domain_weight": 4,
            "motif_weight": 3,
            "gene_context_weight": 3,
            "alphafold_weight": 2,
        },
        "logging": {"level": "INFO", "save_log": True},
    }
    config_path = tmp_path / "config.directory.yaml"
    config_path.write_text(yaml.safe_dump(config_data), encoding="utf-8")

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("main.StartupChecker", PassingStartupChecker)

    main(["--config", str(config_path), "--check-only"])

    assert (tmp_path / "data" / "temp" / "combined" / "target.combined.faa").exists()


def test_annotation_targets_select_no_hit_pfam_when_enabled() -> None:
    """No_hit records should be selected for Pfam only when enabled."""
    no_hit_record = ProteinRecord(protein_id="no_hit_1")
    record_sheets = {
        "Candidates": {},
        "Positive_all_sources": {},
        "Negative_unmatched": {},
        "No_hit": {"no_hit_1": no_hit_record},
        "Negative_hit": {},
    }
    targets = {
        "candidates": AnnotationTargetConfig(True, True, True, True),
        "positive_all_sources": AnnotationTargetConfig(True, True, True, True),
        "negative_unmatched": AnnotationTargetConfig(True, False, False, False),
        "no_hit": AnnotationTargetConfig(True, True, False, False),
        "negative_hit": AnnotationTargetConfig(True, False, False, False),
    }

    selected = _records_enabled_for_annotation(record_sheets, targets, "pfam")

    assert selected == {"no_hit_1": no_hit_record}


def test_annotation_targets_skip_no_hit_pfam_when_disabled() -> None:
    """No_hit Pfam should remain off under the safe default behavior."""
    no_hit_record = ProteinRecord(protein_id="no_hit_1")
    record_sheets = {
        "Candidates": {},
        "Positive_all_sources": {},
        "Negative_unmatched": {},
        "No_hit": {"no_hit_1": no_hit_record},
        "Negative_hit": {},
    }
    targets = {
        "candidates": AnnotationTargetConfig(True, True, True, True),
        "positive_all_sources": AnnotationTargetConfig(True, True, True, True),
        "negative_unmatched": AnnotationTargetConfig(True, False, False, False),
        "no_hit": AnnotationTargetConfig(True, False, False, False),
        "negative_hit": AnnotationTargetConfig(True, False, False, False),
    }

    selected = _records_enabled_for_annotation(record_sheets, targets, "pfam")

    assert selected == {}


def test_annotation_targets_deduplicate_records_shared_by_sheets() -> None:
    """A record in Candidates and Positive_all_sources should be annotated once."""
    shared_record = ProteinRecord(protein_id="candidate_1")
    record_sheets = {
        "Candidates": {"candidate_1": shared_record},
        "Positive_all_sources": {"candidate_1": shared_record},
        "Negative_unmatched": {},
        "No_hit": {},
        "Negative_hit": {},
    }
    targets = {
        "candidates": AnnotationTargetConfig(True, True, True, True),
        "positive_all_sources": AnnotationTargetConfig(True, True, True, True),
        "negative_unmatched": AnnotationTargetConfig(True, False, False, False),
        "no_hit": AnnotationTargetConfig(True, False, False, False),
        "negative_hit": AnnotationTargetConfig(True, False, False, False),
    }

    selected = _records_enabled_for_annotation(record_sheets, targets, "pfam")

    assert list(selected) == ["candidate_1"]
    assert selected["candidate_1"] is shared_record


def test_annotation_targets_missing_subkey_is_safe_default_off() -> None:
    """Missing target objects should not accidentally enable annotation."""
    candidate_record = ProteinRecord(protein_id="candidate_1")
    record_sheets = {
        "Candidates": {"candidate_1": candidate_record},
        "Positive_all_sources": {},
        "Negative_unmatched": {},
        "No_hit": {},
        "Negative_hit": {},
    }
    targets = {
        "candidates": SimpleNamespace(gff=True),
    }

    selected = _records_enabled_for_annotation(record_sheets, targets, "pfam")

    assert selected == {}
