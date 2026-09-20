"""Tests for run provenance (core/provenance.py)."""

from __future__ import annotations

import dataclasses
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

from config import load_config
from core import provenance as provenance_module
from core.constants import APP_VERSION
from core.provenance import (
    CONFIG_HASH_LENGTH,
    RunProvenance,
    code_version_text,
    collect_run_provenance,
    compute_config_hash,
    effective_scoring_parameters,
    provenance_sidecar_filename,
    provenance_sidecar_path,
    write_provenance_sidecar,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
GIT = shutil.which("git")


def fake_config(**scoring: Path | None) -> SimpleNamespace:
    return SimpleNamespace(interaction_scoring=SimpleNamespace(**scoring))


def git(cwd: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@example.com", *args],
        cwd=cwd,
        check=True,
        capture_output=True,
    )


# ---------------------------------------------------------------------------
# git provenance
# ---------------------------------------------------------------------------


def test_collect_without_a_git_checkout_returns_none_not_an_exception(tmp_path: Path) -> None:
    result = collect_run_provenance(fake_config(), tmp_path)

    assert result.git_commit is None
    assert result.git_dirty is None
    assert result.app_version == APP_VERSION
    assert len(result.config_hash) == CONFIG_HASH_LENGTH


def test_collect_when_git_binary_is_missing_returns_none(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def no_git(*args, **kwargs):
        raise FileNotFoundError("git")

    monkeypatch.setattr(provenance_module.subprocess, "run", no_git)

    result = collect_run_provenance(fake_config(), tmp_path)

    assert (result.git_commit, result.git_dirty) == (None, None)


@pytest.mark.skipif(GIT is None, reason="git not installed")
def test_collect_reads_commit_and_dirty_flag_from_a_real_repository(tmp_path: Path) -> None:
    (tmp_path / "tracked.txt").write_text("one\n")
    git(tmp_path, "init", "-q")
    git(tmp_path, "add", "tracked.txt")
    git(tmp_path, "commit", "-q", "-m", "first")

    clean = collect_run_provenance(fake_config(), tmp_path)
    assert clean.git_commit is not None and len(clean.git_commit) == 7
    assert clean.git_dirty is False

    (tmp_path / "untracked.txt").write_text("x\n")  # untracked files do not count
    assert collect_run_provenance(fake_config(), tmp_path).git_dirty is False

    (tmp_path / "tracked.txt").write_text("two\n")
    dirty = collect_run_provenance(fake_config(), tmp_path)
    assert dirty.git_dirty is True
    assert dirty.git_commit == clean.git_commit


@pytest.mark.skipif(GIT is None, reason="git not installed")
def test_a_subdirectory_of_some_other_repository_is_not_reported_as_its_commit(tmp_path: Path) -> None:
    """An unpacked release zip sitting inside an unrelated repo must not borrow that repo's SHA."""
    (tmp_path / "f.txt").write_text("x\n")
    git(tmp_path, "init", "-q")
    git(tmp_path, "add", "f.txt")
    git(tmp_path, "commit", "-q", "-m", "first")
    nested = tmp_path / "unpacked"
    nested.mkdir()

    result = collect_run_provenance(fake_config(), nested)

    assert (result.git_commit, result.git_dirty) == (None, None)


def test_code_version_text_formats_commit_dirty_and_unknown() -> None:
    now = datetime(2026, 9, 20, 12, 0)

    def make(commit: str | None, dirty: bool | None) -> RunProvenance:
        return RunProvenance("5.0", commit, dirty, "abcd1234", None, now)

    assert code_version_text(make("1a2b3c4", False)) == "5.0 (git 1a2b3c4)"
    assert code_version_text(make("1a2b3c4", True)) == "5.0 (git 1a2b3c4+dirty)"
    assert code_version_text(make(None, None)) == "5.0 (git unknown)"


# ---------------------------------------------------------------------------
# config hash
# ---------------------------------------------------------------------------


def test_config_hash_is_deterministic_and_tracks_file_content(tmp_path: Path) -> None:
    config_yaml = tmp_path / "config.yaml"
    config_yaml.write_text("blast:\n  evalue: 1e-5\n")
    inputs = [("config.yaml", config_yaml)]

    first = compute_config_hash(inputs)
    assert compute_config_hash(inputs) == first

    config_yaml.write_text("blast:\n  evalue: 1e-3\n")
    assert compute_config_hash(inputs) != first


def test_config_hash_covers_the_external_scoring_yaml_files(tmp_path: Path) -> None:
    config_yaml = tmp_path / "config.yaml"
    config_yaml.write_text("a: 1\n")
    engine = tmp_path / "engine.yaml"
    engine.write_text("caps: 1\n")
    rules = tmp_path / "rules.yaml"
    rules.write_text("rules: 1\n")
    family_map = tmp_path / "family.yaml"
    family_map.write_text("map: 1\n")
    config = fake_config(
        scoring_engine_config=engine,
        functional_complementarity_ruleset=rules,
        domain_family_map_path=family_map,
    )

    def hash_now() -> str:
        return collect_run_provenance(config, tmp_path, config_path=config_yaml).config_hash

    baseline = hash_now()
    for path in (engine, rules, family_map):
        original = path.read_text()
        path.write_text(original + "# changed\n")
        assert hash_now() != baseline, path.name
        path.write_text(original)
        assert hash_now() == baseline


def test_config_hash_ignores_line_endings_but_not_labels_or_placeholders(tmp_path: Path) -> None:
    lf = tmp_path / "lf.yaml"
    lf.write_bytes(b"a: 1\nb: 2\n")
    crlf = tmp_path / "crlf.yaml"
    crlf.write_bytes(b"a: 1\r\nb: 2\r\n")

    assert compute_config_hash([("config.yaml", lf)]) == compute_config_hash([("config.yaml", crlf)])
    assert compute_config_hash([("config.yaml", lf)]) != compute_config_hash([("other", lf)])
    # default (no path) vs a missing file are different states
    assert compute_config_hash([("x", None)]) != compute_config_hash([("x", tmp_path / "missing.yaml")])


def test_reference_genome_check_result_is_recorded_as_given(tmp_path: Path) -> None:
    for value in (True, False, None):
        result = collect_run_provenance(fake_config(), tmp_path, reference_genome_check_passed=value)
        assert result.reference_genome_check_passed is value


def test_generated_at_defaults_to_now_and_can_be_pinned(tmp_path: Path) -> None:
    pinned = datetime(2026, 1, 2, 3, 4)

    assert collect_run_provenance(fake_config(), tmp_path, generated_at=pinned).generated_at == pinned
    assert isinstance(collect_run_provenance(fake_config(), tmp_path).generated_at, datetime)


# ---------------------------------------------------------------------------
# sidecar
# ---------------------------------------------------------------------------


def test_sidecar_is_named_after_the_excel_stem_in_the_same_directory(tmp_path: Path) -> None:
    excel = tmp_path / "out" / "MA_4115_2.xlsx"

    assert provenance_sidecar_filename(excel) == "MA_4115_2.run_provenance.yaml"
    assert provenance_sidecar_path(excel) == tmp_path / "out" / "MA_4115_2.run_provenance.yaml"


def _config_copy(tmp_path: Path, name: str, evalue: str) -> Path:
    """A copy of the repo's config.yaml with one value changed -- an uncommitted local edit."""
    text = (REPO_ROOT / "config.yaml").read_text(encoding="utf-8")
    assert "evalue: 1e-5" in text
    path = tmp_path / name
    path.write_text(text.replace("evalue: 1e-5", f"evalue: {evalue}", 1), encoding="utf-8")
    return path


def test_sidecar_records_the_effective_config_of_an_uncommitted_edit(tmp_path: Path) -> None:
    """The scenario the feature exists for: two runs differing by one local config edit."""
    results = []
    for name, evalue in (("a.yaml", "1e-5"), ("b.yaml", "1e-3")):
        config_path = _config_copy(tmp_path, name, evalue)
        config = load_config(config_path, initialize=False)
        provenance = collect_run_provenance(config, tmp_path, config_path=config_path)
        sidecar = write_provenance_sidecar(
            provenance, config, tmp_path / f"{name}.xlsx", config_path=config_path
        )
        results.append((provenance, yaml.safe_load(sidecar.read_text(encoding="utf-8"))))

    (prov_a, side_a), (prov_b, side_b) = results
    assert prov_a.config_hash != prov_b.config_hash
    assert side_a["effective_config"]["blast"]["evalue"] == pytest.approx(1e-5)
    assert side_b["effective_config"]["blast"]["evalue"] == pytest.approx(1e-3)
    assert side_b["run_provenance"]["config_hash"] == prov_b.config_hash
    assert side_b["config_file"].endswith("b.yaml")


def test_sidecar_is_plain_yaml_with_paths_and_tuples_converted(tmp_path: Path) -> None:
    config_path = _config_copy(tmp_path, "c.yaml", "1e-5")
    config = load_config(config_path, initialize=False)
    provenance = collect_run_provenance(config, tmp_path, reference_genome_check_passed=False)

    sidecar = write_provenance_sidecar(provenance, config, tmp_path / "run.xlsx")
    text = sidecar.read_text(encoding="utf-8")

    assert "!!python" not in text  # nothing that safe_load would refuse
    payload = yaml.safe_load(text)
    assert payload["run_provenance"]["reference_genome_check_passed"] is False
    assert isinstance(payload["effective_config"]["paths"]["output_excel"], str)
    assert isinstance(payload["effective_config"]["interaction_scoring"]["query_proteins"], list)
    assert payload["config_file"] is None


# ---------------------------------------------------------------------------
# effective scoring parameters (built-in code defaults are part of the hash)
# ---------------------------------------------------------------------------


def _hash_with_default_config(tmp_path: Path) -> str:
    return collect_run_provenance(fake_config(scoring_engine_config=None), tmp_path).config_hash


def test_config_hash_changes_when_a_built_in_tier_threshold_default_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """No YAML file is involved here: only the code default differs."""
    from analysis import scoring_engine_config as engine_config

    before = _hash_with_default_config(tmp_path)
    assert _hash_with_default_config(tmp_path) == before

    edited = dataclasses.replace(
        engine_config.DEFAULT_SCORING_ENGINE_CONFIG,
        tiers=dataclasses.replace(engine_config.DEFAULT_SCORING_ENGINE_CONFIG.tiers, tier3_min_score=25.0),
    )
    monkeypatch.setattr(engine_config, "DEFAULT_SCORING_ENGINE_CONFIG", edited)

    assert _hash_with_default_config(tmp_path) != before


def test_config_hash_changes_when_a_component_weight_changes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """V2_COMPONENT_WEIGHTS is not YAML-configurable; it must still be fingerprinted."""
    from analysis import interaction_scoring

    before = _hash_with_default_config(tmp_path)
    monkeypatch.setitem(interaction_scoring.V2_COMPONENT_WEIGHTS, "coexpression_gse64349", 0.5)

    assert _hash_with_default_config(tmp_path) != before


def test_effective_scoring_parameters_report_defaults_and_file_overrides(tmp_path: Path) -> None:
    defaults = effective_scoring_parameters(fake_config(scoring_engine_config=None))
    assert defaults["scoring_engine"]["tiers"]["tier3_min_score"] == 35.0
    assert defaults["scoring_engine"]["category_caps"]["interaction"] == 70.0
    assert defaults["v2_component_weights"]["coexpression_gse64349"] == pytest.approx(1 / 3)

    override = tmp_path / "engine.yaml"
    override.write_text("category_caps:\n  source_classification: 30\ntiers:\n  tier3_min_score: 40\n")
    overridden = effective_scoring_parameters(fake_config(scoring_engine_config=override))
    assert overridden["scoring_engine"]["tiers"]["tier3_min_score"] == 40.0


def test_effective_scoring_parameters_survive_an_unparsable_engine_file(tmp_path: Path) -> None:
    broken = tmp_path / "engine.yaml"
    broken.write_text("tiers: [not, a, mapping\n")

    params = effective_scoring_parameters(fake_config(scoring_engine_config=broken))
    result = collect_run_provenance(fake_config(scoring_engine_config=broken), tmp_path)

    assert params["scoring_engine"] == "<unavailable>"
    assert len(result.config_hash) == CONFIG_HASH_LENGTH


def test_compute_config_hash_accepts_in_memory_bytes(tmp_path: Path) -> None:
    assert compute_config_hash([("x", b"1")]) == compute_config_hash([("x", b"1")])
    assert compute_config_hash([("x", b"1")]) != compute_config_hash([("x", b"2")])
    assert compute_config_hash([("x", b"<default>")]) == compute_config_hash([("x", None)])


def test_sidecar_records_the_effective_scoring_parameters(tmp_path: Path) -> None:
    config = fake_config(scoring_engine_config=None)
    provenance = collect_run_provenance(config, tmp_path)

    sidecar = write_provenance_sidecar(provenance, config, tmp_path / "run.xlsx")
    payload = yaml.safe_load(sidecar.read_text(encoding="utf-8"))

    assert payload["effective_scoring_parameters"]["scoring_engine"]["tiers"]["tier3_min_score"] == 35.0
    assert "coexpression_gse64349" in payload["effective_scoring_parameters"]["v2_component_weights"]
