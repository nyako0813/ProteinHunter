"""Tests for configuration validation."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
import yaml

from config import load_config
from core.exceptions import ConfigError


def valid_config_data() -> dict[str, Any]:
    """Return a complete valid config dictionary for tests."""
    return {
        "project": {
            "name": "ProteinHunter",
            "version": "5.0",
        },
        "paths": {
            "target_fasta": "./data/input/target.faa",
            "positive_fasta": "./data/databases/positive.faa",
            "negative_fasta": "./data/databases/negative.faa",
            "gff": "./data/input/genome.gff",
            "output_excel": "./data/output/results.xlsx",
            "cache_dir": "./.cache",
            "log_dir": "./logs",
        },
        "blast": {
            "evalue": 1e-5,
            "max_target_seqs": 10,
            "threads": "auto",
        },
        "annotation": {
            "enable_cdd": True,
            "enable_pfam": True,
            "enable_alphafold": True,
            "enable_uniprot": True,
            "enable_gene_context": False,
            "cdd_threads": 1,
            "pfam_threads": 1,
            "alphafold_threads": 1,
        },
        "cache": {
            "enabled": True,
            "overwrite": False,
        },
        "score": {
            "blast_weight": 5,
            "domain_weight": 4,
            "motif_weight": 3,
            "gene_context_weight": 3,
            "alphafold_weight": 2,
        },
        "logging": {
            "level": "INFO",
            "save_log": True,
        },
    }


def write_config(tmp_path: Path, data: dict[str, Any]) -> Path:
    """Write a temporary YAML config file."""
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return config_path


def test_valid_config_passes(tmp_path: Path) -> None:
    """A complete valid config should load successfully."""
    config_path = write_config(tmp_path, valid_config_data())

    cfg = load_config(config_path, initialize=False)

    assert cfg.project_name == "ProteinHunter"
    assert cfg.input_mode == "file"
    assert cfg.blast.evalue == 1e-5
    assert cfg.blast.threads >= 1
    assert cfg.paths.gff_file == Path("./data/input/genome.gff")
    assert cfg.annotation.pfam_evalue_threshold == 1e-5


def test_optional_pfam_evalue_threshold_can_be_set(tmp_path: Path) -> None:
    """Pfam e-value threshold should be configurable when present."""
    data = valid_config_data()
    data["annotation"]["pfam_evalue_threshold"] = 1e-20
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.annotation.pfam_evalue_threshold == 1e-20


def test_missing_pfam_evalue_threshold_uses_default(tmp_path: Path) -> None:
    """Missing Pfam e-value threshold should remain backward compatible."""
    data = valid_config_data()
    data["annotation"].pop("pfam_evalue_threshold", None)
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.annotation.pfam_evalue_threshold == 1e-5


def test_missing_annotation_targets_uses_safe_defaults(tmp_path: Path) -> None:
    """Missing annotation_targets should preserve the default annotation plan."""
    data = valid_config_data()
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.annotation_targets["candidates"].pfam is True
    assert cfg.annotation_targets["positive_all_sources"].pfam is True
    assert cfg.annotation_targets["no_hit"].gff is True
    assert cfg.annotation_targets["no_hit"].pfam is False
    assert cfg.annotation_targets["negative_unmatched"].uniprot is False
    assert cfg.annotation_targets["negative_hit"].alphafold is False


def test_annotation_targets_can_enable_no_hit_pfam(tmp_path: Path) -> None:
    """Per-sheet annotation targets should override safe defaults."""
    data = valid_config_data()
    data["annotation_targets"] = {
        "no_hit": {
            "pfam": True,
        },
    }
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.annotation_targets["no_hit"].gff is True
    assert cfg.annotation_targets["no_hit"].pfam is True


def test_annotation_targets_missing_subkeys_keep_defaults(tmp_path: Path) -> None:
    """Partial per-sheet settings should keep unspecified default values."""
    data = valid_config_data()
    data["annotation_targets"] = {
        "candidates": {
            "pfam": False,
        },
        "negative_hit": {
            "uniprot": True,
        },
    }
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.annotation_targets["candidates"].gff is True
    assert cfg.annotation_targets["candidates"].pfam is False
    assert cfg.annotation_targets["candidates"].uniprot is True
    assert cfg.annotation_targets["negative_hit"].gff is True
    assert cfg.annotation_targets["negative_hit"].pfam is False
    assert cfg.annotation_targets["negative_hit"].uniprot is True


def test_invalid_annotation_target_boolean_raises_config_error(
    tmp_path: Path,
) -> None:
    """annotation_targets values should be true or false."""
    data = valid_config_data()
    data["annotation_targets"] = {
        "no_hit": {
            "pfam": "yes",
        },
    }
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="annotation_targets.no_hit.pfam"):
        load_config(config_path, initialize=False)


def test_missing_optional_gff_does_not_fail_config_validation(tmp_path: Path) -> None:
    """Optional GFF path should not be required."""
    data = valid_config_data()
    del data["paths"]["gff"]
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.paths.gff_file is None


def test_optional_gff_file_key_is_supported(tmp_path: Path) -> None:
    """paths.gff_file should be accepted as the preferred optional GFF key."""
    data = valid_config_data()
    del data["paths"]["gff"]
    data["paths"]["gff_file"] = "./data/input/custom.gff"
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.paths.gff_file == Path("./data/input/custom.gff")


def test_directory_input_mode_accepts_directory_paths(tmp_path: Path) -> None:
    """Directory mode should require source directories instead of FASTA files."""
    data = valid_config_data()
    data["input_mode"] = "directory"
    del data["paths"]["target_fasta"]
    del data["paths"]["positive_fasta"]
    del data["paths"]["negative_fasta"]
    data["paths"]["target_dir"] = "./data/databases/target"
    data["paths"]["positive_dir"] = "./data/databases/positive"
    data["paths"]["negative_dir"] = "./data/databases/negative"
    config_path = write_config(tmp_path, data)

    cfg = load_config(config_path, initialize=False)

    assert cfg.input_mode == "directory"
    assert cfg.paths.target_fasta is None
    assert cfg.paths.target_dir == Path("./data/databases/target")


def test_invalid_input_mode_raises_config_error(tmp_path: Path) -> None:
    """input_mode should be either file or directory."""
    data = valid_config_data()
    data["input_mode"] = "folders"
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="input_mode"):
        load_config(config_path, initialize=False)


def test_directory_mode_missing_directory_key_raises_config_error(
    tmp_path: Path,
) -> None:
    """Directory mode should clearly report missing directory paths."""
    data = valid_config_data()
    data["input_mode"] = "directory"
    data["paths"]["target_dir"] = "./data/databases/target"
    data["paths"]["positive_dir"] = "./data/databases/positive"
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="paths.negative_dir"):
        load_config(config_path, initialize=False)


def test_missing_paths_key_raises_config_error(tmp_path: Path) -> None:
    """Missing required path keys should raise ConfigError."""
    data = valid_config_data()
    del data["paths"]["target_fasta"]
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="paths.target_fasta"):
        load_config(config_path, initialize=False)


def test_empty_path_value_raises_config_error(tmp_path: Path) -> None:
    """Empty required path values should raise ConfigError."""
    data = valid_config_data()
    data["paths"]["positive_fasta"] = ""
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="paths.positive_fasta"):
        load_config(config_path, initialize=False)


@pytest.mark.parametrize("bad_evalue", [0, -1, "not-a-number"])
def test_invalid_blast_evalue_raises_config_error(
    tmp_path: Path,
    bad_evalue: object,
) -> None:
    """blast.evalue must be a positive number."""
    data = valid_config_data()
    data["blast"]["evalue"] = bad_evalue
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="blast.evalue"):
        load_config(config_path, initialize=False)


@pytest.mark.parametrize("bad_max_target_seqs", [0, -1, "ten", True])
def test_invalid_blast_max_target_seqs_raises_config_error(
    tmp_path: Path,
    bad_max_target_seqs: object,
) -> None:
    """blast.max_target_seqs must be a positive integer."""
    data = valid_config_data()
    data["blast"]["max_target_seqs"] = bad_max_target_seqs
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="blast.max_target_seqs"):
        load_config(config_path, initialize=False)


@pytest.mark.parametrize("bad_threads", [0, -1, "many", False])
def test_invalid_blast_threads_raises_config_error(
    tmp_path: Path,
    bad_threads: object,
) -> None:
    """blast.threads must be auto or a positive integer."""
    data = valid_config_data()
    data["blast"]["threads"] = bad_threads
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="blast.threads"):
        load_config(config_path, initialize=False)


def test_invalid_annotation_boolean_raises_config_error(tmp_path: Path) -> None:
    """Annotation enable flags must be booleans."""
    data = valid_config_data()
    data["annotation"]["enable_cdd"] = "yes"
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="annotation.enable_cdd"):
        load_config(config_path, initialize=False)


def test_invalid_cache_boolean_raises_config_error(tmp_path: Path) -> None:
    """Cache flags must be booleans."""
    data = valid_config_data()
    data["cache"]["enabled"] = "true"
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="cache.enabled"):
        load_config(config_path, initialize=False)


def test_invalid_logging_level_raises_config_error(tmp_path: Path) -> None:
    """logging.level must be one of the supported log levels."""
    data = valid_config_data()
    data["logging"]["level"] = "VERBOSE"
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="logging.level"):
        load_config(config_path, initialize=False)


def test_invalid_logging_save_log_raises_config_error(tmp_path: Path) -> None:
    """logging.save_log must be a boolean."""
    data = deepcopy(valid_config_data())
    data["logging"]["save_log"] = "yes"
    config_path = write_config(tmp_path, data)

    with pytest.raises(ConfigError, match="logging.save_log"):
        load_config(config_path, initialize=False)
