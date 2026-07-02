"""Tests for main.py command-line config handling."""

from __future__ import annotations

from pathlib import Path

import yaml

from main import build_arg_parser


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
