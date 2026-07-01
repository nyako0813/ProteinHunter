"""Tests for main.py command-line config handling."""

from __future__ import annotations

from main import build_arg_parser


def test_main_config_arg_defaults_to_config_yaml() -> None:
    """No --config argument should keep the production config default."""
    args = build_arg_parser().parse_args([])

    assert args.config == "config.yaml"


def test_main_config_arg_accepts_custom_config() -> None:
    """--config should allow a custom YAML config path."""
    args = build_arg_parser().parse_args(["--config", "config.demo.yaml"])

    assert args.config == "config.demo.yaml"
