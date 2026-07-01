# AGENTS.md

## Project Overview

ProteinHunter_v5 is a Python 3.12 / Linux-oriented research tool for identifying candidate uncharacterized proteins related to cyano-to-thioamide or sulfur-related enzymatic functions.

The project is intended for WSL Ubuntu 24.04 and Linux servers.

## Current Status

Implemented or partially implemented:

- `requirements.txt`
- `config.yaml`
- `config.py`
- `core/constants.py`
- `core/startup.py`
- `core/logger.py` partially implemented, but currently needs cleanup/rewrite

Known issue:

- `core/logger.py` contains invalid ordering:
  - `from __future__ import annotations` is not at the top
  - `latest`, `logfile`, and `return logger` are incorrectly placed outside `create_logger()`
  - It should be rewritten as a clean logger module

## Development Rules

Use:

- Python 3.12+
- Linux/WSL paths
- `pathlib.Path`
- type hints
- readable docstrings
- beginner-friendly error messages
- modular design

Avoid:

- hard-coded Windows paths
- editing system Python
- relying on Anaconda
- committing `.venv`, `.cache`, `logs`, or generated outputs

## Desired Architecture

```text
ProteinHunter_v5/
├── main.py
├── config.py
├── config.yaml
├── requirements.txt
├── AGENTS.md
├── core/
│   ├── constants.py
│   ├── startup.py
│   ├── logger.py
│   ├── models.py
│   ├── cache.py
│   └── exceptions.py
├── blast/
├── annotation/
├── analysis/
├── output/
├── tests/
├── data/
│   ├── input/
│   ├── output/
│   ├── databases/
│   └── temp/
├── logs/
└── .cache/