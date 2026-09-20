"""Run provenance: which code and which effective configuration produced a run.

Two artifacts come out of this module (see claude/run_provenance_design.md):

* a short fingerprint (:class:`RunProvenance`) that ``output/excel.py`` and
  ``output/word_report.py`` print in the workbook Index sheet and the report
  title page -- git commit, dirty flag, config hash, reference-genome check
  result -- so a reader can tell at a glance whether two outputs came from the
  same setup;
* a sidecar YAML next to the Excel workbook (``<stem>.run_provenance.yaml``)
  holding the full effective ``Config`` plus the effective scoring parameters
  (built-in code defaults included), so a temporary, uncommitted edit of
  ``config.yaml`` can still be reconstructed from the saved outputs.

Everything here is best-effort: a missing ``git`` binary, a source tree that is
not a git checkout, or an unreadable YAML file never stops a run; the affected
field just becomes ``None`` (or a placeholder in the config hash).
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any

import yaml

from core.constants import APP_VERSION

SIDECAR_SUFFIX = ".run_provenance.yaml"
CONFIG_HASH_LENGTH = 8
GIT_TIMEOUT_SECONDS = 10

#: Shown in the Excel Index / Word title page when the startup reference-genome
#: check (core/reference_genomes.py) reported findings.
REFERENCE_GENOME_WARNING_TEXT = "WARNING - see run log for details (core/reference_genomes.py)"


@dataclass(frozen=True)
class RunProvenance:
    app_version: str
    git_commit: str | None  # 7-character short SHA; None when unavailable
    git_dirty: bool | None  # uncommitted changes to tracked files; None when unavailable
    config_hash: str  # first CONFIG_HASH_LENGTH hex digits of a SHA-256, see compute_config_hash
    reference_genome_check_passed: bool | None  # None = not run / could not be verified
    generated_at: datetime


def provenance_sidecar_filename(output_excel: str | Path) -> str:
    """Sidecar filename for a run whose Excel workbook is ``output_excel``."""
    return f"{Path(output_excel).stem}{SIDECAR_SUFFIX}"


def provenance_sidecar_path(output_excel: str | Path) -> Path:
    """Sidecar path: same directory as the Excel workbook, named after its stem."""
    output_path = Path(output_excel)
    return output_path.with_name(provenance_sidecar_filename(output_path))


def code_version_text(provenance: RunProvenance) -> str:
    """e.g. ``5.0 (git 1a2b3c4+dirty)``; ``unknown`` when git provenance is unavailable."""
    dirty = "+dirty" if provenance.git_dirty else ""
    return f"{provenance.app_version} (git {provenance.git_commit or 'unknown'}{dirty})"


# ---------------------------------------------------------------------------
# git
# ---------------------------------------------------------------------------


def _run_git(args: list[str], cwd: Path) -> str | None:
    """Run ``git`` and return stripped stdout, or None on any failure."""
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=GIT_TIMEOUT_SECONDS,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout.strip()


def _git_state(repo_root: Path) -> tuple[str | None, bool | None]:
    """Return ``(short_sha, dirty)`` for ``repo_root``, or ``(None, None)``.

    ``repo_root`` must itself be the top level of the git checkout: an
    unpacked release zip that happens to sit inside some other repository
    would otherwise report that repository's commit.
    """
    toplevel = _run_git(["rev-parse", "--show-toplevel"], repo_root)
    if toplevel is None or Path(toplevel).resolve() != repo_root.resolve():
        return None, None
    commit = _run_git(["rev-parse", "--short=7", "HEAD"], repo_root)
    if not commit:
        return None, None
    status = _run_git(["status", "--porcelain", "--untracked-files=no"], repo_root)
    return commit, (None if status is None else bool(status))


# ---------------------------------------------------------------------------
# config hash
# ---------------------------------------------------------------------------


def effective_scoring_parameters(config: Any) -> dict[str, Any]:
    """Scoring parameters as they are *in effect*, including built-in code defaults.

    Hashing only the YAML files misses everything that lives in code: with no
    ``scoring_engine_config`` file the tier thresholds and category caps come
    from ``TierThresholds``/``DEFAULT_CATEGORY_CAPS``, and the per-component
    weights (``V2_COMPONENT_WEIGHTS``) are never configurable at all, so a
    change to any of them would leave the fingerprint unchanged. This resolves
    the engine config the way the pipeline does (file if configured, else
    defaults) and adds the module-level component weights.

    Best-effort like the rest of this module: a section that cannot be
    resolved becomes the placeholder ``"<unavailable>"`` instead of raising.
    """
    scoring = getattr(config, "interaction_scoring", None)
    try:
        from analysis.scoring_engine_config import load_scoring_engine_config

        engine: Any = _to_plain(load_scoring_engine_config(getattr(scoring, "scoring_engine_config", None)))
    except Exception:  # noqa: BLE001 - provenance must never stop a run
        engine = "<unavailable>"
    try:
        from analysis.interaction_scoring import V2_COMPONENT_WEIGHTS

        weights: Any = {str(name): float(weight) for name, weight in V2_COMPONENT_WEIGHTS.items()}
    except Exception:  # noqa: BLE001
        weights = "<unavailable>"
    return {"scoring_engine": engine, "v2_component_weights": weights}


def _config_hash_inputs(config: Any, config_path: Path | None) -> list[tuple[str, Path | bytes | None]]:
    """The files and effective parameters whose content defines scoring behavior, labelled.

    config.yaml itself, plus the external YAML files that change scores:
    the scoring-engine caps, the functional-complementarity rule set and the
    domain-family map; plus :func:`effective_scoring_parameters`, which covers
    the built-in code defaults those files fall back to. Data files
    (STRING/GEO caches, the PIH bundle, BLAST databases) are deliberately not
    hashed: they are data, not configuration, and reading them would cost real
    run time. Reference genomes are checked separately by
    core/reference_genomes.py.
    """
    scoring = getattr(config, "interaction_scoring", None)
    effective = json.dumps(effective_scoring_parameters(config), sort_keys=True, separators=(",", ":"))
    return [
        ("config.yaml", config_path),
        ("scoring_engine_config", getattr(scoring, "scoring_engine_config", None)),
        ("functional_complementarity_ruleset", getattr(scoring, "functional_complementarity_ruleset", None)),
        ("domain_family_map", getattr(scoring, "domain_family_map_path", None)),
        ("effective_scoring_parameters", effective.encode("utf-8")),
    ]


def compute_config_hash(inputs: list[tuple[str, Path | bytes | None]]) -> str:
    """SHA-256 (hex) over the labelled contents of ``inputs``.

    An input is a file path or, for values computed in memory, raw bytes. A
    ``None`` path (the built-in default is in use), a missing file and an
    unreadable file each contribute a distinct placeholder, so switching
    between them changes the hash. CRLF is normalized to LF first so the same
    file checked out on Windows and WSL hashes identically.
    """
    digest = hashlib.sha256()
    for label, source in inputs:
        digest.update(label.encode("utf-8") + b"\0")
        if source is None:
            digest.update(b"<default>")
        elif isinstance(source, bytes):
            digest.update(source)
        else:
            try:
                digest.update(Path(source).read_bytes().replace(b"\r\n", b"\n"))
            except FileNotFoundError:
                digest.update(b"<missing>")
            except OSError:
                digest.update(b"<unreadable>")
        digest.update(b"\0")
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# collect / write
# ---------------------------------------------------------------------------


def collect_run_provenance(
    config: Any,
    repo_root: Path,
    *,
    config_path: Path | None = None,
    reference_genome_check_passed: bool | None = None,
    generated_at: datetime | None = None,
) -> RunProvenance:
    """Gather the fingerprint for one run. Never raises.

    ``reference_genome_check_passed`` is the result of the startup check in
    core/reference_genomes.py (True: all genomes match the manifest, False:
    findings were logged, None: not run or not verifiable); it is recorded as
    is, not recomputed here.
    """
    commit, dirty = _git_state(Path(repo_root))
    return RunProvenance(
        app_version=APP_VERSION,
        git_commit=commit,
        git_dirty=dirty,
        config_hash=compute_config_hash(_config_hash_inputs(config, config_path))[:CONFIG_HASH_LENGTH],
        reference_genome_check_passed=reference_genome_check_passed,
        generated_at=generated_at or datetime.now(),
    )


def _to_plain(value: Any) -> Any:
    """Convert a config object tree into YAML-safe plain data (Path -> str, dataclass -> dict, ...)."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Enum):
        return _to_plain(value.value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, datetime):
        return value.isoformat(timespec="seconds")
    if dataclasses.is_dataclass(value) and not isinstance(value, type):
        return {field.name: _to_plain(getattr(value, field.name)) for field in dataclasses.fields(value)}
    if isinstance(value, dict):
        return {str(key): _to_plain(item) for key, item in value.items()}
    if isinstance(value, (set, frozenset)):
        return sorted((_to_plain(item) for item in value), key=str)
    if isinstance(value, (list, tuple)):
        return [_to_plain(item) for item in value]
    if hasattr(value, "__dict__"):
        return _to_plain(vars(value))
    return str(value)


def write_provenance_sidecar(
    provenance: RunProvenance,
    config: Any,
    output_path: str | Path,
    *,
    config_path: Path | None = None,
) -> Path:
    """Write ``<output_path stem>.run_provenance.yaml`` and return its path.

    ``output_path`` is the Excel workbook path. The file holds the
    fingerprint plus the full effective configuration (after defaults and
    presets were applied), meant to be read by a person and applied by hand;
    nothing reads it back automatically, so it can never overwrite a config.
    """
    sidecar_path = provenance_sidecar_path(output_path)
    payload: dict[str, Any] = {
        "run_provenance": _to_plain(provenance),
        "config_file": None if config_path is None else str(config_path),
        "effective_config": _to_plain(config),
        "effective_scoring_parameters": effective_scoring_parameters(config),
    }
    sidecar_path.parent.mkdir(parents=True, exist_ok=True)
    sidecar_path.write_text(
        yaml.safe_dump(payload, sort_keys=False, allow_unicode=True, default_flow_style=False),
        encoding="utf-8",
    )
    return sidecar_path
