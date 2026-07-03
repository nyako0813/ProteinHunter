"""
Protein Hunter v5
Configuration loader

Author: OpenAI + nyako
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import multiprocessing
import yaml

from core.exceptions import ConfigError


# ==========================================================
# Dataclasses
# ==========================================================

@dataclass
class PathConfig:
    target_fasta: Path | None
    positive_fasta: Path | None
    negative_fasta: Path | None
    target_dir: Path | None
    positive_dir: Path | None
    negative_dir: Path | None
    gff: Path | None
    gff_file: Path | None
    output_excel: Path
    cache_dir: Path
    log_dir: Path


@dataclass
class BlastConfig:
    evalue: float
    max_target_seqs: int
    threads: int


@dataclass
class AnnotationConfig:
    enable_cdd: bool
    enable_pfam: bool
    enable_alphafold: bool
    enable_uniprot: bool
    enable_gene_context: bool

    cdd_threads: int
    pfam_threads: int
    alphafold_threads: int
    pfam_evalue_threshold: float = 1e-5


@dataclass(frozen=True)
class AnnotationTargetConfig:
    """Annotation switches for one Excel classification sheet."""

    gff: bool
    pfam: bool
    uniprot: bool
    alphafold: bool


@dataclass(frozen=True)
class AnnotationTargetsConfig:
    """Annotation switches for all Excel classification sheets."""

    candidates: AnnotationTargetConfig
    positive_all_sources: AnnotationTargetConfig
    no_hit: AnnotationTargetConfig
    negative_unmatched: AnnotationTargetConfig
    negative_hit: AnnotationTargetConfig

    def items(self):
        """Return sheet target names and settings in workbook order."""
        return (
            ("candidates", self.candidates),
            ("positive_all_sources", self.positive_all_sources),
            ("no_hit", self.no_hit),
            ("negative_unmatched", self.negative_unmatched),
            ("negative_hit", self.negative_hit),
        )

    def get(self, key: str, default: object = None) -> object:
        """Dictionary-like access for existing call sites."""
        return getattr(self, key, default)

    def __getitem__(self, key: str) -> AnnotationTargetConfig:
        """Dictionary-like access for tests and backward compatibility."""
        return getattr(self, key)


@dataclass
class CacheConfig:
    enabled: bool
    overwrite: bool


@dataclass
class ScoreConfig:
    blast_weight: int
    domain_weight: int
    motif_weight: int
    gene_context_weight: int
    alphafold_weight: int


@dataclass
class LoggingConfig:
    level: str
    save_log: bool


@dataclass
class Config:

    project_name: str
    version: str
    input_mode: str

    paths: PathConfig
    blast: BlastConfig
    annotation: AnnotationConfig
    annotation_targets: AnnotationTargetsConfig
    cache: CacheConfig
    score: ScoreConfig
    logging: LoggingConfig


# ==========================================================
# Config Loader
# ==========================================================

CONFIG_FILE = Path(__file__).parent / "config.yaml"

ANNOTATION_TARGET_DEFAULTS: dict[str, AnnotationTargetConfig] = {
    "candidates": AnnotationTargetConfig(
        gff=True,
        pfam=True,
        uniprot=True,
        alphafold=True,
    ),
    "positive_all_sources": AnnotationTargetConfig(
        gff=True,
        pfam=True,
        uniprot=True,
        alphafold=True,
    ),
    "no_hit": AnnotationTargetConfig(
        gff=True,
        pfam=False,
        uniprot=False,
        alphafold=False,
    ),
    "negative_unmatched": AnnotationTargetConfig(
        gff=True,
        pfam=False,
        uniprot=False,
        alphafold=False,
    ),
    "negative_hit": AnnotationTargetConfig(
        gff=True,
        pfam=False,
        uniprot=False,
        alphafold=False,
    ),
}


def _auto_threads(value: object) -> int:

    if value == "auto":
        return max(1, multiprocessing.cpu_count() - 2)

    return int(value)


def load_config(config_file: str | Path = CONFIG_FILE, initialize: bool = True) -> Config:
    """Load and validate a ProteinHunter configuration file."""
    config_path = Path(config_file)

    if not config_path.exists():
        raise ConfigError(f"Configuration file was not found: {config_path}")

    with open(config_path, encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    validate_config(raw)

    optional_gff = _optional_path(raw["paths"].get("gff_file"))
    if optional_gff is None:
        optional_gff = _optional_path(raw["paths"].get("gff"))

    paths = PathConfig(
        target_fasta=_optional_path(raw["paths"].get("target_fasta")),
        positive_fasta=_optional_path(raw["paths"].get("positive_fasta")),
        negative_fasta=_optional_path(raw["paths"].get("negative_fasta")),
        target_dir=_optional_path(raw["paths"].get("target_dir")),
        positive_dir=_optional_path(raw["paths"].get("positive_dir")),
        negative_dir=_optional_path(raw["paths"].get("negative_dir")),
        gff=optional_gff,
        gff_file=optional_gff,
        output_excel=Path(raw["paths"]["output_excel"]),
        cache_dir=Path(raw["paths"]["cache_dir"]),
        log_dir=Path(raw["paths"]["log_dir"]),
    )

    blast = BlastConfig(
        evalue=float(raw["blast"]["evalue"]),
        max_target_seqs=int(raw["blast"]["max_target_seqs"]),
        threads=_auto_threads(raw["blast"]["threads"]),
    )

    annotation = AnnotationConfig(
        enable_cdd=raw["annotation"]["enable_cdd"],
        enable_pfam=raw["annotation"]["enable_pfam"],
        enable_alphafold=raw["annotation"]["enable_alphafold"],
        enable_uniprot=raw["annotation"]["enable_uniprot"],
        enable_gene_context=raw["annotation"].get("enable_gene_context", False),
        cdd_threads=int(raw["annotation"].get("cdd_threads", 1)),
        pfam_threads=int(raw["annotation"].get("pfam_threads", 1)),
        alphafold_threads=int(raw["annotation"].get("alphafold_threads", 1)),
        pfam_evalue_threshold=float(
            raw["annotation"].get("pfam_evalue_threshold", 1e-5)
        ),
    )
    annotation_targets = _load_annotation_targets(raw.get("annotation_targets"))

    cache = CacheConfig(**raw["cache"])

    score = ScoreConfig(**raw["score"])

    logging = LoggingConfig(**raw["logging"])

    cfg = Config(
        project_name=raw["project"]["name"],
        version=raw["project"]["version"],
        input_mode=raw.get("input_mode", "file"),
        paths=paths,
        blast=blast,
        annotation=annotation,
        annotation_targets=annotation_targets,
        cache=cache,
        score=score,
        logging=logging,
    )

    if initialize:
        initialize_directories(cfg)

    return cfg


def validate_config(raw: object) -> None:
    """Validate raw YAML config data and raise beginner-friendly errors."""
    if not isinstance(raw, dict):
        raise ConfigError("config.yaml must contain a YAML mapping at the top level.")

    _validate_input_mode(raw)
    _validate_paths_section(raw)
    _validate_blast_section(raw)
    _validate_annotation_section(raw)
    _validate_annotation_targets_section(raw)
    _validate_cache_section(raw)
    _validate_logging_section(raw)


def _section(raw: dict[object, object], name: str) -> dict[object, object]:
    """Return a required config section."""
    section = raw.get(name)
    if not isinstance(section, dict):
        raise ConfigError(f"config.yaml is missing the '{name}' section.")

    return section


def _require_key(section: dict[object, object], section_name: str, key: str) -> object:
    """Return a required config value."""
    if key not in section:
        raise ConfigError(f"config.yaml is missing '{section_name}.{key}'.")

    return section[key]


def _validate_paths_section(raw: dict[object, object]) -> None:
    """Validate required path settings."""
    paths = _section(raw, "paths")
    input_mode = raw.get("input_mode", "file")
    required_keys = (
        "output_excel",
        "cache_dir",
        "log_dir",
    )
    if input_mode == "file":
        required_keys = (
            "target_fasta",
            "positive_fasta",
            "negative_fasta",
            *required_keys,
        )
    else:
        required_keys = (
            "target_dir",
            "positive_dir",
            "negative_dir",
            *required_keys,
        )

    for key in required_keys:
        value = _require_key(paths, "paths", key)
        if not isinstance(value, str) or not value.strip():
            raise ConfigError(
                f"config.yaml value 'paths.{key}' must be a non-empty string."
            )

    optional_path_keys = (
        "target_fasta",
        "positive_fasta",
        "negative_fasta",
        "target_dir",
        "positive_dir",
        "negative_dir",
        "gff",
        "gff_file",
    )
    for key in optional_path_keys:
        value = paths.get(key)
        if value is not None and not isinstance(value, str):
            raise ConfigError(
                f"config.yaml value 'paths.{key}' must be a string when provided."
            )


def _validate_input_mode(raw: dict[object, object]) -> None:
    """Validate optional input mode setting."""
    input_mode = raw.get("input_mode", "file")
    if input_mode not in {"file", "directory"}:
        raise ConfigError("config.yaml value 'input_mode' must be 'file' or 'directory'.")


def _validate_blast_section(raw: dict[object, object]) -> None:
    """Validate BLAST settings."""
    blast = _section(raw, "blast")
    evalue = _require_key(blast, "blast", "evalue")
    max_target_seqs = _require_key(blast, "blast", "max_target_seqs")
    threads = _require_key(blast, "blast", "threads")

    try:
        evalue_number = float(evalue)
    except (TypeError, ValueError) as exc:
        raise ConfigError("config.yaml value 'blast.evalue' must be a positive number.") from exc

    if evalue_number <= 0:
        raise ConfigError("config.yaml value 'blast.evalue' must be greater than 0.")

    if not _is_positive_int(max_target_seqs):
        raise ConfigError(
            "config.yaml value 'blast.max_target_seqs' must be a positive integer."
        )

    if threads != "auto" and not _is_positive_int(threads):
        raise ConfigError(
            "config.yaml value 'blast.threads' must be 'auto' or a positive integer."
        )


def _validate_annotation_section(raw: dict[object, object]) -> None:
    """Validate annotation toggles."""
    annotation = _section(raw, "annotation")
    for key in ("enable_cdd", "enable_pfam", "enable_alphafold", "enable_uniprot"):
        value = _require_key(annotation, "annotation", key)
        if not isinstance(value, bool):
            raise ConfigError(
                f"config.yaml value 'annotation.{key}' must be true or false."
            )

    threshold = annotation.get("pfam_evalue_threshold", 1e-5)
    try:
        threshold_number = float(threshold)
    except (TypeError, ValueError) as exc:
        raise ConfigError(
            "config.yaml value 'annotation.pfam_evalue_threshold' must be a positive number."
        ) from exc

    if threshold_number <= 0:
        raise ConfigError(
            "config.yaml value 'annotation.pfam_evalue_threshold' must be greater than 0."
        )


def _validate_annotation_targets_section(raw: dict[object, object]) -> None:
    """Validate optional per-sheet annotation target switches."""
    targets = raw.get("annotation_targets")
    if targets is None:
        return
    if not isinstance(targets, dict):
        raise ConfigError("config.yaml value 'annotation_targets' must be a mapping.")

    valid_sheets = set(ANNOTATION_TARGET_DEFAULTS)
    valid_annotations = {"gff", "pfam", "uniprot", "alphafold"}
    for sheet_name, sheet_targets in targets.items():
        if not isinstance(sheet_name, str) or not sheet_name.strip():
            raise ConfigError(
                "config.yaml annotation_targets sheet names must be non-empty strings."
            )
        sheet_key = sheet_name.strip().lower()
        if sheet_key not in valid_sheets:
            raise ConfigError(
                "config.yaml value "
                f"'annotation_targets.{sheet_name}' is not supported. "
                "Use candidates, positive_all_sources, no_hit, "
                "negative_unmatched, or negative_hit."
            )
        if not isinstance(sheet_targets, dict):
            raise ConfigError(
                f"config.yaml value 'annotation_targets.{sheet_name}' must be a mapping."
            )
        for annotation_name, enabled in sheet_targets.items():
            if annotation_name not in valid_annotations:
                raise ConfigError(
                    "config.yaml value "
                    f"'annotation_targets.{sheet_name}.{annotation_name}' is not supported. "
                    "Use gff, pfam, uniprot, or alphafold."
                )
            if not isinstance(enabled, bool):
                raise ConfigError(
                    "config.yaml value "
                    f"'annotation_targets.{sheet_name}.{annotation_name}' "
                    "must be true or false."
                )


def _load_annotation_targets(
    raw_targets: object,
) -> AnnotationTargetsConfig:
    """Load per-sheet annotation target settings with safe defaults."""
    targets = dict(ANNOTATION_TARGET_DEFAULTS)
    if not isinstance(raw_targets, dict):
        return AnnotationTargetsConfig(**targets)

    for sheet_name, raw_sheet_targets in raw_targets.items():
        if not isinstance(sheet_name, str) or not isinstance(raw_sheet_targets, dict):
            continue

        key = sheet_name.strip().lower()
        if key not in targets:
            continue
        default = targets.get(
            key,
            AnnotationTargetConfig(
                gff=False,
                pfam=False,
                uniprot=False,
                alphafold=False,
            ),
        )
        targets[key] = AnnotationTargetConfig(
            gff=bool(raw_sheet_targets.get("gff", default.gff)),
            pfam=bool(raw_sheet_targets.get("pfam", default.pfam)),
            uniprot=bool(raw_sheet_targets.get("uniprot", default.uniprot)),
            alphafold=bool(raw_sheet_targets.get("alphafold", default.alphafold)),
        )

    return AnnotationTargetsConfig(**targets)


def _validate_cache_section(raw: dict[object, object]) -> None:
    """Validate cache settings."""
    cache = _section(raw, "cache")
    for key in ("enabled", "overwrite"):
        value = _require_key(cache, "cache", key)
        if not isinstance(value, bool):
            raise ConfigError(f"config.yaml value 'cache.{key}' must be true or false.")


def _validate_logging_section(raw: dict[object, object]) -> None:
    """Validate logging settings."""
    logging_section = _section(raw, "logging")
    level = _require_key(logging_section, "logging", "level")
    save_log = _require_key(logging_section, "logging", "save_log")
    valid_levels = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}

    if not isinstance(level, str) or level.upper() not in valid_levels:
        raise ConfigError(
            "config.yaml value 'logging.level' must be one of: "
            "DEBUG, INFO, WARNING, ERROR, CRITICAL."
        )

    if not isinstance(save_log, bool):
        raise ConfigError("config.yaml value 'logging.save_log' must be true or false.")


def _is_positive_int(value: object) -> bool:
    """Return True when value is an integer greater than zero."""
    if isinstance(value, bool):
        return False

    try:
        number = int(value)
    except (TypeError, ValueError):
        return False

    return number > 0 and str(value).strip() == str(number)


def _optional_path(value: object) -> Path | None:
    """Return a Path for an optional non-empty string value."""
    if not isinstance(value, str) or not value.strip():
        return None

    return Path(value)


# ==========================================================
# Directory initialization
# ==========================================================

def initialize_directories(cfg: Config):

    cfg.paths.cache_dir.mkdir(parents=True, exist_ok=True)

    cfg.paths.log_dir.mkdir(parents=True, exist_ok=True)

    cfg.paths.output_excel.parent.mkdir(
        parents=True,
        exist_ok=True,
    )


# ==========================================================
# Global instance
# ==========================================================

CONFIG = load_config()
