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
    target_fasta: Path
    positive_fasta: Path
    negative_fasta: Path
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

    paths: PathConfig
    blast: BlastConfig
    annotation: AnnotationConfig
    cache: CacheConfig
    score: ScoreConfig
    logging: LoggingConfig


# ==========================================================
# Config Loader
# ==========================================================

CONFIG_FILE = Path(__file__).parent / "config.yaml"


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
        target_fasta=Path(raw["paths"]["target_fasta"]),
        positive_fasta=Path(raw["paths"]["positive_fasta"]),
        negative_fasta=Path(raw["paths"]["negative_fasta"]),
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

    cache = CacheConfig(**raw["cache"])

    score = ScoreConfig(**raw["score"])

    logging = LoggingConfig(**raw["logging"])

    cfg = Config(
        project_name=raw["project"]["name"],
        version=raw["project"]["version"],
        paths=paths,
        blast=blast,
        annotation=annotation,
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

    _validate_paths_section(raw)
    _validate_blast_section(raw)
    _validate_annotation_section(raw)
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
    required_keys = (
        "target_fasta",
        "positive_fasta",
        "negative_fasta",
        "output_excel",
        "cache_dir",
        "log_dir",
    )

    for key in required_keys:
        value = _require_key(paths, "paths", key)
        if not isinstance(value, str) or not value.strip():
            raise ConfigError(
                f"config.yaml value 'paths.{key}' must be a non-empty string."
            )

    for key in ("gff", "gff_file"):
        value = paths.get(key)
        if value is not None and not isinstance(value, str):
            raise ConfigError(
                f"config.yaml value 'paths.{key}' must be a string when provided."
            )


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
