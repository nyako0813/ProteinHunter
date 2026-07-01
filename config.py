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


# ==========================================================
# Dataclasses
# ==========================================================

@dataclass
class PathConfig:
    target_fasta: Path
    positive_fasta: Path
    negative_fasta: Path
    gff: Path
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


def _auto_threads(value):

    if value == "auto":
        return max(1, multiprocessing.cpu_count() - 2)

    return int(value)


def load_config() -> Config:

    if not CONFIG_FILE.exists():
        raise FileNotFoundError(CONFIG_FILE)

    with open(CONFIG_FILE, encoding="utf-8") as f:
        raw = yaml.safe_load(f)

    paths = PathConfig(
        target_fasta=Path(raw["paths"]["target_fasta"]),
        positive_fasta=Path(raw["paths"]["positive_fasta"]),
        negative_fasta=Path(raw["paths"]["negative_fasta"]),
        gff=Path(raw["paths"]["gff"]),
        output_excel=Path(raw["paths"]["output_excel"]),
        cache_dir=Path(raw["paths"]["cache_dir"]),
        log_dir=Path(raw["paths"]["log_dir"]),
    )

    blast = BlastConfig(
        evalue=float(raw["blast"]["evalue"]),
        max_target_seqs=int(raw["blast"]["max_target_seqs"]),
        threads=_auto_threads(raw["blast"]["threads"]),
    )

    annotation = AnnotationConfig(**raw["annotation"])

    cache = CacheConfig(**raw["cache"])

    score = ScoreConfig(**raw["score"])

    logging = LoggingConfig(**raw["logging"])

    return Config(
        project_name=raw["project"]["name"],
        version=raw["project"]["version"],
        paths=paths,
        blast=blast,
        annotation=annotation,
        cache=cache,
        score=score,
        logging=logging,
    )


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

initialize_directories(CONFIG)