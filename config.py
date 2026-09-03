"""
Protein Hunter v5
Configuration loader

Author: OpenAI + nyako
"""

from __future__ import annotations

from dataclasses import dataclass, field
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
    candidates_relaxed: AnnotationTargetConfig
    positive_all_sources: AnnotationTargetConfig
    no_hit: AnnotationTargetConfig
    negative_unmatched: AnnotationTargetConfig
    negative_hit: AnnotationTargetConfig

    def items(self):
        """Return sheet target names and settings in workbook order."""
        return (
            ("candidates", self.candidates),
            ("candidates_relaxed", self.candidates_relaxed),
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


@dataclass(frozen=True)
class OrthologThresholdConfig:
    min_identity: float
    min_query_coverage: float
    max_evalue: float


@dataclass(frozen=True)
class OrthologFilterConfig:
    negative_exclusion_mode: str
    strong: OrthologThresholdConfig
    medium: OrthologThresholdConfig
    weak: OrthologThresholdConfig


@dataclass(frozen=True)
class InteractionQueryConfig:
    protein_id: str
    old_locus_tag: str
    sequence: str


@dataclass(frozen=True)
class InteractionScoringWeightsConfig:
    candidate_priority: float
    gene_neighborhood: float
    co_occurrence: float
    domain_complementarity: float
    alphafold_readiness: float
    # legacy_additive counterpart of v2's external_ppi_evidence category
    # cap (Phase 6a M4, analysis/scoring_engine_config.py::DEFAULT_CATEGORY_CAPS).
    # Same provisional value (15.0), not derived from a fit.
    external_ppi: float = 15.0


@dataclass(frozen=True)
class InteractionAlphaFoldConfig:
    enabled: bool
    max_pair_total_length: int


@dataclass(frozen=True)
class InteractionNeighborhoodConfig:
    enabled: bool
    max_distance_bp: int
    max_rows_per_query: int


@dataclass(frozen=True)
class InteractionEvidenceDetailConfig:
    """Scope for the optional Interaction_Evidence_Detail sheet.

    ``no_hit`` is excluded by default: it is typically the largest candidate
    source and, for scoring model v2, its per-pair ranks/scores are heavily
    tied (most no_hit candidates share one saturated score), so a
    component-level breakdown for it adds volume without much analytical
    value. Set ``include_no_hit: true`` to include it anyway.
    """

    include_no_hit: bool = False


@dataclass(frozen=True)
class InteractionScoringConfig:
    enabled: bool
    query_proteins: tuple[InteractionQueryConfig, ...]
    query_fasta: Path | None
    candidate_sources: dict[str, bool]
    max_candidates_per_query: int
    include_sequences_in_excel: bool
    scoring_weights: InteractionScoringWeightsConfig
    alphafold: InteractionAlphaFoldConfig
    neighborhood: InteractionNeighborhoodConfig
    # scoring model v2 (evidence-based, category-capped scoring). Defaults
    # keep every existing run on the original "legacy_additive" behavior.
    scoring_model: str = "legacy_additive"
    scoring_engine_config: Path | None = None
    functional_complementarity_ruleset: Path | None = None
    # Optional path to a ProteinInteractionHunter candidate_evidence_bundle
    # .jsonl file (scoring_model: v2_evidence_based only). PIH is never
    # imported; this is a plain-file bridge. See
    # analysis/pih_evidence_bridge.py.
    pih_evidence_bundle: Path | None = None
    evidence_detail_sheet: InteractionEvidenceDetailConfig = field(
        default_factory=InteractionEvidenceDetailConfig
    )
    # Which score determines candidate_rank / sheet ordering within each
    # Interaction_* source (design spec section 22, Phase 5 M5).
    # "interaction_priority_score" (default) preserves the exact pre-Phase-5
    # ranking behavior. "interaction_score" ranks by query-specific evidence
    # only (see analysis/interaction_scoring.py::INTERACTION_SCORE_COMPONENT_NAMES);
    # interaction_priority_score/evidence_tier/priority_group columns keep
    # their original values and meaning either way -- only candidate_rank
    # and row order change.
    ranking_metric: str = "interaction_priority_score"
    # NCBI taxonomy id for STRING (string-db.org) PPI evidence
    # (scoring_model: v2_evidence_based only; Phase 6a, see
    # claude/phase6_external_evidence_design.md). Must be STRING's own
    # taxid for the exact strain, which is not always the species-level
    # NCBI taxid -- e.g. M. acetivorans's species taxid 2214 returns
    # nothing from STRING; the strain-level taxid 188937 ("... C2A") is
    # what actually has data. Leave unset (None) to disable STRING
    # evidence entirely -- every pair is then MISSING for the
    # external_ppi_evidence category and genomic_context's
    # string_neighborhood component.
    string_ppi_ncbi_taxon_id: int | None = None
    # Public GEO coexpression evidence (scoring_model: v2_evidence_based
    # only; Phase 6b, see claude/phase6b_coexpression_design.md). Enables
    # both GSE77738 and GSE64349 bridges together (analysis/coexpression_bridge.py)
    # -- the datasets/accessions themselves are not configurable, unlike
    # STRING's taxon id, so a single on/off flag is enough. False by
    # default: enabling it downloads two GEO supplementary files (a few MB
    # each) on first use.
    geo_coexpression_enabled: bool = False


VALID_INTERACTION_SCORING_MODELS: tuple[str, ...] = ("legacy_additive", "v2_evidence_based")

VALID_INTERACTION_RANKING_METRICS: tuple[str, ...] = (
    "interaction_priority_score",
    "interaction_score",
)


@dataclass
class Config:

    project_name: str
    version: str
    input_mode: str

    paths: PathConfig
    blast: BlastConfig
    annotation: AnnotationConfig
    annotation_targets: AnnotationTargetsConfig
    ortholog_filter: OrthologFilterConfig
    interaction_scoring: InteractionScoringConfig
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
    "candidates_relaxed": AnnotationTargetConfig(
        gff=True,
        pfam=True,
        uniprot=False,
        alphafold=False,
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

ORTHOLOG_FILTER_DEFAULT = OrthologFilterConfig(
    negative_exclusion_mode="any_hit",
    strong=OrthologThresholdConfig(
        min_identity=40.0,
        min_query_coverage=70.0,
        max_evalue=1e-5,
    ),
    medium=OrthologThresholdConfig(
        min_identity=30.0,
        min_query_coverage=70.0,
        max_evalue=1e-5,
    ),
    weak=OrthologThresholdConfig(
        min_identity=25.0,
        min_query_coverage=50.0,
        max_evalue=1e-3,
    ),
)

INTERACTION_CANDIDATE_SOURCE_DEFAULTS: dict[str, bool] = {
    "candidates": True,
    "candidates_relaxed": True,
    "positive_all_sources": False,
    "negative_unmatched": False,
    "no_hit": True,
    "negative_hit": False,
    "negative_strong_hit": False,
    "negative_medium_hit": False,
    "negative_weak_hit": False,
}

INTERACTION_SCORING_WEIGHTS_DEFAULT = InteractionScoringWeightsConfig(
    candidate_priority=30.0,
    gene_neighborhood=25.0,
    co_occurrence=20.0,
    domain_complementarity=15.0,
    alphafold_readiness=10.0,
)

INTERACTION_ALPHAFOLD_DEFAULT = InteractionAlphaFoldConfig(
    enabled=False,
    max_pair_total_length=2500,
)

INTERACTION_NEIGHBORHOOD_DEFAULT = InteractionNeighborhoodConfig(
    enabled=True,
    max_distance_bp=100000,
    max_rows_per_query=200,
)

INTERACTION_EVIDENCE_DETAIL_DEFAULT = InteractionEvidenceDetailConfig(
    include_no_hit=False,
)

INTERACTION_SCORING_DEFAULT = InteractionScoringConfig(
    enabled=False,
    query_proteins=(),
    query_fasta=None,
    candidate_sources=dict(INTERACTION_CANDIDATE_SOURCE_DEFAULTS),
    max_candidates_per_query=200,
    include_sequences_in_excel=False,
    scoring_weights=INTERACTION_SCORING_WEIGHTS_DEFAULT,
    alphafold=INTERACTION_ALPHAFOLD_DEFAULT,
    neighborhood=INTERACTION_NEIGHBORHOOD_DEFAULT,
    evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
    ranking_metric="interaction_priority_score",
)


def _default_interaction_scoring() -> InteractionScoringConfig:
    """Return a fresh disabled interaction scoring config."""
    return InteractionScoringConfig(
        enabled=False,
        query_proteins=(),
        query_fasta=None,
        candidate_sources=dict(INTERACTION_CANDIDATE_SOURCE_DEFAULTS),
        max_candidates_per_query=200,
        include_sequences_in_excel=False,
        scoring_weights=INTERACTION_SCORING_WEIGHTS_DEFAULT,
        alphafold=INTERACTION_ALPHAFOLD_DEFAULT,
        neighborhood=INTERACTION_NEIGHBORHOOD_DEFAULT,
        evidence_detail_sheet=INTERACTION_EVIDENCE_DETAIL_DEFAULT,
        ranking_metric="interaction_priority_score",
    )


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
    ortholog_filter = _load_ortholog_filter(raw.get("ortholog_filter"))
    interaction_scoring = _load_interaction_scoring(raw.get("interaction_scoring"))

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
        ortholog_filter=ortholog_filter,
        interaction_scoring=interaction_scoring,
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
    _validate_ortholog_filter_section(raw)
    _validate_interaction_scoring_section(raw)
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
                "Use candidates, candidates_relaxed, positive_all_sources, no_hit, "
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


def _validate_ortholog_filter_section(raw: dict[object, object]) -> None:
    """Validate optional ortholog-aware negative-hit settings."""
    section = raw.get("ortholog_filter")
    if section is None:
        return
    if not isinstance(section, dict):
        raise ConfigError("config.yaml value 'ortholog_filter' must be a mapping.")

    mode = section.get("negative_exclusion_mode", "any_hit")
    if mode not in {"any_hit", "strong_only", "strong_or_medium", "none"}:
        raise ConfigError(
            "config.yaml value 'ortholog_filter.negative_exclusion_mode' must be "
            "any_hit, strong_only, strong_or_medium, or none."
        )

    for level in ("strong", "medium", "weak"):
        raw_threshold = section.get(level, {})
        if raw_threshold is None:
            raw_threshold = {}
        if not isinstance(raw_threshold, dict):
            raise ConfigError(
                f"config.yaml value 'ortholog_filter.{level}' must be a mapping."
            )
        for key in ("min_identity", "min_query_coverage", "max_evalue"):
            if key not in raw_threshold:
                continue
            try:
                value = float(raw_threshold[key])
            except (TypeError, ValueError) as exc:
                raise ConfigError(
                    f"config.yaml value 'ortholog_filter.{level}.{key}' "
                    "must be a number."
                ) from exc
            if value < 0:
                raise ConfigError(
                    f"config.yaml value 'ortholog_filter.{level}.{key}' "
                    "must be greater than or equal to 0."
                )


def _load_ortholog_filter(raw_filter: object) -> OrthologFilterConfig:
    """Load ortholog-aware negative-hit settings with backward-compatible defaults."""
    if not isinstance(raw_filter, dict):
        return ORTHOLOG_FILTER_DEFAULT

    def threshold(
        name: str,
        default: OrthologThresholdConfig,
    ) -> OrthologThresholdConfig:
        raw_threshold = raw_filter.get(name, {})
        if not isinstance(raw_threshold, dict):
            raw_threshold = {}

        return OrthologThresholdConfig(
            min_identity=float(
                raw_threshold.get("min_identity", default.min_identity)
            ),
            min_query_coverage=float(
                raw_threshold.get(
                    "min_query_coverage",
                    default.min_query_coverage,
                )
            ),
            max_evalue=float(raw_threshold.get("max_evalue", default.max_evalue)),
        )

    return OrthologFilterConfig(
        negative_exclusion_mode=str(
            raw_filter.get(
                "negative_exclusion_mode",
                ORTHOLOG_FILTER_DEFAULT.negative_exclusion_mode,
            )
        ),
        strong=threshold("strong", ORTHOLOG_FILTER_DEFAULT.strong),
        medium=threshold("medium", ORTHOLOG_FILTER_DEFAULT.medium),
        weak=threshold("weak", ORTHOLOG_FILTER_DEFAULT.weak),
    )


def _validate_interaction_scoring_section(raw: dict[object, object]) -> None:
    """Validate optional lightweight interaction scoring settings."""
    section = raw.get("interaction_scoring")
    if section is None:
        return
    if not isinstance(section, dict):
        raise ConfigError("config.yaml value 'interaction_scoring' must be a mapping.")

    enabled = section.get("enabled", False)
    if not isinstance(enabled, bool):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.enabled' must be true or false."
        )

    query_proteins = section.get("query_proteins", [])
    if query_proteins is None:
        query_proteins = []
    if not isinstance(query_proteins, list):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.query_proteins' must be a list."
        )
    for index, query in enumerate(query_proteins):
        if not isinstance(query, dict):
            raise ConfigError(
                "config.yaml value "
                f"'interaction_scoring.query_proteins[{index}]' must be a mapping."
            )
        for key in ("protein_id", "old_locus_tag", "sequence"):
            value = query.get(key, "")
            if value is not None and not isinstance(value, str):
                raise ConfigError(
                    "config.yaml value "
                    f"'interaction_scoring.query_proteins[{index}].{key}' "
                    "must be a string."
                )

    query_fasta = section.get("query_fasta", "")
    if query_fasta is not None and not isinstance(query_fasta, str):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.query_fasta' must be a string."
        )

    candidate_sources = section.get("candidate_sources", {})
    if candidate_sources is None:
        candidate_sources = {}
    if not isinstance(candidate_sources, dict):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.candidate_sources' "
            "must be a mapping."
        )
    valid_sources = set(INTERACTION_CANDIDATE_SOURCE_DEFAULTS)
    for source_name, enabled_source in candidate_sources.items():
        if source_name not in valid_sources:
            raise ConfigError(
                "config.yaml value "
                f"'interaction_scoring.candidate_sources.{source_name}' "
                "is not supported."
            )
        if not isinstance(enabled_source, bool):
            raise ConfigError(
                "config.yaml value "
                f"'interaction_scoring.candidate_sources.{source_name}' "
                "must be true or false."
            )

    max_candidates = section.get("max_candidates_per_query", 200)
    if not _is_positive_int(max_candidates):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.max_candidates_per_query' "
            "must be a positive integer."
        )

    include_sequences = section.get("include_sequences_in_excel", False)
    if not isinstance(include_sequences, bool):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.include_sequences_in_excel' "
            "must be true or false."
        )

    scoring_weights = section.get("scoring_weights", {})
    if scoring_weights is None:
        scoring_weights = {}
    if not isinstance(scoring_weights, dict):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.scoring_weights' "
            "must be a mapping."
        )
    for key in (
        "candidate_priority",
        "gene_neighborhood",
        "co_occurrence",
        "domain_complementarity",
        "alphafold_readiness",
        "external_ppi",
    ):
        if key not in scoring_weights:
            continue
        try:
            float(scoring_weights[key])
        except (TypeError, ValueError) as exc:
            raise ConfigError(
                "config.yaml value "
                f"'interaction_scoring.scoring_weights.{key}' must be a number."
            ) from exc

    alphafold = section.get("alphafold", {})
    if alphafold is None:
        alphafold = {}
    if not isinstance(alphafold, dict):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.alphafold' must be a mapping."
        )
    alphafold_enabled = alphafold.get("enabled", False)
    if not isinstance(alphafold_enabled, bool):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.alphafold.enabled' "
            "must be true or false."
        )
    max_pair_total_length = alphafold.get("max_pair_total_length", 2500)
    if not _is_positive_int(max_pair_total_length):
        raise ConfigError(
            "config.yaml value "
            "'interaction_scoring.alphafold.max_pair_total_length' "
            "must be a positive integer."
        )

    neighborhood = section.get("neighborhood", {})
    if neighborhood is None:
        neighborhood = {}
    if not isinstance(neighborhood, dict):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.neighborhood' must be a mapping."
        )
    neighborhood_enabled = neighborhood.get("enabled", True)
    if not isinstance(neighborhood_enabled, bool):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.neighborhood.enabled' "
            "must be true or false."
        )
    max_distance_bp = neighborhood.get("max_distance_bp", 100000)
    if not _is_positive_int(max_distance_bp):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.neighborhood.max_distance_bp' "
            "must be a positive integer."
        )
    max_rows_per_query = neighborhood.get("max_rows_per_query", 200)
    if not _is_positive_int(max_rows_per_query):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.neighborhood.max_rows_per_query' "
            "must be a positive integer."
        )

    scoring_model = section.get("scoring_model", "legacy_additive")
    if scoring_model not in VALID_INTERACTION_SCORING_MODELS:
        raise ConfigError(
            "config.yaml value 'interaction_scoring.scoring_model' must be one of "
            f"{VALID_INTERACTION_SCORING_MODELS}, got {scoring_model!r}."
        )

    scoring_engine_config = section.get("scoring_engine_config")
    if scoring_engine_config is not None and not isinstance(scoring_engine_config, str):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.scoring_engine_config' "
            "must be a string path."
        )

    functional_complementarity_ruleset = section.get("functional_complementarity_ruleset")
    if functional_complementarity_ruleset is not None and not isinstance(
        functional_complementarity_ruleset, str
    ):
        raise ConfigError(
            "config.yaml value "
            "'interaction_scoring.functional_complementarity_ruleset' "
            "must be a string path."
        )

    pih_evidence_bundle = section.get("pih_evidence_bundle")
    if pih_evidence_bundle is not None and not isinstance(pih_evidence_bundle, str):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.pih_evidence_bundle' "
            "must be a string path."
        )

    evidence_detail_sheet = section.get("evidence_detail_sheet", {})
    if evidence_detail_sheet is None:
        evidence_detail_sheet = {}
    if not isinstance(evidence_detail_sheet, dict):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.evidence_detail_sheet' "
            "must be a mapping."
        )
    include_no_hit = evidence_detail_sheet.get("include_no_hit", False)
    if not isinstance(include_no_hit, bool):
        raise ConfigError(
            "config.yaml value "
            "'interaction_scoring.evidence_detail_sheet.include_no_hit' "
            "must be true or false."
        )

    ranking_metric = section.get("ranking_metric", "interaction_priority_score")
    if ranking_metric not in VALID_INTERACTION_RANKING_METRICS:
        raise ConfigError(
            "config.yaml value 'interaction_scoring.ranking_metric' must be "
            f"one of {VALID_INTERACTION_RANKING_METRICS}, got {ranking_metric!r}."
        )

    string_ppi_ncbi_taxon_id = section.get("string_ppi_ncbi_taxon_id")
    if string_ppi_ncbi_taxon_id is not None and not _is_positive_int(string_ppi_ncbi_taxon_id):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.string_ppi_ncbi_taxon_id' "
            "must be a positive integer or left unset."
        )

    geo_coexpression_enabled = section.get("geo_coexpression_enabled", False)
    if not isinstance(geo_coexpression_enabled, bool):
        raise ConfigError(
            "config.yaml value 'interaction_scoring.geo_coexpression_enabled' "
            "must be true or false."
        )


def _load_interaction_scoring(raw_scoring: object) -> InteractionScoringConfig:
    """Load optional lightweight interaction scoring settings."""
    if not isinstance(raw_scoring, dict):
        return _default_interaction_scoring()

    raw_queries = raw_scoring.get("query_proteins", [])
    if not isinstance(raw_queries, list):
        raw_queries = []
    query_proteins: list[InteractionQueryConfig] = []
    for raw_query in raw_queries:
        if not isinstance(raw_query, dict):
            continue
        query_proteins.append(
            InteractionQueryConfig(
                protein_id=str(raw_query.get("protein_id") or ""),
                old_locus_tag=str(raw_query.get("old_locus_tag") or ""),
                sequence=str(raw_query.get("sequence") or ""),
            )
        )

    raw_candidate_sources = raw_scoring.get("candidate_sources", {})
    candidate_sources = dict(INTERACTION_CANDIDATE_SOURCE_DEFAULTS)
    if isinstance(raw_candidate_sources, dict):
        for source_name, enabled in raw_candidate_sources.items():
            if source_name in candidate_sources:
                candidate_sources[source_name] = bool(enabled)

    raw_weights = raw_scoring.get("scoring_weights", {})
    if not isinstance(raw_weights, dict):
        raw_weights = {}
    scoring_weights = InteractionScoringWeightsConfig(
        candidate_priority=float(
            raw_weights.get(
                "candidate_priority",
                INTERACTION_SCORING_WEIGHTS_DEFAULT.candidate_priority,
            )
        ),
        gene_neighborhood=float(
            raw_weights.get(
                "gene_neighborhood",
                INTERACTION_SCORING_WEIGHTS_DEFAULT.gene_neighborhood,
            )
        ),
        co_occurrence=float(
            raw_weights.get(
                "co_occurrence",
                INTERACTION_SCORING_WEIGHTS_DEFAULT.co_occurrence,
            )
        ),
        domain_complementarity=float(
            raw_weights.get(
                "domain_complementarity",
                INTERACTION_SCORING_WEIGHTS_DEFAULT.domain_complementarity,
            )
        ),
        alphafold_readiness=float(
            raw_weights.get(
                "alphafold_readiness",
                INTERACTION_SCORING_WEIGHTS_DEFAULT.alphafold_readiness,
            )
        ),
        external_ppi=float(
            raw_weights.get(
                "external_ppi",
                INTERACTION_SCORING_WEIGHTS_DEFAULT.external_ppi,
            )
        ),
    )

    raw_alphafold = raw_scoring.get("alphafold", {})
    if not isinstance(raw_alphafold, dict):
        raw_alphafold = {}
    alphafold = InteractionAlphaFoldConfig(
        enabled=bool(raw_alphafold.get("enabled", INTERACTION_ALPHAFOLD_DEFAULT.enabled)),
        max_pair_total_length=int(
            raw_alphafold.get(
                "max_pair_total_length",
                INTERACTION_ALPHAFOLD_DEFAULT.max_pair_total_length,
            )
        ),
    )

    raw_neighborhood = raw_scoring.get("neighborhood", {})
    if not isinstance(raw_neighborhood, dict):
        raw_neighborhood = {}
    neighborhood = InteractionNeighborhoodConfig(
        enabled=bool(
            raw_neighborhood.get(
                "enabled",
                INTERACTION_NEIGHBORHOOD_DEFAULT.enabled,
            )
        ),
        max_distance_bp=int(
            raw_neighborhood.get(
                "max_distance_bp",
                INTERACTION_NEIGHBORHOOD_DEFAULT.max_distance_bp,
            )
        ),
        max_rows_per_query=int(
            raw_neighborhood.get(
                "max_rows_per_query",
                INTERACTION_NEIGHBORHOOD_DEFAULT.max_rows_per_query,
            )
        ),
    )

    raw_evidence_detail_sheet = raw_scoring.get("evidence_detail_sheet", {})
    if not isinstance(raw_evidence_detail_sheet, dict):
        raw_evidence_detail_sheet = {}
    evidence_detail_sheet = InteractionEvidenceDetailConfig(
        include_no_hit=bool(
            raw_evidence_detail_sheet.get(
                "include_no_hit",
                INTERACTION_EVIDENCE_DETAIL_DEFAULT.include_no_hit,
            )
        ),
    )

    return InteractionScoringConfig(
        enabled=bool(raw_scoring.get("enabled", False)),
        query_proteins=tuple(query_proteins),
        query_fasta=_optional_path(raw_scoring.get("query_fasta")),
        candidate_sources=candidate_sources,
        max_candidates_per_query=int(raw_scoring.get("max_candidates_per_query", 200)),
        include_sequences_in_excel=bool(
            raw_scoring.get("include_sequences_in_excel", False)
        ),
        scoring_weights=scoring_weights,
        alphafold=alphafold,
        neighborhood=neighborhood,
        scoring_model=str(raw_scoring.get("scoring_model", "legacy_additive")),
        scoring_engine_config=_optional_path(raw_scoring.get("scoring_engine_config")),
        functional_complementarity_ruleset=_optional_path(
            raw_scoring.get("functional_complementarity_ruleset")
        ),
        pih_evidence_bundle=_optional_path(raw_scoring.get("pih_evidence_bundle")),
        evidence_detail_sheet=evidence_detail_sheet,
        ranking_metric=str(
            raw_scoring.get("ranking_metric", "interaction_priority_score")
        ),
        string_ppi_ncbi_taxon_id=(
            int(raw_scoring["string_ppi_ncbi_taxon_id"])
            if raw_scoring.get("string_ppi_ncbi_taxon_id") is not None
            else None
        ),
        geo_coexpression_enabled=bool(raw_scoring.get("geo_coexpression_enabled", False)),
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
