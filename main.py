"""
Protein Hunter v5
Main entry point.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from collections.abc import Sequence
from typing import Any

from core.startup import StartupChecker
from core.models import ProteinRecord


def _require_path(path: Path | None, config_key: str) -> Path:
    """Return a configured path or raise a friendly config error."""
    if path is None:
        from core.exceptions import ConfigError

        raise ConfigError(f"config.yaml is missing '{config_key}'.")

    return path


SHEET_TO_ANNOTATION_TARGET: dict[str, str] = {
    "Candidates": "candidates",
    "Candidates_relaxed": "candidates_relaxed",
    "Positive_all_sources": "positive_all_sources",
    "Negative_unmatched": "negative_unmatched",
    "No_hit": "no_hit",
    "Negative_hit": "negative_hit",
    "Negative_strong_hit": "negative_hit",
    "Negative_medium_hit": "negative_hit",
    "Negative_weak_hit": "negative_hit",
}


def _classification_record_sheets(
    blast_classification: Any,
) -> dict[str, dict[str, ProteinRecord]]:
    """Return Excel classification sheets that can receive annotations."""
    return {
        "Candidates": blast_classification.positive_only_records,
        "Candidates_relaxed": blast_classification.candidates_relaxed_records,
        "Positive_all_sources": blast_classification.positive_all_sources_records,
        "Negative_unmatched": blast_classification.negative_unmatched_records,
        "No_hit": blast_classification.no_hit_records,
        "Negative_hit": blast_classification.negative_hit_records,
        "Negative_strong_hit": blast_classification.negative_strong_hit_records,
        "Negative_medium_hit": blast_classification.negative_medium_hit_records,
        "Negative_weak_hit": blast_classification.negative_weak_hit_records,
    }


def _records_enabled_for_annotation(
    record_sheets: dict[str, dict[str, ProteinRecord]],
    annotation_targets: dict[str, Any],
    annotation_name: str,
) -> dict[str, ProteinRecord]:
    """Return de-duplicated records enabled for one annotation source."""
    selected: dict[str, ProteinRecord] = {}

    for sheet_name, records in record_sheets.items():
        target_key = SHEET_TO_ANNOTATION_TARGET.get(sheet_name)
        if target_key is None:
            continue

        target = annotation_targets.get(target_key)
        if target is None or not bool(getattr(target, annotation_name, False)):
            continue

        for protein_id, record in records.items():
            selected.setdefault(protein_id, record)

    return selected


def _records_for_annotation_step(
    blast_classification: Any,
    annotation_targets: Any,
    annotation_name: str,
) -> dict[str, ProteinRecord]:
    """Return de-duplicated classification records enabled for one annotation step."""
    return _records_enabled_for_annotation(
        _classification_record_sheets(blast_classification),
        annotation_targets,
        annotation_name,
    )


def _build_cdd_target_records(
    positive_only_records: dict[str, ProteinRecord],
    config: Any,
    blast_classification: Any,
) -> dict[str, ProteinRecord]:
    """Return the full CDD annotation target set for this run.

    Always includes ``positive_only_records`` ("Candidates"), the existing
    default scope, plus (when interaction_scoring is enabled) every record
    interaction_scoring can actually score or query -- see
    analysis/interaction_scoring.py::resolve_cdd_annotation_targets. This
    is a plain dict merge; it deliberately does not mutate
    ``positive_only_records`` itself, so the caller's own reference (used
    later for candidate scoring and the Excel "Candidates" sheet) keeps its
    original size regardless of how many extra records CDD ends up
    annotating.
    """
    from analysis.interaction_scoring import resolve_cdd_annotation_targets

    target_records = dict(positive_only_records)
    target_records.update(resolve_cdd_annotation_targets(config, blast_classification))
    return target_records


def build_arg_parser() -> argparse.ArgumentParser:
    """Build the command-line parser for Protein Hunter."""
    parser = argparse.ArgumentParser(
        description="Run the ProteinHunter_v5 analysis pipeline.",
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to the YAML configuration file. Defaults to config.yaml.",
    )
    parser.add_argument(
        "--check-only",
        action="store_true",
        help=(
            "Validate startup, config, and input FASTA files, then stop before "
            "BLAST and annotation."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> None:
    """Run Protein Hunter."""
    args = build_arg_parser().parse_args(argv)

    checker = StartupChecker()

    if not checker.run():
        raise SystemExit(1)

    from core.logger import logger

    try:
        from annotation.domain_annotator import (
            annotate_records_cdd,
            annotate_records_pfam,
        )
        from annotation.record_annotator import (
            annotate_records_alphafold,
            annotate_records_uniprot,
        )
        from annotation.gff import (
            annotate_records_with_gff_locus_tags,
            load_gff_locus_map,
        )
        from analysis.blast_pipeline import run_blast_classification_pipeline
        from analysis.input_summary import format_input_summary, summarize_input_fastas
        from analysis.interaction_scoring import run_interaction_scoring
        from analysis.scoring import get_sorted_records, score_records
        from config import load_config, redundant_negative_hit_sources
        from core.cache import JsonCache
        from core.fasta_sources import DirectoryFastaResult, prepare_directory_fasta
        from output.excel import write_classification_workbook

        config_path = Path(args.config)
        config = load_config(config_path)

        logger.info("Protein Hunter started")
        logger.success("Startup check passed")
        logger.info(f"Using config file: {config_path}")

        blast_work_dir = Path("data") / "temp" / "blast"
        directory_results: dict[str, DirectoryFastaResult] = {}
        source_counts: dict[str, int] = {}

        if config.input_mode == "directory":
            directory_results = {
                "target": prepare_directory_fasta(
                    _require_path(config.paths.target_dir, "paths.target_dir"),
                    "target",
                ),
                "positive": prepare_directory_fasta(
                    _require_path(config.paths.positive_dir, "paths.positive_dir"),
                    "positive",
                ),
                "negative": prepare_directory_fasta(
                    _require_path(config.paths.negative_dir, "paths.negative_dir"),
                    "negative",
                ),
            }
            target_fasta = directory_results["target"].combined_fasta
            positive_fasta = directory_results["positive"].combined_fasta
            negative_fasta = directory_results["negative"].combined_fasta
            source_counts = {
                "target_sources": len(directory_results["target"].source_labels),
                "positive_sources": len(directory_results["positive"].source_labels),
                "negative_sources": len(directory_results["negative"].source_labels),
            }
            positive_source_labels = directory_results["positive"].source_labels
        else:
            target_fasta = _require_path(config.paths.target_fasta, "paths.target_fasta")
            positive_fasta = _require_path(
                config.paths.positive_fasta,
                "paths.positive_fasta",
            )
            negative_fasta = _require_path(
                config.paths.negative_fasta,
                "paths.negative_fasta",
            )
            positive_source_labels = ("positive_fasta",)

        with logger.section("Configuration"):
            logger.info(f"Input mode: {config.input_mode}")
            if config.input_mode == "directory":
                logger.info(f"Target directory: {config.paths.target_dir}")
                logger.info(f"Positive directory: {config.paths.positive_dir}")
                logger.info(f"Negative directory: {config.paths.negative_dir}")
                for category, title in (
                    ("target", "Target sources"),
                    ("positive", "Positive sources"),
                    ("negative", "Negative sources"),
                ):
                    result = directory_results[category]
                    logger.info(f"{title}:")
                    for source_label in result.source_labels:
                        logger.info(f"- {source_label}")
                    if result.skipped_folders:
                        logger.warning(
                            f"Skipped {category} folders without protein.faa: "
                            f"{', '.join(result.skipped_folders)}"
                        )
                    if result.multiple_file_labels:
                        logger.warning(
                            f"Multiple protein.faa files found under {category} "
                            f"source folders: {', '.join(result.multiple_file_labels)}"
                        )
                    if result.duplicate_ids:
                        logger.warning(
                            f"Duplicate FASTA IDs in {category} sources: "
                            f"{', '.join(result.duplicate_ids)}"
                        )
                logger.info(f"Combined target FASTA: {target_fasta}")
                logger.info(f"Combined positive FASTA: {positive_fasta}")
                logger.info(f"Combined negative FASTA: {negative_fasta}")
            else:
                logger.info(f"Target FASTA: {target_fasta}")
                logger.info(f"Positive FASTA: {positive_fasta}")
                logger.info(f"Negative FASTA: {negative_fasta}")
            logger.info(f"Excel output: {config.paths.output_excel}")
            logger.info(f"BLAST work directory: {blast_work_dir}")
            logger.info(f"Cache directory: {config.paths.cache_dir}")
            if config.paths.gff_file is not None:
                logger.info(f"Optional GFF file: {config.paths.gff_file}")
            else:
                logger.info("Optional GFF file: not configured")
            logger.info(
                "Negative exclusion mode: "
                f"{config.ortholog_filter.negative_exclusion_mode}"
            )
            for level in ("strong", "medium", "weak"):
                threshold = getattr(config.ortholog_filter, level)
                logger.info(
                    f"{level} negative threshold: "
                    f"identity >= {threshold.min_identity}, "
                    f"query coverage >= {threshold.min_query_coverage}, "
                    f"e-value <= {threshold.max_evalue}"
                )
            logger.info(
                "Interaction scoring: "
                f"{'enabled' if config.interaction_scoring.enabled else 'disabled'}"
            )
            redundant_sources = redundant_negative_hit_sources(
                config.interaction_scoring.candidate_sources
            )
            if redundant_sources:
                logger.warning(
                    "interaction_scoring.candidate_sources.negative_hit is enabled "
                    f"together with {', '.join(redundant_sources)} -- these are "
                    "subsets of negative_hit, so matching candidates will be scored "
                    "and reported twice (once per bucket) with different "
                    "candidate_rank values for the same pair. Usually only one of "
                    "them should be enabled."
                )

            input_summary = summarize_input_fastas(
                target_fasta=target_fasta,
                positive_fasta=positive_fasta,
                negative_fasta=negative_fasta,
                source_counts=source_counts,
            )
            logger.info("Input FASTA summary:")
            for line in format_input_summary(input_summary):
                logger.info(line)

        if args.check_only:
            logger.summary()
            logger.success(
                "Check-only mode completed successfully. "
                "No BLAST or annotation was run."
            )
            return

        with logger.section("BLAST candidate search"):
            with logger.timer("BLAST candidate pipeline"):
                blast_classification = run_blast_classification_pipeline(
                    target_fasta=target_fasta,
                    positive_fasta=positive_fasta,
                    negative_fasta=negative_fasta,
                    work_dir=blast_work_dir,
                    evalue=config.blast.evalue,
                    max_target_seqs=config.blast.max_target_seqs,
                    threads=config.blast.threads,
                    positive_source_labels=positive_source_labels,
                    ortholog_filter=config.ortholog_filter,
                )
                records = blast_classification.positive_only_records

            logger.info(f"Total target proteins: {len(blast_classification.all_records)}")
            logger.info(f"BLAST positive-only candidates: {len(records)}")
            logger.info(
                "Relaxed candidates: "
                f"{len(blast_classification.candidates_relaxed_records)}"
            )
            logger.info(
                "Positive all-source candidates: "
                f"{len(blast_classification.positive_all_sources_records)}"
            )
            logger.info(
                "Negative-unmatched proteins: "
                f"{len(blast_classification.negative_unmatched_records)}"
            )
            logger.info(f"No-hit proteins: {len(blast_classification.no_hit_records)}")
            logger.info(
                f"Negative-hit proteins: {len(blast_classification.negative_hit_records)}"
            )
            logger.info(
                "Strong negative-hit proteins: "
                f"{len(blast_classification.negative_strong_hit_records)}"
            )
            logger.info(
                "Medium negative-hit proteins: "
                f"{len(blast_classification.negative_medium_hit_records)}"
            )
            logger.info(
                "Weak negative-hit proteins: "
                f"{len(blast_classification.negative_weak_hit_records)}"
            )

        cache = JsonCache(config.paths.cache_dir)

        with logger.section("Annotation target settings"):
            for target_name, target in config.annotation_targets.items():
                enabled_sources = [
                    source
                    for source in ("gff", "pfam", "uniprot", "alphafold")
                    if bool(getattr(target, source, False))
                ]
                logger.info(
                    f"{target_name}: "
                    f"{', '.join(enabled_sources) if enabled_sources else 'none'}"
                )

        with logger.section("GFF old locus tag annotation"):
            gff_path = config.paths.gff_file
            gff_records = _records_for_annotation_step(
                blast_classification,
                config.annotation_targets,
                "gff",
            )
            if gff_path is None:
                logger.info("No optional GFF file is configured; skipping GFF annotation.")
            elif not gff_path.exists():
                logger.info(f"Optional GFF file was not found: {gff_path}")
                logger.info("Skipping GFF annotation; the pipeline will continue.")
            elif not gff_records:
                logger.info("No Excel classification sheets are enabled for GFF annotation.")
            else:
                logger.info(f"GFF annotation is enabled: {gff_path}")
                logger.info(f"Records enabled for GFF annotation: {len(gff_records)}")
                with logger.timer("GFF old locus tag annotation"):
                    gff_mapping = load_gff_locus_map(gff_path)
                    updated_records = annotate_records_with_gff_locus_tags(
                        gff_records,
                        gff_mapping,
                    )

                logger.info(
                    f"GFF protein_id to locus tag mappings loaded: {len(gff_mapping)}"
                )
                logger.info(f"Records updated from GFF: {updated_records}")
                if gff_mapping and updated_records == 0:
                    selected_examples = list(gff_records.keys())[:5]
                    gff_key_examples = list(gff_mapping.keys())[:5]
                    logger.warning(
                        "GFF mappings were loaded but no selected records matched."
                    )
                    logger.warning(
                        "Example selected protein IDs: "
                        f"{', '.join(selected_examples) if selected_examples else 'none'}"
                    )
                    logger.warning(
                        "Example GFF protein ID keys: "
                        f"{', '.join(gff_key_examples) if gff_key_examples else 'none'}"
                    )

        with logger.section("CDD domain annotation"):
            if config.annotation.enable_cdd:
                logger.info("CDD annotation is enabled.")
                cdd_target_records = _build_cdd_target_records(
                    records, config, blast_classification
                )
                extra_target_count = len(cdd_target_records) - len(records)
                if extra_target_count > 0:
                    logger.info(
                        "Interaction scoring is enabled: extending CDD annotation to "
                        f"{len(cdd_target_records)} records total "
                        f"({extra_target_count} beyond the Candidates bucket)."
                    )

                with logger.timer("CDD domain annotation"):
                    annotate_records_cdd(
                        cdd_target_records,
                        cache=cache,
                    )

                total_domain_hits = sum(
                    len(record.domains) for record in cdd_target_records.values()
                )
                logger.info(f"Records after CDD annotation: {len(cdd_target_records)}")
                logger.info(f"Total CDD/domain hits: {total_domain_hits}")
                logger.info(
                    "Individual CDD annotation failures are saved in each record's notes."
                )
            else:
                logger.info("CDD annotation is disabled in config.yaml; skipping it.")

        with logger.section("Pfam domain annotation"):
            pfam_records = _records_for_annotation_step(
                blast_classification,
                config.annotation_targets,
                "pfam",
            )
            if config.annotation.enable_pfam and pfam_records:
                logger.info("Pfam annotation is enabled.")
                logger.info(f"Records enabled for Pfam annotation: {len(pfam_records)}")
                with logger.timer("Pfam domain annotation"):
                    annotate_records_pfam(
                        pfam_records,
                        cache=cache,
                        evalue_threshold=config.annotation.pfam_evalue_threshold,
                    )

                total_domain_hits = sum(
                    len(record.domains) for record in records.values()
                )
                logger.info(f"Records after Pfam annotation: {len(records)}")
                logger.info(f"Total domain hits after Pfam: {total_domain_hits}")
                logger.info(
                    "Individual Pfam annotation failures are saved in each record's notes."
                )
            elif config.annotation.enable_pfam:
                logger.info("No Excel classification sheets are enabled for Pfam annotation.")
            else:
                logger.info("Pfam annotation is disabled in config.yaml; skipping it.")

        with logger.section("UniProt and AlphaFold annotation"):
            uniprot_records = _records_for_annotation_step(
                blast_classification,
                config.annotation_targets,
                "uniprot",
            )
            if config.annotation.enable_uniprot and uniprot_records:
                logger.info("UniProt annotation is enabled.")
                logger.info(
                    f"Records enabled for UniProt annotation: {len(uniprot_records)}"
                )
                with logger.timer("UniProt annotation"):
                    annotate_records_uniprot(
                        uniprot_records,
                        cache=cache,
                    )

                logger.info(f"Records after UniProt annotation: {len(records)}")
                logger.info(
                    "Individual UniProt annotation failures are saved in each record's notes."
                )
            elif config.annotation.enable_uniprot:
                logger.info(
                    "No Excel classification sheets are enabled for UniProt annotation."
                )
            else:
                logger.info(
                    "UniProt annotation is disabled in config.yaml; skipping it."
                )

            alphafold_records = _records_for_annotation_step(
                blast_classification,
                config.annotation_targets,
                "alphafold",
            )
            if config.annotation.enable_alphafold and alphafold_records:
                logger.info("AlphaFold annotation is enabled.")
                logger.info(
                    "Records enabled for AlphaFold annotation: "
                    f"{len(alphafold_records)}"
                )
                with logger.timer("AlphaFold annotation"):
                    annotate_records_alphafold(
                        alphafold_records,
                        cache=cache,
                    )

                logger.info(f"Records after AlphaFold annotation: {len(records)}")
                logger.info(
                    "If UniProt accessions are missing, AlphaFold skip notes are saved "
                    "in each record's notes."
                )
            elif config.annotation.enable_alphafold:
                logger.info(
                    "No Excel classification sheets are enabled for AlphaFold annotation."
                )
            else:
                logger.info(
                    "AlphaFold annotation is disabled in config.yaml; skipping it."
                )

        with logger.section("Candidate scoring"):
            with logger.timer("Candidate scoring"):
                records = score_records(records)
                sorted_records = get_sorted_records(records, descending=True)
                records = {record.protein_id: record for record in sorted_records}

            logger.info(f"Records scored: {len(records)}")
            if sorted_records:
                top_candidate = sorted_records[0]
                top_score = (
                    top_candidate.score.total_score if top_candidate.score else 0.0
                )
                logger.info(f"Top candidate: {top_candidate.protein_id}")
                logger.info(f"Top candidate score: {top_score}")
            else:
                logger.info("No candidates were available for scoring.")

        interaction_result = None
        with logger.section("Interaction candidate ranking"):
            if not config.interaction_scoring.enabled:
                logger.info("Interaction scoring is disabled; skipping it.")
            else:
                logger.info("Interaction scoring is enabled.")
                logger.info(
                    "Interaction candidate sources: "
                    + ", ".join(
                        source_name
                        for source_name, enabled in (
                            config.interaction_scoring.candidate_sources.items()
                        )
                        if enabled
                    )
                )
                interaction_result = run_interaction_scoring(
                    config,
                    blast_classification,
                )
                if interaction_result is None:
                    logger.info("No interaction scoring output was created.")
                else:
                    for warning in interaction_result.warnings:
                        logger.warning(warning)
                    logger.info(
                        "Interaction queries: "
                        f"{len(interaction_result.query_rows)}"
                    )
                    logger.info(
                        "Interaction sheets: "
                        f"{len(interaction_result.source_rows)}"
                    )

        with logger.section("Excel output"):
            with logger.timer("Write Excel output"):
                excel_path = write_classification_workbook(
                    config=config,
                    blast_classification=blast_classification,
                    output_path=config.paths.output_excel,
                    interaction_result=interaction_result,
                )

            logger.info(f"Final annotated candidate count: {len(records)}")
            logger.info(
                "Excel sheets written: 01_Index, 02_Final_Score, 03_Candidate_Overview, "
                "04_Score_Breakdown, 05_Sequence_Evidence, 06_Functional_Domain_Evidence, "
                "07_Evolutionary_Evidence, 08_Genomic_Context, 09_Interaction_Evidence, "
                "10_Negative_Evidence, 11_Raw_Audit, 12_Reserved"
            )
            logger.info(f"Excel file written to: {excel_path}")

        logger.summary()
        logger.success("Protein Hunter finished successfully")

    except Exception as exc:
        if hasattr(logger, "exception"):
            logger.exception(exc)
        else:
            logger.error(str(exc))

        logger.summary()
        raise


if __name__ == "__main__":
    main()
