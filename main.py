"""
Protein Hunter v5
Main entry point.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from collections.abc import Sequence

from core.startup import StartupChecker


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
        from analysis.blast_pipeline import run_blast_candidate_pipeline
        from analysis.input_summary import format_input_summary, summarize_input_fastas
        from analysis.scoring import get_sorted_records, score_records
        from config import load_config
        from core.cache import JsonCache
        from output.excel import write_records_to_excel

        config_path = Path(args.config)
        config = load_config(config_path)

        logger.info("Protein Hunter started")
        logger.success("Startup check passed")
        logger.info(f"Using config file: {config_path}")

        blast_work_dir = Path("data") / "temp" / "blast"

        with logger.section("Configuration"):
            logger.info(f"Target FASTA: {config.paths.target_fasta}")
            logger.info(f"Positive FASTA: {config.paths.positive_fasta}")
            logger.info(f"Negative FASTA: {config.paths.negative_fasta}")
            logger.info(f"Excel output: {config.paths.output_excel}")
            logger.info(f"BLAST work directory: {blast_work_dir}")
            logger.info(f"Cache directory: {config.paths.cache_dir}")
            if config.paths.gff_file is not None:
                logger.info(f"Optional GFF file: {config.paths.gff_file}")
            else:
                logger.info("Optional GFF file: not configured")

            input_summary = summarize_input_fastas(
                target_fasta=config.paths.target_fasta,
                positive_fasta=config.paths.positive_fasta,
                negative_fasta=config.paths.negative_fasta,
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
                records = run_blast_candidate_pipeline(
                    target_fasta=config.paths.target_fasta,
                    positive_fasta=config.paths.positive_fasta,
                    negative_fasta=config.paths.negative_fasta,
                    work_dir=blast_work_dir,
                    evalue=config.blast.evalue,
                    max_target_seqs=config.blast.max_target_seqs,
                    threads=config.blast.threads,
                )

            logger.info(f"BLAST positive-only candidates: {len(records)}")

        cache = JsonCache(config.paths.cache_dir)

        with logger.section("GFF old locus tag annotation"):
            gff_path = config.paths.gff_file
            if gff_path is None:
                logger.info("No optional GFF file is configured; skipping GFF annotation.")
            elif not gff_path.exists():
                logger.info(f"Optional GFF file was not found: {gff_path}")
                logger.info("Skipping GFF annotation; the pipeline will continue.")
            else:
                logger.info(f"GFF annotation is enabled: {gff_path}")
                with logger.timer("GFF old locus tag annotation"):
                    gff_mapping = load_gff_locus_map(gff_path)
                    updated_records = annotate_records_with_gff_locus_tags(
                        records,
                        gff_mapping,
                    )

                logger.info(
                    f"GFF protein_id to locus tag mappings loaded: {len(gff_mapping)}"
                )
                logger.info(f"Candidate records updated from GFF: {updated_records}")

        with logger.section("CDD domain annotation"):
            if config.annotation.enable_cdd:
                logger.info("CDD annotation is enabled.")
                with logger.timer("CDD domain annotation"):
                    records = annotate_records_cdd(
                        records,
                        cache=cache,
                    )

                total_domain_hits = sum(
                    len(record.domains) for record in records.values()
                )
                logger.info(f"Records after CDD annotation: {len(records)}")
                logger.info(f"Total CDD/domain hits: {total_domain_hits}")
                logger.info(
                    "Individual CDD annotation failures are saved in each record's notes."
                )
            else:
                logger.info("CDD annotation is disabled in config.yaml; skipping it.")

        with logger.section("Pfam domain annotation"):
            if config.annotation.enable_pfam:
                logger.info("Pfam annotation is enabled.")
                with logger.timer("Pfam domain annotation"):
                    records = annotate_records_pfam(
                        records,
                        cache=cache,
                    )

                total_domain_hits = sum(
                    len(record.domains) for record in records.values()
                )
                logger.info(f"Records after Pfam annotation: {len(records)}")
                logger.info(f"Total domain hits after Pfam: {total_domain_hits}")
                logger.info(
                    "Individual Pfam annotation failures are saved in each record's notes."
                )
            else:
                logger.info("Pfam annotation is disabled in config.yaml; skipping it.")

        with logger.section("UniProt and AlphaFold annotation"):
            if config.annotation.enable_uniprot:
                logger.info("UniProt annotation is enabled.")
                with logger.timer("UniProt annotation"):
                    records = annotate_records_uniprot(
                        records,
                        cache=cache,
                    )

                logger.info(f"Records after UniProt annotation: {len(records)}")
                logger.info(
                    "Individual UniProt annotation failures are saved in each record's notes."
                )
            else:
                logger.info(
                    "UniProt annotation is disabled in config.yaml; skipping it."
                )

            if config.annotation.enable_alphafold:
                logger.info("AlphaFold annotation is enabled.")
                with logger.timer("AlphaFold annotation"):
                    records = annotate_records_alphafold(
                        records,
                        cache=cache,
                    )

                logger.info(f"Records after AlphaFold annotation: {len(records)}")
                logger.info(
                    "If UniProt accessions are missing, AlphaFold skip notes are saved "
                    "in each record's notes."
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

        with logger.section("Excel output"):
            with logger.timer("Write Excel output"):
                excel_path = write_records_to_excel(
                    records=records,
                    output_path=config.paths.output_excel,
                )

            logger.info(f"Final annotated candidate count: {len(records)}")
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
