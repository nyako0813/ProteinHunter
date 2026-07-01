"""
Protein Hunter v5
Main entry point.
"""

from __future__ import annotations

from pathlib import Path

from core.startup import StartupChecker


def main() -> None:
    """Run Protein Hunter."""

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
        from analysis.blast_pipeline import run_blast_candidate_pipeline
        from analysis.input_summary import format_input_summary, summarize_input_fastas
        from analysis.scoring import get_sorted_records, score_records
        from config import CONFIG
        from core.cache import JsonCache
        from output.excel import write_records_to_excel

        logger.info("Protein Hunter started")
        logger.success("Startup check passed")

        blast_work_dir = Path("data") / "temp" / "blast"

        with logger.section("Configuration"):
            logger.info(f"Target FASTA: {CONFIG.paths.target_fasta}")
            logger.info(f"Positive FASTA: {CONFIG.paths.positive_fasta}")
            logger.info(f"Negative FASTA: {CONFIG.paths.negative_fasta}")
            logger.info(f"Excel output: {CONFIG.paths.output_excel}")
            logger.info(f"BLAST work directory: {blast_work_dir}")
            logger.info(f"Cache directory: {CONFIG.paths.cache_dir}")

            input_summary = summarize_input_fastas(
                target_fasta=CONFIG.paths.target_fasta,
                positive_fasta=CONFIG.paths.positive_fasta,
                negative_fasta=CONFIG.paths.negative_fasta,
            )
            logger.info("Input FASTA summary:")
            for line in format_input_summary(input_summary):
                logger.info(line)

        with logger.section("BLAST candidate search"):
            with logger.timer("BLAST candidate pipeline"):
                records = run_blast_candidate_pipeline(
                    target_fasta=CONFIG.paths.target_fasta,
                    positive_fasta=CONFIG.paths.positive_fasta,
                    negative_fasta=CONFIG.paths.negative_fasta,
                    work_dir=blast_work_dir,
                    evalue=CONFIG.blast.evalue,
                    max_target_seqs=CONFIG.blast.max_target_seqs,
                    threads=CONFIG.blast.threads,
                )

            logger.info(f"BLAST positive-only candidates: {len(records)}")

        cache = JsonCache(CONFIG.paths.cache_dir)

        with logger.section("CDD domain annotation"):
            if CONFIG.annotation.enable_cdd:
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
            if CONFIG.annotation.enable_pfam:
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
            if CONFIG.annotation.enable_uniprot:
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

            if CONFIG.annotation.enable_alphafold:
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
                    output_path=CONFIG.paths.output_excel,
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
