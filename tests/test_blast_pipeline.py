"""Tests for the high-level BLAST candidate pipeline."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from analysis.blast_pipeline import (
    run_blast_candidate_pipeline,
    run_blast_classification_pipeline,
)
from core.models import BlastHit


def make_hit(query_id: str, source: str) -> BlastHit:
    """Create a fake BLAST hit for pipeline tests."""
    return BlastHit(
        query_id=query_id,
        subject_id=f"{source}_subject",
        percent_identity=85.0,
        alignment_length=100,
        evalue=1e-20,
        bitscore=50.0,
        source=source,
    )


def test_blast_candidate_pipeline_calls_positive_and_negative_blast(
    tmp_path: Path,
) -> None:
    """The pipeline should run positive and negative BLAST with correct labels."""
    target_fasta = tmp_path / "targets.faa"
    positive_fasta = tmp_path / "positive.faa"
    negative_fasta = tmp_path / "negative.faa"
    positive_hit = make_hit("protein_1", "positive")
    negative_hit = make_hit("protein_2", "negative")

    with (
        patch(
            "analysis.blast_pipeline.read_fasta_as_components",
            return_value=(
                ["protein_1", "protein_2"],
                {"protein_1": "protein_1 desc", "protein_2": "protein_2 desc"},
                {"protein_1": "MSTN", "protein_2": "AAAA"},
            ),
        ) as read_mock,
        patch(
            "analysis.blast_pipeline.run_blast_pipeline",
            side_effect=[[positive_hit], [negative_hit]],
        ) as blast_mock,
    ):
        records = run_blast_candidate_pipeline(
            target_fasta=target_fasta,
            positive_fasta=positive_fasta,
            negative_fasta=negative_fasta,
            work_dir=tmp_path / "work",
            evalue=1e-10,
            max_target_seqs=5,
            threads=4,
        )

    read_mock.assert_called_once_with(target_fasta)
    assert blast_mock.call_count == 2
    assert blast_mock.call_args_list[0].kwargs == {
        "query_fasta": target_fasta,
        "subject_fasta": positive_fasta,
        "work_dir": (tmp_path / "work").resolve() / "positive",
        "db_name": "positive",
        "source": "positive",
        "evalue": 1e-10,
        "max_target_seqs": 5,
        "threads": 4,
    }
    assert blast_mock.call_args_list[1].kwargs == {
        "query_fasta": target_fasta,
        "subject_fasta": negative_fasta,
        "work_dir": (tmp_path / "work").resolve() / "negative",
        "db_name": "negative",
        "source": "negative",
        "evalue": 1e-10,
        "max_target_seqs": 5,
        "threads": 4,
    }
    assert set(records) == {"protein_1"}
    assert records["protein_1"].description == "protein_1 desc"
    assert records["protein_1"].sequence == "MSTN"
    assert records["protein_1"].positive_hits == [positive_hit]
    assert records["protein_1"].negative_hits == []


def test_blast_candidate_pipeline_returns_filtered_positive_only_records(
    tmp_path: Path,
) -> None:
    """The pipeline should return only records with positive hits and no negatives."""
    positive_hits = [
        make_hit("positive_only", "positive"),
        make_hit("mixed", "positive"),
    ]
    negative_hits = [
        make_hit("mixed", "negative"),
        make_hit("negative_only", "negative"),
    ]

    with (
        patch(
            "analysis.blast_pipeline.read_fasta_as_components",
            return_value=(
                ["positive_only", "mixed", "negative_only", "no_hits"],
                {
                    "positive_only": "positive only",
                    "mixed": "mixed",
                    "negative_only": "negative only",
                    "no_hits": "no hits",
                },
                {
                    "positive_only": "AAAA",
                    "mixed": "CCCC",
                    "negative_only": "GGGG",
                    "no_hits": "TTTT",
                },
            ),
        ),
        patch(
            "analysis.blast_pipeline.run_blast_pipeline",
            side_effect=[positive_hits, negative_hits],
        ),
    ):
        records = run_blast_candidate_pipeline(
            target_fasta=tmp_path / "targets.faa",
            positive_fasta=tmp_path / "positive.faa",
            negative_fasta=tmp_path / "negative.faa",
            work_dir=tmp_path / "work",
        )

    assert set(records) == {"positive_only"}
    assert records["positive_only"].positive_hits == [positive_hits[0]]
    assert records["positive_only"].negative_hits == []


def test_blast_classification_pipeline_returns_all_classification_groups(
    tmp_path: Path,
) -> None:
    """Targets should be split into Candidates, No_hit, and Negative_hit groups."""
    positive_hits = [
        make_hit("A_positive_only", "positive"),
        make_hit("D_both", "positive"),
    ]
    negative_hits = [
        make_hit("C_negative_only", "negative"),
        make_hit("D_both", "negative"),
    ]

    with (
        patch(
            "analysis.blast_pipeline.read_fasta_as_components",
            return_value=(
                ["A_positive_only", "B_no_hits", "C_negative_only", "D_both"],
                {
                    "A_positive_only": "A desc",
                    "B_no_hits": "B desc",
                    "C_negative_only": "C desc",
                    "D_both": "D desc",
                },
                {
                    "A_positive_only": "AAAA",
                    "B_no_hits": "BBBB",
                    "C_negative_only": "CCCC",
                    "D_both": "DDDD",
                },
            ),
        ),
        patch(
            "analysis.blast_pipeline.run_blast_pipeline",
            side_effect=[positive_hits, negative_hits],
        ),
    ):
        result = run_blast_classification_pipeline(
            target_fasta=tmp_path / "targets.faa",
            positive_fasta=tmp_path / "positive.faa",
            negative_fasta=tmp_path / "negative.faa",
            work_dir=tmp_path / "work",
        )

    assert set(result.all_records) == {
        "A_positive_only",
        "B_no_hits",
        "C_negative_only",
        "D_both",
    }
    assert set(result.positive_only_records) == {"A_positive_only"}
    assert set(result.negative_unmatched_records) == {
        "A_positive_only",
        "B_no_hits",
    }
    assert set(result.no_hit_records) == {"B_no_hits"}
    assert set(result.negative_hit_records) == {"C_negative_only", "D_both"}
