"""Tests for the startup reference-genome check (core/reference_genomes.py)."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from core.reference_genomes import (
    DEFAULT_MANIFEST_PATH,
    ReferenceGenomeRecord,
    check_reference_genomes,
    describe_reference_genomes,
    load_manifest,
    majority_species,
)

REPO_ROOT = Path(__file__).resolve().parent.parent


def write_genome(root: Path, category: str, label: str, accession: str, organism: str, n: int = 3, salt: str = "") -> None:
    folder = root / category / label / "ncbi_dataset" / "data" / accession
    folder.mkdir(parents=True)
    lines = [f">WP_{i}{salt}.1 protein {i} [{organism}]\nMKV{'A' * i}\n" for i in range(1, n + 1)]
    (folder / "protein.faa").write_text("".join(lines))


def rec(category: str, label: str, accession: str, organism: str, md5: str) -> ReferenceGenomeRecord:
    return ReferenceGenomeRecord(category, label, accession, organism, 10, md5)


MANIFEST = {
    "positive": {"Pos_a": {"accession": "GCF_1", "organism": "Pos aaa"}},
    "negative": {"Neg_b": {"accession": "GCF_2", "organism": "Neg bbb"}},
}


def test_consistent_genomes_produce_no_findings() -> None:
    records = [rec("positive", "Pos_a", "GCF_1", "Pos aaa", "m1"), rec("negative", "Neg_b", "GCF_2", "Neg bbb", "m2")]

    assert check_reference_genomes(records, MANIFEST) == []


def test_identical_content_in_two_folders_is_flagged_even_without_manifest() -> None:
    records = [rec("positive", "Pos_a", "GCF_1", "Pos aaa", "same"), rec("negative", "Neg_b", "GCF_1", "Pos aaa", "same")]

    findings = check_reference_genomes(records, None)

    assert any("identical protein.faa content" in f and "positive/Pos_a" in f and "negative/Neg_b" in f for f in findings)
    assert any("both a positive" in f for f in findings)


def test_wrong_organism_in_labelled_folder_is_flagged() -> None:
    """The original mistake: folder Neg_b holds a different species than expected."""
    records = [rec("positive", "Pos_a", "GCF_1", "Pos aaa", "m1"), rec("negative", "Neg_b", "GCF_2", "Pos aaa", "m2")]

    findings = check_reference_genomes(records, MANIFEST)

    assert any("negative/Neg_b" in f and "'Pos aaa'" in f and "expected 'Neg bbb'" in f for f in findings)


def test_wrong_accession_is_flagged() -> None:
    records = [rec("positive", "Pos_a", "GCF_9", "Pos aaa", "m1"), rec("negative", "Neg_b", "GCF_2", "Neg bbb", "m2")]

    findings = check_reference_genomes(records, MANIFEST)

    assert findings == ["positive/Pos_a holds assembly GCF_9, expected GCF_1"]


def test_missing_expected_genome_is_flagged() -> None:
    """A machine that lacks a genome (e.g. one that was git-ignored) must say so."""
    findings = check_reference_genomes([rec("positive", "Pos_a", "GCF_1", "Pos aaa", "m1")], MANIFEST)

    assert findings == ["negative/Neg_b is expected (GCF_2) but has no protein.faa on this machine"]


def test_unexpected_genome_is_flagged() -> None:
    records = [
        rec("positive", "Pos_a", "GCF_1", "Pos aaa", "m1"),
        rec("negative", "Neg_b", "GCF_2", "Neg bbb", "m2"),
        rec("negative", "Neg_extra", "GCF_3", "Neg ccc", "m3"),
    ]

    findings = check_reference_genomes(records, MANIFEST)

    assert len(findings) == 1
    assert "negative/Neg_extra" in findings[0] and "not listed" in findings[0]


def test_describe_reference_genomes_reads_real_folder_layout(tmp_path: Path) -> None:
    write_genome(tmp_path, "positive", "Pos_a", "GCF_1", "Pos aaa", n=3)
    write_genome(tmp_path, "negative", "Neg_b", "GCF_2", "Neg bbb", n=5)

    records = describe_reference_genomes({"positive": tmp_path / "positive", "negative": tmp_path / "negative"})

    assert [(r.category, r.label, r.accession, r.organism, r.protein_count) for r in records] == [
        ("positive", "Pos_a", "GCF_1", "Pos aaa", 3),
        ("negative", "Neg_b", "GCF_2", "Neg bbb", 5),
    ]
    assert check_reference_genomes(records, MANIFEST) == []


def test_end_to_end_duplicate_genome_copied_into_negative_folder(tmp_path: Path) -> None:
    write_genome(tmp_path, "positive", "Pos_a", "GCF_1", "Pos aaa")
    write_genome(tmp_path, "negative", "Neg_b", "GCF_1", "Pos aaa")  # same genome, wrong folder

    records = describe_reference_genomes({"positive": tmp_path / "positive", "negative": tmp_path / "negative"})
    findings = check_reference_genomes(records, MANIFEST)

    assert any("identical protein.faa content" in f for f in findings)
    assert any("negative/Neg_b" in f and "expected 'Neg bbb'" in f for f in findings)


def test_majority_species_uses_trailing_organism_tag(tmp_path: Path) -> None:
    fasta = tmp_path / "p.faa"
    fasta.write_text(">a x [Foo bar]\nM\n>b y [Foo bar]\nM\n>c z [Foo]\nM\n")

    assert majority_species(fasta) == ("Foo bar", 3)


@pytest.mark.skipif(
    not (REPO_ROOT / "data" / "databases" / "positive").is_dir()
    or not (REPO_ROOT / "data" / "databases" / "negative").is_dir(),
    reason="reference genomes not present",
)
def test_committed_reference_genomes_match_committed_manifest() -> None:
    """The repo's own data must satisfy the repo's own manifest (catches a bad or machine-local genome)."""
    records = describe_reference_genomes(
        {
            "positive": REPO_ROOT / "data" / "databases" / "positive",
            "negative": REPO_ROOT / "data" / "databases" / "negative",
        }
    )

    assert check_reference_genomes(records, load_manifest(REPO_ROOT / DEFAULT_MANIFEST_PATH)) == []


@pytest.mark.skipif(shutil.which("git") is None or not (REPO_ROOT / ".git").exists(), reason="needs a git checkout")
def test_manifest_genomes_are_not_gitignored() -> None:
    """A manifest genome that .gitignore hides exists on one machine only.

    negative/thermoplasma_acidophilum was git-ignored (2026-09-09..09-20), so a
    fresh clone silently ran with one fewer negative reference genome.
    """
    manifest = load_manifest(REPO_ROOT / DEFAULT_MANIFEST_PATH)
    ignored = []
    for category, genomes in manifest.items():
        for label, spec in genomes.items():
            path = f"data/databases/{category}/{label}/ncbi_dataset/data/{spec['accession']}/protein.faa"
            result = subprocess.run(["git", "check-ignore", "-q", path], cwd=REPO_ROOT)
            if result.returncode == 0:
                ignored.append(path)

    assert ignored == []


# ---------------------------------------------------------------------------
# main._log_reference_genome_check -> RunProvenance.reference_genome_check_passed
# ---------------------------------------------------------------------------


class RecordingLogger:
    def __init__(self) -> None:
        self.warnings: list[str] = []

    def info(self, message: str) -> None:
        pass

    def warning(self, message: str) -> None:
        self.warnings.append(message)


def _run_main_check(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, manifest: dict | None) -> tuple[bool | None, RecordingLogger]:
    import yaml

    import core.reference_genomes as reference_genomes
    from main import _log_reference_genome_check

    manifest_path = tmp_path / "manifest.yaml"
    if manifest is not None:
        manifest_path.write_text(yaml.safe_dump(manifest))
    monkeypatch.setattr(reference_genomes, "DEFAULT_MANIFEST_PATH", manifest_path)
    logger = RecordingLogger()
    result = _log_reference_genome_check(
        logger, {"positive": tmp_path / "positive", "negative": tmp_path / "negative"}
    )
    return result, logger


def test_main_check_returns_true_when_genomes_match_the_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    write_genome(tmp_path, "positive", "Pos_a", "GCF_1", "Pos aaa")
    write_genome(tmp_path, "negative", "Neg_b", "GCF_2", "Neg bbb", salt="b")

    result, logger = _run_main_check(tmp_path, monkeypatch, {"version": "v1", **MANIFEST})

    assert result is True
    assert logger.warnings == []


def test_main_check_returns_false_when_findings_were_logged(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    write_genome(tmp_path, "positive", "Pos_a", "GCF_1", "Pos aaa")
    write_genome(tmp_path, "negative", "Neg_b", "GCF_1", "Pos aaa")  # the original mix-up

    result, logger = _run_main_check(tmp_path, monkeypatch, {"version": "v1", **MANIFEST})

    assert result is False
    assert logger.warnings


def test_main_check_returns_none_when_nothing_could_be_verified(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    write_genome(tmp_path, "positive", "Pos_a", "GCF_1", "Pos aaa")
    write_genome(tmp_path, "negative", "Neg_b", "GCF_2", "Neg bbb", salt="b")

    result, logger = _run_main_check(tmp_path, monkeypatch, None)  # no manifest, no findings

    assert result is None
    assert any("manifest not found" in w for w in logger.warnings)
