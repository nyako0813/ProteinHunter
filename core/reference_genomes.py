"""Startup check: do the BLAST reference genome folders match what is expected?

The pipeline only sees the combined FASTA per category, so a wrong, missing,
duplicated, or machine-local reference genome changes every classification
without raising anything. ``config/reference_genomes.v1.yaml`` records the
expected composition; ``check_reference_genomes`` compares the folders on disk
against it and against each other, and ``describe_reference_genomes`` returns
one identifying record per genome (accession, organism, protein count, md5)
suitable for logging next to a run's results.
"""

from __future__ import annotations

import hashlib
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import yaml

from core.fasta_sources import discover_fasta_sources

DEFAULT_MANIFEST_PATH = Path("config") / "reference_genomes.v1.yaml"

_ORGANISM_TAG = re.compile(r"\[([^\]]+)\]\s*$")


@dataclass(frozen=True)
class ReferenceGenomeRecord:
    category: str
    label: str
    accession: str
    organism: str
    protein_count: int
    md5: str


def majority_species(fasta_path: Path) -> tuple[str, int]:
    """Return the most common 'Genus species' in the FASTA's trailing [organism] tags, and the header count."""
    counts: Counter[str] = Counter()
    total = 0
    with fasta_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            total += 1
            match = _ORGANISM_TAG.search(line)
            if match:
                counts[" ".join(match.group(1).split()[:2])] += 1
    return (counts.most_common(1)[0][0] if counts else "", total)


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def describe_reference_genomes(directories: dict[str, Path]) -> list[ReferenceGenomeRecord]:
    """One record per protein.faa found under each category directory ({'positive': dir, ...})."""
    records = []
    for category, directory in directories.items():
        sources, _skipped, _multiple = discover_fasta_sources(directory, category)
        for source in sources:
            organism, count = majority_species(source.path)
            records.append(
                ReferenceGenomeRecord(
                    category=category,
                    label=source.label,
                    accession=source.path.parent.name,
                    organism=organism,
                    protein_count=count,
                    md5=_md5(source.path),
                )
            )
    return records


def load_manifest(path: Path) -> dict[str, dict[str, dict[str, str]]]:
    with path.open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return {category: dict(data.get(category) or {}) for category in ("positive", "negative")}


def check_reference_genomes(
    records: list[ReferenceGenomeRecord],
    manifest: dict[str, dict[str, dict[str, str]]] | None,
) -> list[str]:
    """Return human-readable findings; empty means everything is consistent.

    Checks that need no manifest (the same file or the same species used as
    both a positive and a negative reference, or twice) always run; the
    manifest checks (missing/unexpected folder, accession, organism) run only
    when a manifest is given.
    """
    findings: list[str] = []

    by_md5: dict[str, list[str]] = {}
    for record in records:
        by_md5.setdefault(record.md5, []).append(f"{record.category}/{record.label}")
    for folders in by_md5.values():
        if len(folders) > 1:
            findings.append(f"identical protein.faa content in {', '.join(folders)}")

    species_by_category: dict[str, dict[str, str]] = {}
    for record in records:
        species_by_category.setdefault(record.category, {})[record.organism] = record.label
    positive = species_by_category.get("positive", {})
    negative = species_by_category.get("negative", {})
    for species in sorted((set(positive) & set(negative)) - {""}):
        findings.append(
            f"{species} is used as both a positive ({positive[species]}) and a negative ({negative[species]}) "
            "reference genome"
        )

    if manifest is None:
        return findings

    present = {(record.category, record.label): record for record in records}
    for category in ("positive", "negative"):
        expected = manifest.get(category, {})
        for label, spec in expected.items():
            record = present.get((category, label))
            if record is None:
                findings.append(
                    f"{category}/{label} is expected ({spec.get('accession')}) but has no protein.faa on this machine"
                )
                continue
            if record.accession != spec.get("accession"):
                findings.append(
                    f"{category}/{label} holds assembly {record.accession}, expected {spec.get('accession')}"
                )
            if record.organism.lower() != str(spec.get("organism", "")).lower():
                findings.append(
                    f"{category}/{label} FASTA says '{record.organism}', expected '{spec.get('organism')}'"
                )
        for record in records:
            if record.category == category and record.label not in expected:
                findings.append(
                    f"{category}/{record.label} ({record.accession}, {record.organism}) is not listed in the "
                    "reference genome manifest"
                )
    return findings
