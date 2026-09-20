"""Guards on the committed BLAST reference genomes under data/databases.

Regression for a real mistake: the folder labelled ``Sulfolobus_solfataricus``
(negative) actually held a byte-identical copy of the *positive*
``Methanococcus_maripaludis`` genome (GCF_002945325.1). Because
ortholog_filter treats any negative hit as evidence against a candidate, a
positive genome doubling as a negative one silently inflated every
negative_hit_strength and emptied the positive-only candidate pool -- with no
error anywhere, since the pipeline only ever sees the combined FASTA. These
tests read the tracked files directly so such a mix-up fails loudly instead.
"""

from __future__ import annotations

import csv
import hashlib
import re
from collections import Counter
from pathlib import Path

import pytest

DATABASES = Path(__file__).resolve().parent.parent / "data" / "databases"
CATEGORIES = ("positive", "negative")
_ORGANISM = re.compile(r"\[([^\]]+)\]\s*$")


def _genome_files() -> list[tuple[str, str, Path]]:
    """(category, folder label, protein.faa path) for every reference genome present."""
    found = []
    for category in CATEGORIES:
        for faa in sorted((DATABASES / category).glob("*/ncbi_dataset/data/*/protein.faa")):
            label = faa.relative_to(DATABASES / category).parts[0]
            found.append((category, label, faa))
    return found


GENOMES = _genome_files()
pytestmark = pytest.mark.skipif(not GENOMES, reason="no reference genomes present under data/databases")


def _majority_species(faa: Path) -> str:
    """Most common 'Genus species' in the FASTA headers' trailing [organism] tag."""
    counts: Counter[str] = Counter()
    with faa.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                match = _ORGANISM.search(line)
                if match:
                    counts[" ".join(match.group(1).split()[:2])] += 1
    assert counts, f"{faa}: no [organism] tags found in FASTA headers"
    return counts.most_common(1)[0][0]


def _md5(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()


@pytest.mark.parametrize(("category", "label", "faa"), GENOMES, ids=lambda v: str(v) if isinstance(v, str) else None)
def test_folder_label_species_matches_fasta_content(category: str, label: str, faa: Path) -> None:
    """The folder's species epithet must be the FASTA's organism.

    Only the epithet (last word) is compared, so a genus reclassification
    (Sulfolobus -> Saccharolobus solfataricus) is fine, but a different
    species under the label is not.
    """
    epithet = label.split("_")[-1].lower()
    species = _majority_species(faa)

    assert epithet in species.lower(), (
        f"{category}/{label} is labelled '{epithet}' but {faa.name} is '{species}'"
    )


def test_no_species_is_both_positive_and_negative() -> None:
    species_by_category: dict[str, dict[str, str]] = {c: {} for c in CATEGORIES}
    for category, label, faa in GENOMES:
        species_by_category[category][_majority_species(faa)] = label

    overlap = set(species_by_category["positive"]) & set(species_by_category["negative"])

    assert not overlap, {
        species: (species_by_category["positive"][species], species_by_category["negative"][species])
        for species in overlap
    }


def test_no_two_reference_genome_folders_share_identical_protein_fasta() -> None:
    by_hash: dict[str, list[str]] = {}
    for category, label, faa in GENOMES:
        by_hash.setdefault(_md5(faa), []).append(f"{category}/{label}")

    duplicates = [folders for folders in by_hash.values() if len(folders) > 1]

    assert not duplicates, duplicates


@pytest.mark.parametrize(("category", "label", "faa"), GENOMES, ids=lambda v: str(v) if isinstance(v, str) else None)
def test_data_summary_matches_folder_when_present(category: str, label: str, faa: Path) -> None:
    summary = faa.parent.parent / "data_summary.tsv"
    if not summary.exists():
        pytest.skip("no data_summary.tsv for this genome")
    with summary.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    accession = faa.parent.name
    refseq_rows = [row for row in rows if row["Assembly Accession"] == accession]
    assert refseq_rows, f"{summary}: no row for {accession} (folder {category}/{label})"
    organism = refseq_rows[0]["Organism Scientific Name"].lower()
    assert label.split("_")[-1].lower() in organism, (
        f"{category}/{label}: data_summary.tsv says '{organism}' for {accession}"
    )
