"""Directory-based FASTA source helpers for ProteinHunter."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

from core.exceptions import FileValidationError


@dataclass(frozen=True)
class FastaSource:
    label: str
    path: Path


@dataclass(frozen=True)
class DirectoryFastaResult:
    category: str
    directory: Path
    combined_fasta: Path
    sources: tuple[FastaSource, ...]
    source_labels: tuple[str, ...]
    skipped_folders: tuple[str, ...]
    multiple_file_labels: tuple[str, ...]
    duplicate_ids: tuple[str, ...]
    record_count: int


def discover_fasta_sources(
    directory: str | Path,
    category: str,
) -> tuple[list[FastaSource], list[str], list[str]]:
    """Find protein.faa files recursively under immediate source folders."""
    root = Path(directory)

    if not root.exists():
        raise FileValidationError(f"The {category} FASTA directory was not found: {root}")

    if not root.is_dir():
        raise FileValidationError(f"The {category} FASTA path is not a directory: {root}")

    sources: list[FastaSource] = []
    skipped: list[str] = []
    multiple_file_labels: list[str] = []

    for child in sorted((path for path in root.iterdir() if path.is_dir()), key=lambda path: path.name):
        protein_fastas = sorted(child.rglob("protein.faa"))
        protein_fastas = [path for path in protein_fastas if path.is_file()]
        if not protein_fastas:
            skipped.append(child.name)
            continue

        if len(protein_fastas) > 1:
            multiple_file_labels.append(child.name)

        for protein_faa in protein_fastas:
            sources.append(FastaSource(label=child.name, path=protein_faa))

    if not sources:
        raise FileValidationError(
            f"No valid protein.faa files were found for {category} under {root}."
        )

    return sources, skipped, multiple_file_labels


def combine_fasta_sources(
    sources: list[FastaSource],
    output_path: str | Path,
    category: str,
) -> tuple[Path, tuple[str, ...], int]:
    """Write a combined FASTA while preserving record IDs."""
    combined_path = Path(output_path)
    combined_path.parent.mkdir(parents=True, exist_ok=True)

    seen_ids: set[str] = set()
    duplicate_ids: list[str] = []
    total_records = 0

    with combined_path.open("w", encoding="utf-8") as handle:
        for source in sources:
            records = list(SeqIO.parse(source.path, "fasta"))
            if not records:
                raise FileValidationError(
                    f"The {category} FASTA source contains no readable records: {source.path}"
                )

            for record in records:
                if record.id in seen_ids and record.id not in duplicate_ids:
                    duplicate_ids.append(record.id)
                seen_ids.add(record.id)
                record.description = f"{record.description} [source={source.label}]"
                SeqIO.write(record, handle, "fasta")
                total_records += 1

    if total_records == 0:
        raise FileValidationError(
            f"No readable FASTA records were found for {category}."
        )

    return combined_path, tuple(duplicate_ids), total_records


def prepare_directory_fasta(
    directory: str | Path,
    category: str,
    combined_dir: str | Path = Path("data") / "temp" / "combined",
) -> DirectoryFastaResult:
    """Discover directory sources and create one combined FASTA for a category."""
    sources, skipped, multiple_file_labels = discover_fasta_sources(directory, category)
    combined_path = Path(combined_dir) / f"{category}.combined.faa"
    combined_fasta, duplicate_ids, record_count = combine_fasta_sources(
        sources=sources,
        output_path=combined_path,
        category=category,
    )

    return DirectoryFastaResult(
        category=category,
        directory=Path(directory),
        combined_fasta=combined_fasta,
        sources=tuple(sources),
        source_labels=tuple(dict.fromkeys(source.label for source in sources)),
        skipped_folders=tuple(skipped),
        multiple_file_labels=tuple(multiple_file_labels),
        duplicate_ids=duplicate_ids,
        record_count=record_count,
    )


__all__: tuple[str, ...] = (
    "DirectoryFastaResult",
    "FastaSource",
    "combine_fasta_sources",
    "discover_fasta_sources",
    "prepare_directory_fasta",
)
