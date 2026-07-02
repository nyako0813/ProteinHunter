"""Tests for GFF old/locus tag annotation helpers."""

from __future__ import annotations

from pathlib import Path

from annotation.gff import (
    annotate_records_with_gff_locus_tags,
    load_gff_locus_map,
    normalize_protein_id,
)
from core.models import ProteinRecord


def write_gff(path: Path, attribute_text: str, feature_type: str = "CDS") -> Path:
    """Write one small GFF row with the provided attributes."""
    path.write_text(
        f"seqid\tRefSeq\t{feature_type}\t1\t100\t.\t+\t0\t{attribute_text}\n",
        encoding="utf-8",
    )
    return path


def test_gff_protein_id_and_locus_tag_maps_correctly(tmp_path: Path) -> None:
    """protein_id plus locus_tag should map a protein to MA_####."""
    gff = write_gff(
        tmp_path / "genome.gff",
        "protein_id=WP_011020109.1;locus_tag=MA_0050;product=ATPase",
    )

    mapping = load_gff_locus_map(gff)

    assert mapping["WP_011020109.1"] == "MA_0050"


def test_gff_dbxref_and_old_locus_tag_maps_correctly(tmp_path: Path) -> None:
    """Dbxref Genbank protein IDs should be normalized and mapped."""
    gff = write_gff(
        tmp_path / "genome.gff",
        "Dbxref=Genbank:WP_011020109.1;old_locus_tag=MA_0050",
    )

    mapping = load_gff_locus_map(gff)

    assert mapping["WP_011020109.1"] == "MA_0050"


def test_gff_id_and_parent_maps_correctly(tmp_path: Path) -> None:
    """ID=cds-* should map using the MA tag found in Parent."""
    gff = write_gff(
        tmp_path / "genome.gff",
        "ID=cds-WP_011020109.1;Parent=gene-MA_0050;product=ATPase",
    )

    mapping = load_gff_locus_map(gff)

    assert mapping["WP_011020109.1"] == "MA_0050"


def test_gff_annotation_updates_matching_records() -> None:
    """GFF mapping should set old_locus_tag on matching ProteinRecord objects."""
    records = {
        "WP_011020109.1": ProteinRecord(protein_id="WP_011020109.1"),
        "WP_missing": ProteinRecord(protein_id="WP_missing"),
    }

    updated = annotate_records_with_gff_locus_tags(
        records,
        {"WP_011020109.1": "MA_0050"},
    )

    assert updated == 1
    assert records["WP_011020109.1"].old_locus_tag == "MA_0050"
    assert records["WP_missing"].old_locus_tag is None


def test_normalize_protein_id_supports_common_gff_prefixes() -> None:
    """Common GFF prefixes should normalize to the protein accession."""
    assert normalize_protein_id("WP_011020109.1") == "WP_011020109.1"
    assert normalize_protein_id("cds-WP_011020109.1") == "WP_011020109.1"
    assert normalize_protein_id("Genbank:WP_011020109.1") == "WP_011020109.1"
