"""Tests for GFF old/locus tag annotation helpers."""

from __future__ import annotations

from pathlib import Path

from annotation.gff import (
    annotate_records_with_gff_locus_tags,
    load_gff_locus_map,
    normalize_protein_id,
    parse_gff_attributes,
)
from core.models import ProteinRecord


def write_gff(path: Path, attribute_text: str, feature_type: str = "CDS") -> Path:
    """Write one small GFF row with the provided attributes."""
    path.write_text(
        f"seqid\tRefSeq\t{feature_type}\t1\t100\t.\t+\t0\t{attribute_text}\n",
        encoding="utf-8",
    )
    return path


def write_gff_lines(path: Path, lines: list[str]) -> Path:
    """Write GFF rows exactly as provided."""
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
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


def test_gff_gene_parent_old_locus_tag_maps_cds_protein_id(
    tmp_path: Path,
) -> None:
    """CDS Parent should connect protein_id to old_locus_tag on the gene row."""
    gff = write_gff_lines(
        tmp_path / "genome.gff",
        [
            (
                "NC_003552.1\tRefSeq\tgene\t59888\t61279\t.\t+\t.\t"
                "ID=gene-MA_RS00255;Dbxref=GeneID:1471942;Name=MA_RS00255;"
                "gbkey=Gene;gene_biotype=protein_coding;locus_tag=MA_RS00255;"
                "old_locus_tag=MA0050%2CMA_0050"
            ),
            (
                "NC_003552.1\tProtein Homology\tCDS\t59888\t61279\t.\t+\t0\t"
                "ID=cds-WP_011020109.1;Parent=gene-MA_RS00255;"
                "Dbxref=GenBank:WP_011020109.1,GeneID:1471942;"
                "Name=WP_011020109.1;gbkey=CDS;locus_tag=MA_RS00255;"
                "product=ATP-binding protein;protein_id=WP_011020109.1;"
                "transl_table=11"
            ),
        ],
    )

    mapping = load_gff_locus_map(gff)

    assert mapping["WP_011020109.1"] == "MA_0050"
    assert mapping["WP_011020109"] == "MA_0050"


def test_gff_url_encoded_comma_is_decoded() -> None:
    """URL-encoded commas in GFF attributes should become separate values."""
    attributes = parse_gff_attributes("old_locus_tag=MA0050%2CMA_0050")

    assert attributes["old_locus_tag"] == ["MA0050", "MA_0050"]


def test_gff_compact_ma_tag_normalizes_to_old_locus_tag(tmp_path: Path) -> None:
    """Compact MA0050 tags should normalize to MA_0050."""
    gff = write_gff(
        tmp_path / "genome.gff",
        "protein_id=WP_011020109.1;old_locus_tag=MA0050",
    )

    mapping = load_gff_locus_map(gff)

    assert mapping["WP_011020109.1"] == "MA_0050"


def test_gff_parent_old_locus_tag_overrides_cds_locus_tag(
    tmp_path: Path,
) -> None:
    """Parent gene old_locus_tag should override a CDS locus_tag."""
    gff = write_gff_lines(
        tmp_path / "genome.gff",
        [
            (
                "seqid\tRefSeq\tgene\t1\t100\t.\t+\t.\t"
                "ID=gene-MA_RS00255;old_locus_tag=MA0050%2CMA_0050"
            ),
            (
                "seqid\tRefSeq\tCDS\t1\t100\t.\t+\t0\t"
                "protein_id=WP_011020109.1;Parent=gene-MA_RS00255;"
                "locus_tag=MA_RS00255"
            ),
        ],
    )

    mapping = load_gff_locus_map(gff)

    assert mapping["WP_011020109.1"] == "MA_0050"


def test_gff_direct_cds_old_locus_tag_still_works(tmp_path: Path) -> None:
    """Direct CDS old_locus_tag should continue to map protein IDs."""
    gff = write_gff(
        tmp_path / "genome.gff",
        "protein_id=WP_011020109.1;old_locus_tag=MA_0050;locus_tag=MA_RS00255",
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


def test_gff_annotation_matches_unversioned_record_ids() -> None:
    """Versioned GFF IDs should also match records without the version suffix."""
    records = {
        "WP_011020109": ProteinRecord(protein_id="WP_011020109"),
    }

    updated = annotate_records_with_gff_locus_tags(
        records,
        {
            "WP_011020109.1": "MA_0050",
            "WP_011020109": "MA_0050",
        },
    )

    assert updated == 1
    assert records["WP_011020109"].old_locus_tag == "MA_0050"


def test_gff_annotation_ignores_directory_source_description_suffix() -> None:
    """A FASTA description suffix should not prevent protein ID matching."""
    records = {
        "WP_011020109.1": ProteinRecord(
            protein_id="WP_011020109.1",
            description="WP_011020109.1 [source=Sulfolobus_solfataricus]",
        ),
    }

    updated = annotate_records_with_gff_locus_tags(
        records,
        {"WP_011020109.1": "MA_0050"},
    )

    assert updated == 1
    assert records["WP_011020109.1"].old_locus_tag == "MA_0050"


def test_normalize_protein_id_supports_common_gff_prefixes() -> None:
    """Common GFF prefixes should normalize to the protein accession."""
    assert normalize_protein_id("WP_011020109.1") == "WP_011020109.1"
    assert normalize_protein_id("cds-WP_011020109.1") == "WP_011020109.1"
    assert normalize_protein_id("Genbank:WP_011020109.1") == "WP_011020109.1"
    assert normalize_protein_id("GenBank:WP_011020109.1") == "WP_011020109.1"
    assert normalize_protein_id("RefSeq:WP_011020109.1") == "WP_011020109.1"
    assert normalize_protein_id("NCBI_GP:WP_011020109.1") == "WP_011020109.1"
