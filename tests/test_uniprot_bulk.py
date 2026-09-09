"""Tests for annotation/uniprot_bulk.py."""

from __future__ import annotations

import json
from pathlib import Path

from annotation.uniprot_bulk import UniProtDomainInfo, load_uniprot_bulk_domain_map


def _sample_payload() -> dict[str, object]:
    return {
        "results": [
            {
                "primaryAccession": "Q8TIN1",
                "organism": {"scientificName": "Methanosarcina acetivorans"},
                "proteinExistence": "4: Predicted",
                "genes": [{"orderedLocusNames": [{"value": "MA_4115"}]}],
                "uniProtKBCrossReferences": [
                    {
                        "database": "Pfam",
                        "id": "PF24167",
                        "properties": [{"key": "EntryName", "value": "DUF7411"}],
                    },
                    {
                        "database": "InterPro",
                        "id": "IPR055834",
                        "properties": [{"key": "EntryName", "value": "DUF7411"}],
                    },
                    {"database": "SUPFAM", "id": "SSF52402", "properties": []},
                    {"database": "GO", "id": "GO:0005524", "properties": []},
                ],
            },
            {
                "primaryAccession": "Q8TTR5",
                "organism": {"scientificName": "Methanosarcina acetivorans"},
                "proteinExistence": "4: Predicted",
                "genes": [{"orderedLocusNames": [{"value": "MA_0361"}]}],
                "uniProtKBCrossReferences": [
                    {"database": "Pfam", "id": "PF12724", "properties": []},
                ],
            },
            {
                # No orderedLocusNames -- must not appear in the resulting map.
                "primaryAccession": "Q0000X",
                "genes": [{"name": {"value": "unnamed"}}],
                "uniProtKBCrossReferences": [
                    {"database": "Pfam", "id": "PF99999", "properties": []},
                ],
            },
            {
                # Multiple orderedLocusNames -- registered under every tag.
                "primaryAccession": "Q0000Y",
                "genes": [
                    {
                        "orderedLocusNames": [
                            {"value": "MA_9001"},
                            {"value": "MA_9002"},
                        ]
                    }
                ],
                "uniProtKBCrossReferences": [
                    {"database": "InterPro", "id": "IPR000001", "properties": []},
                ],
            },
        ]
    }


def test_load_builds_old_locus_tag_map(tmp_path: Path) -> None:
    path = tmp_path / "uniprot_bulk.json"
    path.write_text(json.dumps(_sample_payload()), encoding="utf-8")

    domain_map = load_uniprot_bulk_domain_map(path)

    assert domain_map["MA_4115"] == UniProtDomainInfo(
        pfam=frozenset({"PF24167"}),
        interpro=frozenset({"IPR055834"}),
        supfam=frozenset({"SSF52402"}),
        protein_existence="4: Predicted",
    )
    assert domain_map["MA_0361"].pfam == frozenset({"PF12724"})


def test_load_ignores_non_domain_databases(tmp_path: Path) -> None:
    path = tmp_path / "uniprot_bulk.json"
    path.write_text(json.dumps(_sample_payload()), encoding="utf-8")

    domain_map = load_uniprot_bulk_domain_map(path)

    # GO:0005524 on MA_4115's entry must not leak into any of pfam/interpro/supfam.
    info = domain_map["MA_4115"]
    assert "GO:0005524" not in info.pfam
    assert "GO:0005524" not in info.interpro
    assert "GO:0005524" not in info.supfam


def test_load_skips_entries_without_old_locus_tag(tmp_path: Path) -> None:
    path = tmp_path / "uniprot_bulk.json"
    path.write_text(json.dumps(_sample_payload()), encoding="utf-8")

    domain_map = load_uniprot_bulk_domain_map(path)

    assert all("PF99999" not in info.pfam for info in domain_map.values())


def test_load_registers_multiple_old_locus_tags(tmp_path: Path) -> None:
    path = tmp_path / "uniprot_bulk.json"
    path.write_text(json.dumps(_sample_payload()), encoding="utf-8")

    domain_map = load_uniprot_bulk_domain_map(path)

    assert domain_map["MA_9001"].interpro == frozenset({"IPR000001"})
    assert domain_map["MA_9002"].interpro == frozenset({"IPR000001"})


def test_missing_file_returns_empty_dict_without_raising(tmp_path: Path) -> None:
    domain_map = load_uniprot_bulk_domain_map(tmp_path / "does_not_exist.json")
    assert domain_map == {}


def test_invalid_json_returns_empty_dict_without_raising(tmp_path: Path) -> None:
    path = tmp_path / "broken.json"
    path.write_text("{not valid json", encoding="utf-8")

    domain_map = load_uniprot_bulk_domain_map(path)

    assert domain_map == {}


def test_unexpected_top_level_shape_returns_empty_dict(tmp_path: Path) -> None:
    path = tmp_path / "wrong_shape.json"
    path.write_text(json.dumps(["not", "a", "mapping"]), encoding="utf-8")

    domain_map = load_uniprot_bulk_domain_map(path)

    assert domain_map == {}
