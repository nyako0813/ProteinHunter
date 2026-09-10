#!/usr/bin/env python3
"""Convert an NCBI GenBank flat file (.gb/.gbk, e.g. from
`efetch.fcgi?db=nuccore&id=<acc>&rettype=gbwithparts&retmode=text`) into the
legacy NCBI Genome Projects .ptt/.rnt/.fna trio that Rockhopper's `-g` option
expects.

Background: Rockhopper (https://cs.wellesley.edu/~btjaden/Rockhopper/) does
NOT accept GFF3. It wants a directory containing exactly one *.fna (genome
sequence, FASTA), one *.ptt (protein/CDS table) and one *.rnt (RNA table),
using the classic NCBI PTT/RNT column format:

    <description> - 1..<seq_len>
    <n> proteins
    Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct

- Location: "<start>..<end>", 1-based inclusive (GenBank convention already;
  Biopython feature.location is 0-based half-open internally so we add 1 to
  the start).
- Strand: '+' or '-'.
- Length: length in AMINO ACIDS for .ptt (i.e. (nt length / 3) - 1, excluding
  the stop codon) and in NUCLEOTIDES for .rnt.
- PID: protein_id (CDS) / GI if present, else the locus_tag as a fallback so
  the column is never empty (Rockhopper does not appear to use this field's
  content, but the column must be present and non-empty).
- Gene: /gene qualifier, else '-'.
- Synonym: /locus_tag qualifier (this is the identifier that round-trips back
  to MA_#### loci via old_locus_tag in the source GFF3).
- Code/COG: not available from a plain GenBank record -> '-'.
- Product: /product qualifier, else '-'.

Usage:
    python genbank_to_ptt_rnt.py input.gb output_dir/ [--replicon-name NAME]

Writes output_dir/<accession>.fna, .ptt, .rnt.
"""
import argparse
import os
import sys

from Bio import SeqIO


def cds_aa_length(feature, record_len):
    """Amino-acid length excluding the stop codon."""
    if "translation" in feature.qualifiers:
        return len(feature.qualifiers["translation"][0])
    nt_len = len(feature.location)
    return max(nt_len // 3 - 1, 0)


def get_qual(feature, name, default="-"):
    vals = feature.qualifiers.get(name)
    if not vals:
        return default
    return vals[0]


def location_str(feature):
    start = int(feature.location.start) + 1  # 0-based -> 1-based
    end = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-"
    return f"{start}..{end}", strand


def build_ptt_rows(record):
    rows = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        loc, strand = location_str(feat)
        length_aa = cds_aa_length(feat, len(record.seq))
        pid = get_qual(feat, "protein_id", get_qual(feat, "locus_tag", "-"))
        gene = get_qual(feat, "gene", "-")
        synonym = get_qual(feat, "locus_tag", "-")
        product = get_qual(feat, "product", "-")
        rows.append([loc, strand, str(length_aa), pid, gene, synonym, "-", "-", product])
    return rows


RNA_TYPES = {"tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA"}


def build_rnt_rows(record):
    rows = []
    for feat in record.features:
        if feat.type not in RNA_TYPES:
            continue
        loc, strand = location_str(feat)
        length_nt = len(feat.location)
        pid = get_qual(feat, "locus_tag", "-")
        gene = get_qual(feat, "gene", "-")
        synonym = get_qual(feat, "locus_tag", "-")
        product = get_qual(feat, "product", "-")
        rows.append([loc, strand, str(length_nt), pid, gene, synonym, "-", "-", product])
    return rows


def write_table(path, description, seq_len, rows, unit_label):
    with open(path, "w") as fh:
        fh.write(f"{description} - 1..{seq_len}\n")
        fh.write(f"{len(rows)} {unit_label}\n")
        fh.write("Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genbank_file")
    ap.add_argument("output_dir")
    ap.add_argument("--replicon-name", default=None, help="Override base filename (default: accession)")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    record = SeqIO.read(args.genbank_file, "genbank")
    base = args.replicon_name or record.id

    fna_path = os.path.join(args.output_dir, f"{base}.fna")
    ptt_path = os.path.join(args.output_dir, f"{base}.ptt")
    rnt_path = os.path.join(args.output_dir, f"{base}.rnt")

    description = record.description or record.id
    seq_len = len(record.seq)

    with open(fna_path, "w") as fh:
        fh.write(f">{record.id} {description}\n")
        seq = str(record.seq)
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")

    ptt_rows = build_ptt_rows(record)
    write_table(ptt_path, description, seq_len, ptt_rows, "proteins")

    rnt_rows = build_rnt_rows(record)
    write_table(rnt_path, description, seq_len, rnt_rows, "RNAs")

    print(f"Wrote {fna_path} ({seq_len} bp)", file=sys.stderr)
    print(f"Wrote {ptt_path} ({len(ptt_rows)} CDS)", file=sys.stderr)
    print(f"Wrote {rnt_path} ({len(rnt_rows)} RNA features)", file=sys.stderr)


if __name__ == "__main__":
    main()
