#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import os
import sys
from collections import Counter


def open_maybe_gzip(path, mode="rt"):
    """
    Automatically open a plain text file or a gzip-compressed file
    based on the file extension.

    Files ending with .gz are opened with gzip.open;
    all other files are opened with the regular open function.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def extract_gene_names_from_info(info_field):
    """
    Extract Gene_Name values from the snpEff ANN field in the INFO column.

    The standard snpEff ANN format is approximately:
    Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|...

    Here, the 4th field (index 3) is used as Gene_Name.

    Returns:
        A set containing all unique gene names associated with this record.
    """
    gene_names = set()

    # INFO fields are separated by ';', so first locate ANN=...
    ann_value = None
    for item in info_field.split(";"):
        if item.startswith("ANN="):
            ann_value = item[4:]
            break

    if ann_value is None or ann_value == "":
        return gene_names

    # Multiple annotations are separated by commas
    ann_entries = ann_value.split(",")

    for entry in ann_entries:
        fields = entry.split("|")
        if len(fields) > 3:
            gene_name = fields[3].strip()
            if gene_name != "":
                gene_names.add(gene_name)

    return gene_names


def count_genes_by_record(vcf_path):
    """
    First pass through the file:
    Count how many VCF records each gene appears in.

    Note:
        If the same gene appears multiple times within one record
        (for example, due to multiple transcripts), it is counted only once.
    """
    gene_counter = Counter()

    with open_maybe_gzip(vcf_path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            info_field = cols[7]
            genes_in_this_record = extract_gene_names_from_info(info_field)

            for gene in genes_in_this_record:
                gene_counter[gene] += 1

    return gene_counter


def filter_vcf_by_gene_count(vcf_path, output_path, min_records_per_gene=2):
    """
    Second pass through the file:
    Retain only records where at least one gene appears in the full file
    at least min_records_per_gene times.

    Output rules:
        - Exclude all metadata lines starting with ##
        - Retain only the #CHROM header line
        - Write output as a plain text file
    """
    gene_counter = count_genes_by_record(vcf_path)

    kept_records = 0
    removed_records = 0

    with open_maybe_gzip(vcf_path, "rt") as fin, open(output_path, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue

            # Keep only the final column header line
            if line.startswith("#CHROM"):
                fout.write(line)
                continue

            # Skip all other metadata lines starting with ##
            if line.startswith("##"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            info_field = cols[7]
            genes_in_this_record = extract_gene_names_from_info(info_field)

            # If no gene can be parsed from this record, remove it by default
            if not genes_in_this_record:
                removed_records += 1
                continue

            # Keep the record if any gene in this line appears at least
            # min_records_per_gene times in the full file
            keep = any(gene_counter[gene] >= min_records_per_gene for gene in genes_in_this_record)

            if keep:
                fout.write(line)
                kept_records += 1
            else:
                removed_records += 1

    return kept_records, removed_records, gene_counter


def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Remove records associated with genes that appear only once in an "
            "snpEff-annotated VCF/VCF.GZ file, while retaining only the #CHROM header."
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input file, supporting .vcf or .vcf.gz"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file, for example: inj_vs_ctrl.final.cas9.results"
    )
    parser.add_argument(
        "--min-records-per-gene",
        type=int,
        default=2,
        help="Minimum number of records required for a gene to be retained (default: 2)"
    )
    return parser.parse_args()


def main():
    """
    Main function:
        - Check whether the input file exists
        - Run the filtering procedure
        - Print summary statistics
    """
    args = parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    kept_records, removed_records, gene_counter = filter_vcf_by_gene_count(
        vcf_path=args.input,
        output_path=args.output,
        min_records_per_gene=args.min_records_per_gene
    )

    recurrent_gene_count = sum(1 for gene, n in gene_counter.items() if n >= args.min_records_per_gene)

    print("Done.")
    print(f"Input file:  {args.input}")
    print(f"Output file: {args.output}")
    print(f"Genes with >= {args.min_records_per_gene} records: {recurrent_gene_count}")
    print(f"Kept records:    {kept_records}")
    print(f"Removed records: {removed_records}")


if __name__ == "__main__":
    main()
