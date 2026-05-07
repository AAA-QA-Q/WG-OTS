#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import gzip
import os
import sys
from collections import Counter


def open_maybe_gzip(path, mode="rt"):
    """
    根据文件后缀自动判断是否为 gzip 文件。
    .gz 用 gzip.open，其它文件用普通 open。
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding="utf-8")


def extract_gene_names_from_info(info_field):
    """
    从 INFO 字段中提取 snpEff ANN 里的 Gene_Name。

    snpEff ANN 标准格式大致为：
    Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|...

    这里取第 4 列（下标 3）作为 Gene_Name。

    返回值：
    - 一个 set，表示该行记录涉及到的所有唯一基因名
    """
    gene_names = set()

    # INFO 以 ; 分隔，先找到 ANN=...
    ann_value = None
    for item in info_field.split(";"):
        if item.startswith("ANN="):
            ann_value = item[4:]
            break

    if ann_value is None or ann_value == "":
        return gene_names

    # 多个注释之间用逗号分隔
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
    第一遍读取文件：
    统计每个基因出现在多少条 VCF 记录中。

    注意：
    - 同一条记录里，同一个基因即使出现多次（不同转录本），也只算 1 次。
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
    第二遍读取文件：
    只保留那些“至少有一个基因在全文件中出现次数 >= min_records_per_gene”的记录。

    输出要求：
    - 不保留 ## 开头的元信息行
    - 只保留一行 #CHROM ... 表头
    - 输出为普通文本文件
    """
    gene_counter = count_genes_by_record(vcf_path)

    kept_records = 0
    removed_records = 0

    with open_maybe_gzip(vcf_path, "rt") as fin, open(output_path, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                continue

            # 只保留最后那一行列名表头
            if line.startswith("#CHROM"):
                fout.write(line)
                continue

            # 其它 ## 元信息行全部跳过
            if line.startswith("##"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            info_field = cols[7]
            genes_in_this_record = extract_gene_names_from_info(info_field)

            # 如果这一行没有解析到基因，默认删除
            if not genes_in_this_record:
                removed_records += 1
                continue

            # 只要这一行中任意一个基因在全文件中出现次数 >= 2，就保留
            keep = any(gene_counter[gene] >= min_records_per_gene for gene in genes_in_this_record)

            if keep:
                fout.write(line)
                kept_records += 1
            else:
                removed_records += 1

    return kept_records, removed_records, gene_counter


def parse_args():
    """
    解析命令行参数。
    """
    parser = argparse.ArgumentParser(
        description=(
            "从 snpEff 注释后的 VCF/VCF.GZ 中，删除只出现 1 次的基因对应记录，"
            "并只保留 #CHROM 表头。"
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="输入文件，支持 .vcf 或 .vcf.gz"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="输出文件，例如 inj_vs_ctrl.final.cas9.results"
    )
    parser.add_argument(
        "--min-records-per-gene",
        type=int,
        default=2,
        help="一个基因至少出现在多少条记录里才保留，默认 2"
    )
    return parser.parse_args()


def main():
    """
    主函数：
    - 检查输入文件
    - 执行过滤
    - 输出统计信息
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
