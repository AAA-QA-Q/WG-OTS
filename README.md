<p align="center">
  <img src="doc/WG-OTS.png" alt="WG-OTS logo" width="320">
</p>

<h1 align="center">WG-OTS</h1>

<p align="center">
  <b>WholeGenome OffTarget Scanner</b><br>
  A Snakemake workflow for whole-genome off-target detection and annotation from WGS data.
</p>

**WG-OTS (WholeGenome OffTarget Scanner)** is a Snakemake workflow for whole-genome off-target detection and annotation from whole-genome sequencing (WGS) data.

This workflow is designed for paired control and injected samples, and aims to identify candidate off-target mutations by integrating read preprocessing, genome alignment, somatic-style variant calling, variant filtering, functional annotation, and recurrent-gene screening.

---

## Overview

WG-OTS takes the following inputs:

- reference genome FASTA
- genome annotation GFF
- control paired-end FASTQ files
- injected paired-end FASTQ files

In the current workflow design:

- **ctrl** is treated as the **normal** sample
- **inj** is treated as the **tumor** sample

DeepSomatic is used in tumor-normal mode to detect candidate variants, after which the workflow keeps PASS indels, filters by minimum depth, annotates variants with snpEff, extracts frameshift variants, and finally retains recurrently hit genes.

---

## Workflow logic

The workflow consists of the following major steps:

### 1. Genome preparation
The input reference FASTA and GFF annotation are processed with `gffread` to generate:

- `genome.cds`
- `genome.pep`
- `genome.gtf`

These files are used in downstream analyses, including custom snpEff database construction.

### 2. Read preprocessing
Raw FASTQ files from control and injected samples are processed by `fastp` to perform quality control and read filtering.


### 3. Read alignment
Filtered reads are aligned to the reference genome using `bwa mem`, converted and sorted with `samtools sort`, and then read-group tags are added with `samtools addreplacerg`.


### 4. DeepSomatic calling
The workflow runs **DeepSomatic** in tumor-normal mode:

- `ctrl` -> normal
- `inj` -> tumor

This step generates the primary somatic-style variant callset in VCF/gVCF format.


### 5. bcftools filtering
The DeepSomatic VCF is filtered in two stages:

1. keep only **PASS** records and **indels**
2. retain only variants with `DP >= MIN_DP`


### 6. snpEff database preparation
A custom snpEff annotation database is built from:

- `genome.gtf`
- `ref.fa`
- `genome.cds`
- `genome.pep`


### 7. snpEff annotation
The filtered indel VCF is annotated using the custom snpEff database.

### 8. Frameshift extraction
Annotated variants containing `frameshift_variant` are extracted.

### 9. Recurrent-gene filtering
A custom Python script is used to retain only genes that appear in multiple independent variant lines, generating the final candidate result file.

This filtering step is motivated by the fact that **Cas9-induced cleavage outcomes are typically heterogeneous**, and true editing events in the same target gene may be represented by multiple nearby indel records rather than a single perfectly uniform mutation pattern.

Therefore, genes supported by only one variant line are considered lower-confidence candidates, because such signals may arise from **sample-specific background variation, sequencing noise, or unrelated stochastic events**. To improve specificity, the workflow removes genes represented by only a single variant record and retains genes with recurrent variant support.

---

## Project structure

```text
WG-OTS/
├── Snakefile
├── config.yaml
├── env.yaml
├── example/
├── script/
│   └── filter_recurrent_genes_from_snpeff_vcf.py
```

---

## Requirements

WG-OTS requires:

- Snakemake
- Conda or Mamba
- Docker
- the DeepSomatic Docker image prepared locally before running

---

## Important note about DeepSomatic

⚠️WG-OTS does **not** automatically install the DeepSomatic Docker image.

⚠️Before running the workflow, users must prepare the image manually, for example by pulling it locally：

```bash
docker pull google/deepsomatic:1.10.0
```

---

## Input configuration

The workflow is configured through `config.yaml`.

Example:

```yaml
# input
REF: "example/genome.fa"
GFF: "example/genomic.gff"

ctrl_fastq1: "example/ctrl.1.fq.gz"
ctrl_fastq2: "example/ctrl.2.fq.gz"
inj_fastq1: "example/inj.1.fq.gz"
inj_fastq2: "example/inj.2.fq.gz"

# bcftools filter
# MIN_DP: minimum depth threshold for variant filtering.
# Recommended setting: sequencing depth × 3/4.
MIN_DP: 8

# snpEff
# SNPEFF_MEM_GB: Java heap memory (in GB) used for building the snpEff annotation database.
# Increase this value if snpEff database construction runs out of memory.
SNPEFF_MEM_GB: 200
```

---

## Run the workflow

```bash
snakemake --cores 20 --software-deployment-method conda
```

If interrupted outputs need to be rebuilt:

```bash
snakemake --cores 20 --software-deployment-method conda --rerun-incomplete
```

---

## Main outputs

Key outputs include:

- `04_deepsomatic/inj_vs_ctrl.deepsomatic.vcf.gz`  
  Primary variant callset generated by DeepSomatic in tumor-normal mode.

- `05_bcftools/inj_vs_ctrl.pass.indel.vcf.gz`  
  DeepSomatic VCF after retaining only `PASS` records and indels.

- `05_bcftools/inj_vs_ctrl.pass.indel.MIN_DPfilter.vcf.gz`  
  PASS indel callset after additional filtering by the minimum depth threshold (`MIN_DP`).

- `06_snpeff/inj_vs_ctrl.pass.indel.MIN_DPfilter.ann.vcf.gz`  
  Depth-filtered indel VCF annotated with the custom snpEff database.

- `06_snpeff/inj_vs_ctrl.pass.indel.MIN_DPfilter.ann.frameshift_variant.vcf.gz`  
  Annotated VCF containing only variants labeled as `frameshift_variant`.

- `06_snpeff/inj_vs_ctrl.final.cas9.results`  
  Final high-confidence candidate output after recurrent-gene filtering, retaining genes supported by multiple variant records.

---

## Suggested citations

If you use WG-OTS in a publication, please cite the software used in this workflow as appropriate.

### Workflow engine
- Snakemake [Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33.](https://f1000research.com/articles/10-33/v1)


### Read preprocessing
- gffread [Pertea G, Pertea M. GFF utilities: GffRead and GffCompare[J]. F1000Research, 2020, 9: ISCB Comm J-304.] (https://f1000research.com/articles/9-304/v1)
- fastp [Shifu Chen. fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and preprocessing. iMeta 4.5 (2025): e70078 https://doi.org/10.1002/imt2.70078](https://onlinelibrary.wiley.com/doi/10.1002/imt2.70078)

### Read alignment
- bwa [Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].](https://arxiv.org/abs/1303.3997)

### SAM/BAM processing
- samtools [Danecek P, Bonfield J K, Liddle J, et al. Twelve years of SAMtools and BCFtools[J]. Gigascience, 2021, 10(2): giab008.](https://doi.org/10.1093/gigascience/giab008)

### BCF/VCF processing
- bcftools [Danecek P, Bonfield J K, Liddle J, et al. Twelve years of SAMtools and BCFtools[J]. Gigascience, 2021, 10(2): giab008.](https://doi.org/10.1093/gigascience/giab008)

### Variant annotation
- Snpeff ["A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.", Cingolani P, Platts A, Wang le L, Coon M, Nguyen T, Wang L, Land SJ, Lu X, Ruden DM. Fly (Austin). 2012 Apr-Jun;6(2):80-92. PMID: 22728672](https://www.tandfonline.com/doi/full/10.4161/fly.19695)

### Somatic variant calling
- DeepSomatic [DeepSomatic: Accurate somatic small variant discovery for multiple sequencing technologies](https://www.biorxiv.org/content/10.1101/2024.08.16.608331v1)

---

## License

MIT License
