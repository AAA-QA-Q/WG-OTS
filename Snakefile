import os

configfile: "config.yaml"

REF = config["REF"]
GFF = config["GFF"]

CTRL_R1 = config["ctrl_fastq1"]
CTRL_R2 = config["ctrl_fastq2"]
INJ_R1  = config["inj_fastq1"]
INJ_R2  = config["inj_fastq2"]

MIN_DP = int(config["MIN_DP"])
SNPEFF_MEM_GB = int(config["SNPEFF_MEM_GB"])

GENOME_DIR = "01_genome"
REF_LINK = f"{GENOME_DIR}/ref.fa"
GFF_LINK = f"{GENOME_DIR}/genome.gff"
GENOME_CDS = f"{GENOME_DIR}/genome.cds"
GENOME_PEP = f"{GENOME_DIR}/genome.pep"
GENOME_GTF = f"{GENOME_DIR}/genome.gtf"

# samtools faidx for the real reference used by DeepSomatic
REF_FAI = REF + ".fai"

FASTP_DIR = "02_fastp"

CLEAN_CTRL_R1 = f"{FASTP_DIR}/ctrl_1.clean.fq.gz"
CLEAN_CTRL_R2 = f"{FASTP_DIR}/ctrl_2.clean.fq.gz"
CLEAN_INJ_R1  = f"{FASTP_DIR}/inj_1.clean.fq.gz"
CLEAN_INJ_R2  = f"{FASTP_DIR}/inj_2.clean.fq.gz"

CTRL_HTML = f"{FASTP_DIR}/ctrl.fastp.html"
CTRL_JSON = f"{FASTP_DIR}/ctrl.fastp.json"
INJ_HTML  = f"{FASTP_DIR}/inj.fastp.html"
INJ_JSON  = f"{FASTP_DIR}/inj.fastp.json"

SAMPLES = {
    "ctrl": {
        "r1": CTRL_R1,
        "r2": CTRL_R2
    },
    "inj": {
        "r1": INJ_R1,
        "r2": INJ_R2
    }
}

BWA_DIR = "03_bwa"

SORTED_CTRL_BAM = f"{BWA_DIR}/ctrl.sorted.bam"
SORTED_INJ_BAM  = f"{BWA_DIR}/inj.sorted.bam"

RG_CTRL_BAM = f"{BWA_DIR}/ctrl.sorted.rg.bam"
RG_INJ_BAM  = f"{BWA_DIR}/inj.sorted.rg.bam"

RG_CTRL_BAI = f"{RG_CTRL_BAM}.bai"
RG_INJ_BAI  = f"{RG_INJ_BAM}.bai"

DEEPSOMATIC_IMAGE = "google/deepsomatic:1.10.0"
DEEPSOMATIC_DIR = "04_deepsomatic"

DEEPSOMATIC_VCF = f"{DEEPSOMATIC_DIR}/inj_vs_ctrl.deepsomatic.vcf.gz"
DEEPSOMATIC_GVCF = f"{DEEPSOMATIC_DIR}/inj_vs_ctrl.deepsomatic.g.vcf.gz"
DEEPSOMATIC_LOG_DONE = f"{DEEPSOMATIC_DIR}/logs/.done"
DEEPSOMATIC_INTERMEDIATE_DONE = f"{DEEPSOMATIC_DIR}/intermediate_results_dir/.done"

BCFTOOLS_DIR = "05_bcftools"

PASS_INDEL_VCF = f"{BCFTOOLS_DIR}/inj_vs_ctrl.pass.indel.vcf.gz"
PASS_INDEL_TBI = f"{PASS_INDEL_VCF}.tbi"

PASS_INDEL_DP_VCF = f"{BCFTOOLS_DIR}/inj_vs_ctrl.pass.indel.MIN_DPfilter.vcf.gz"
PASS_INDEL_DP_TBI = f"{PASS_INDEL_DP_VCF}.tbi"

SNPEFF_DIR = "06_snpeff"
SNPEFF_CONFIG = f"{SNPEFF_DIR}/snpEff.config"
SNPEFF_DATA_DIR = f"{SNPEFF_DIR}/data"
SNPEFF_GENOME = "ann_database"
SNPEFF_GENOME_DIR = f"{SNPEFF_DATA_DIR}/{SNPEFF_GENOME}"

SNPEFF_GENES_GTF = f"{SNPEFF_GENOME_DIR}/genes.gtf"
SNPEFF_SEQUENCES = f"{SNPEFF_GENOME_DIR}/sequences.fa"
SNPEFF_CDS = f"{SNPEFF_GENOME_DIR}/cds.fa"
SNPEFF_PROTEIN = f"{SNPEFF_GENOME_DIR}/protein.fa"

SNPEFF_BUILD_DONE = f"{SNPEFF_GENOME_DIR}/.build.done"

SNPEFF_ANN_VCF = f"{SNPEFF_DIR}/inj_vs_ctrl.pass.indel.MIN_DPfilter.ann.vcf.gz"
SNPEFF_ANN_TBI = f"{SNPEFF_ANN_VCF}.tbi"

FRAMESHIFT_VCF = f"{SNPEFF_DIR}/inj_vs_ctrl.pass.indel.MIN_DPfilter.ann.frameshift_variant.vcf.gz"
FRAMESHIFT_TBI = f"{FRAMESHIFT_VCF}.tbi"

FINAL_CAS9_RESULTS = f"{SNPEFF_DIR}/inj_vs_ctrl.final.cas9.results"


rule all:
    input:
        REF_LINK,
        GFF_LINK,
        GENOME_CDS,
        GENOME_PEP,
        GENOME_GTF,
        REF_FAI,
        CLEAN_CTRL_R1,
        CLEAN_CTRL_R2,
        CLEAN_INJ_R1,
        CLEAN_INJ_R2,
        CTRL_HTML,
        CTRL_JSON,
        INJ_HTML,
        INJ_JSON,
        RG_CTRL_BAM,
        RG_INJ_BAM,
        RG_CTRL_BAI,
        RG_INJ_BAI,
        DEEPSOMATIC_VCF,
        DEEPSOMATIC_GVCF,
        DEEPSOMATIC_LOG_DONE,
        DEEPSOMATIC_INTERMEDIATE_DONE,
        PASS_INDEL_VCF,
        PASS_INDEL_TBI,
        PASS_INDEL_DP_VCF,
        PASS_INDEL_DP_TBI,
        SNPEFF_CONFIG,
        SNPEFF_GENES_GTF,
        SNPEFF_SEQUENCES,
        SNPEFF_CDS,
        SNPEFF_PROTEIN,
        SNPEFF_BUILD_DONE,
        SNPEFF_ANN_VCF,
        SNPEFF_ANN_TBI,
        FRAMESHIFT_VCF,
        FRAMESHIFT_TBI,
        FINAL_CAS9_RESULTS


rule prepare_genome_links:
    input:
        ref=REF,
        gff=GFF
    output:
        ref=REF_LINK,
        gff=GFF_LINK
    shell:
        r"""
        mkdir -p {GENOME_DIR}
        ln -sf "$(realpath {input.ref})" {output.ref}
        ln -sf "$(realpath {input.gff})" {output.gff}
        """


rule gffread_extract:
    input:
        ref=REF_LINK,
        gff=GFF_LINK
    output:
        cds=GENOME_CDS,
        pep=GENOME_PEP,
        gtf=GENOME_GTF
    conda:
        "env.yaml"
    shell:
        r"""
        gffread {input.gff} \
            -g {input.ref} \
            -x {output.cds} \
            -y {output.pep} \
            -T \
            -o {output.gtf}
        """


rule faidx_ref:
    input:
        ref=REF
    output:
        fai=REF_FAI
    conda:
        "env.yaml"
    shell:
        r"""
        samtools faidx {input.ref}
        """


rule fastp_clean:
    input:
        r1=lambda wc: SAMPLES[wc.sample]["r1"],
        r2=lambda wc: SAMPLES[wc.sample]["r2"]
    output:
        r1=f"{FASTP_DIR}/{{sample}}_1.clean.fq.gz",
        r2=f"{FASTP_DIR}/{{sample}}_2.clean.fq.gz",
        html=f"{FASTP_DIR}/{{sample}}.fastp.html",
        json=f"{FASTP_DIR}/{{sample}}.fastp.json"
    threads: 20
    conda:
        "env.yaml"
    shell:
        r"""
        mkdir -p {FASTP_DIR}
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            -h {output.html} \
            -j {output.json} \
            -w {threads}
        """


rule bwa_index:
    input:
        REF_LINK
    output:
        amb=REF_LINK + ".amb",
        ann=REF_LINK + ".ann",
        bwt=REF_LINK + ".bwt",
        pac=REF_LINK + ".pac",
        sa=REF_LINK + ".sa"
    conda:
        "env.yaml"
    shell:
        r"""
        bwa index {input}
        """


rule bwa_mem_sort:
    input:
        ref=REF_LINK,
        amb=REF_LINK + ".amb",
        ann=REF_LINK + ".ann",
        bwt=REF_LINK + ".bwt",
        pac=REF_LINK + ".pac",
        sa=REF_LINK + ".sa",
        r1=f"{FASTP_DIR}/{{sample}}_1.clean.fq.gz",
        r2=f"{FASTP_DIR}/{{sample}}_2.clean.fq.gz"
    output:
        bam=temp(f"{BWA_DIR}/{{sample}}.sorted.bam")
    threads: workflow.cores
    conda:
        "env.yaml"
    params:
        bwa_threads=lambda wc, threads: max(1, threads // 2),
        sort_threads=lambda wc, threads: max(1, threads - (threads // 2))
    shell:
        r"""
        mkdir -p {BWA_DIR}
        bwa mem -t {params.bwa_threads} {input.ref} {input.r1} {input.r2} \
          | samtools sort -@ {params.sort_threads} -m 2G -T {BWA_DIR}/{wildcards.sample}.tmp -o {output.bam}
        """


rule add_read_group:
    input:
        bam=f"{BWA_DIR}/{{sample}}.sorted.bam"
    output:
        bam=f"{BWA_DIR}/{{sample}}.sorted.rg.bam",
        bai=f"{BWA_DIR}/{{sample}}.sorted.rg.bam.bai"
    threads: 10
    conda:
        "env.yaml"
    params:
        rgid=lambda wc: wc.sample,
        rgsm=lambda wc: wc.sample
    shell:
        r"""
        samtools addreplacerg \
          -r ID:{params.rgid} \
          -r SM:{params.rgsm} \
          -o {output.bam} \
          {input.bam}

        samtools index -@ {threads} {output.bam}
        """


rule run_deepsomatic:
    input:
        normal_bam=RG_CTRL_BAM,
        normal_bai=RG_CTRL_BAI,
        tumor_bam=RG_INJ_BAM,
        tumor_bai=RG_INJ_BAI,
        ref=REF,
        fai=REF_FAI
    output:
        vcf=DEEPSOMATIC_VCF,
        gvcf=DEEPSOMATIC_GVCF,
        log_done=DEEPSOMATIC_LOG_DONE,
        intermediate_done=DEEPSOMATIC_INTERMEDIATE_DONE
    threads: workflow.cores
    params:
        image=DEEPSOMATIC_IMAGE,
        bam_dir=lambda wc, input: os.path.abspath(os.path.dirname(input.normal_bam)),
        ref_dir=lambda wc: os.path.dirname(os.path.realpath(REF)),
        ref_base=lambda wc: os.path.basename(os.path.realpath(REF)),
        out_dir=lambda wc: os.path.abspath(DEEPSOMATIC_DIR),
        normal_bam_base=lambda wc, input: os.path.basename(input.normal_bam),
        tumor_bam_base=lambda wc, input: os.path.basename(input.tumor_bam),
        vcf_base=lambda wc, output: os.path.basename(output.vcf),
        gvcf_base=lambda wc, output: os.path.basename(output.gvcf)
    shell:
        r"""
        mkdir -p {params.out_dir}
        mkdir -p {params.out_dir}/logs
        mkdir -p {params.out_dir}/intermediate_results_dir

        docker image inspect {params.image} >/dev/null 2>&1 || {{
            echo "ERROR: Docker image {params.image} not found locally."
            echo "Please prepare it first with docker pull or docker load."
            exit 1
        }}

        docker run --rm \
          -v {params.bam_dir}:/bam \
          -v {params.ref_dir}:/ref \
          -v {params.out_dir}:/out \
          {params.image} \
          run_deepsomatic \
          --model_type=WGS \
          --ref=/ref/{params.ref_base} \
          --reads_normal=/bam/{params.normal_bam_base} \
          --reads_tumor=/bam/{params.tumor_bam_base} \
          --output_vcf=/out/{params.vcf_base} \
          --output_gvcf=/out/{params.gvcf_base} \
          --sample_name_normal=ctrl \
          --sample_name_tumor=inj \
          --num_shards={threads} \
          --logging_dir=/out/logs \
          --intermediate_results_dir=/out/intermediate_results_dir

        touch {output.log_done}
        touch {output.intermediate_done}
        """


rule bcftools_pass_indel:
    input:
        vcf=DEEPSOMATIC_VCF
    output:
        vcf=PASS_INDEL_VCF,
        tbi=PASS_INDEL_TBI
    conda:
        "env.yaml"
    shell:
        r"""
        mkdir -p {BCFTOOLS_DIR}
        bcftools view \
          -f PASS \
          -v indels \
          -Oz \
          -o {output.vcf} \
          {input.vcf}

        bcftools index -t {output.vcf}
        """


rule bcftools_min_dp_filter:
    input:
        vcf=PASS_INDEL_VCF,
        tbi=PASS_INDEL_TBI
    output:
        vcf=PASS_INDEL_DP_VCF,
        tbi=PASS_INDEL_DP_TBI
    conda:
        "env.yaml"
    params:
        min_dp=MIN_DP
    shell:
        r"""
        bcftools view \
          -i 'FORMAT/DP[0]>={params.min_dp}' \
          -Oz \
          -o {output.vcf} \
          {input.vcf}

        bcftools index -t {output.vcf}
        """


rule prepare_snpeff_inputs:
    input:
        gtf=GENOME_GTF,
        ref=REF_LINK,
        cds=GENOME_CDS,
        protein=GENOME_PEP
    output:
        config=SNPEFF_CONFIG,
        genes_gtf=SNPEFF_GENES_GTF,
        sequences=SNPEFF_SEQUENCES,
        cds=SNPEFF_CDS,
        protein=SNPEFF_PROTEIN
    shell:
        r"""
        mkdir -p {SNPEFF_GENOME_DIR}

        ln -sf "$(realpath {input.gtf})" {output.genes_gtf}
        ln -sf "$(realpath {input.ref})" {output.sequences}
        ln -sf "$(realpath {input.cds})" {output.cds}
        ln -sf "$(realpath {input.protein})" {output.protein}

        cat > {output.config} << EOF
data.dir = ./data

{SNPEFF_GENOME}.genome : {SNPEFF_GENOME}
EOF
        """


rule build_snpeff_db:
    input:
        config=SNPEFF_CONFIG,
        genes_gtf=SNPEFF_GENES_GTF,
        sequences=SNPEFF_SEQUENCES,
        cds=SNPEFF_CDS,
        protein=SNPEFF_PROTEIN
    output:
        done=SNPEFF_BUILD_DONE
    conda:
        "env.yaml"
    params:
        genome=SNPEFF_GENOME,
        xmx=lambda wc: f"{SNPEFF_MEM_GB}g"
    shell:
        r"""
        SNPEFF_JAR=$(find "$CONDA_PREFIX" -name "snpEff.jar" | head -n 1)

        if [ -z "$SNPEFF_JAR" ]; then
            echo "ERROR: snpEff.jar not found under $CONDA_PREFIX"
            exit 1
        fi

        java -Xmx{params.xmx} -jar "$SNPEFF_JAR" build \
          -c {input.config} \
          -gtf22 \
          -v {params.genome}

        touch {output.done}
        """


rule snpeff_annotate_vcf:
    input:
        db_done=SNPEFF_BUILD_DONE,
        vcf=PASS_INDEL_DP_VCF,
        tbi=PASS_INDEL_DP_TBI,
        config=SNPEFF_CONFIG
    output:
        vcf=SNPEFF_ANN_VCF,
        tbi=SNPEFF_ANN_TBI
    conda:
        "env.yaml"
    params:
        genome=SNPEFF_GENOME,
        xmx=lambda wc: f"{SNPEFF_MEM_GB}g"
    shell:
        r"""
        SNPEFF_JAR=$(find "$CONDA_PREFIX" -name "snpEff.jar" | head -n 1)

        if [ -z "$SNPEFF_JAR" ]; then
            echo "ERROR: snpEff.jar not found under $CONDA_PREFIX"
            exit 1
        fi

        java -Xmx{params.xmx} -jar "$SNPEFF_JAR" ann \
          -c {input.config} \
          -noStats \
          {params.genome} \
          {input.vcf} \
          | bcftools view -Oz -o {output.vcf}

        bcftools index -t {output.vcf}
        """


rule extract_frameshift_variants:
    input:
        vcf=SNPEFF_ANN_VCF,
        tbi=SNPEFF_ANN_TBI
    output:
        vcf=FRAMESHIFT_VCF,
        tbi=FRAMESHIFT_TBI
    conda:
        "env.yaml"
    shell:
        r"""
        bcftools view \
          -i 'ANN~"frameshift_variant"' \
          -Oz \
          -o {output.vcf} \
          {input.vcf}

        bcftools index -t {output.vcf}
        """


rule filter_recurrent_genes:
    input:
        vcf=FRAMESHIFT_VCF,
        tbi=FRAMESHIFT_TBI
    output:
        FINAL_CAS9_RESULTS
    conda:
        "env.yaml"
    shell:
        r"""
        python3 script/final_result.py \
          -i {input.vcf} \
          -o {output}
        """
