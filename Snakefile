### VARIANT CALLING PIPELINE ###

# author: Lorena Derezanin

# snakemake v.6.15.2.
# conda v.4.11.0
# mamba v.0.15.3 - more robust pkg manager than conda, handles snakemake releases and dependencies better
# run in conda env snek2
# snakemake --use-conda


# add unit tests 
# add config with default params - 2 CPUs, required input params: -1 reads1.fq/fa -2 reads2.fa/fq -ref 
# test pipeline on mice reads mapped to MT - extract from dedup bam



###########################################################################################################################################################

# configfile: "config.yaml"
SAMPLES=["sample"]
EXT=["1", "2"]
REF="MN908947.3"
GFF="Sars_cov_2.ASM985889v3.101.gff3"

rule all:
    input:
        expand("results/stats/{sample}.R{ext}.html", sample=SAMPLES, ext=EXT),
        expand("results/stats/{sample}.R{ext}_fastqc.zip", sample=SAMPLES, ext=EXT),
        expand("results/stats/{sample}_multiqc.html", sample=SAMPLES),
        # expand("results/trimmed_reads/{sample}.1_val_1.fq.gz", sample=SAMPLES),
        # expand("results/trimmed_reads/{sample}.1.fastq.gz_trimming_report.txt", sample=SAMPLES),
        # expand("results/trimmed_reads/{sample}.2_val_2.fq.gz", sample=SAMPLES),
        # expand("results/trimmed_reads/{sample}.2.fastq.gz_trimming_report.txt", sample=SAMPLES),
        # expand("reference/{ref}{ext}", ref=REF, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
        expand("reference/{ref}.fasta.fai", ref=REF),
        # expand("results/mapped/{sample}_srt.bam", sample=SAMPLES),
        # expand("results/dedup/{sample}_dedup.bam", sample=SAMPLES),
        # expand("results/dedup/{sample}_dedup.bam.bai", sample=SAMPLES),
        expand("results/dedup/{sample}.dedup.metrics.txt", sample=SAMPLES),
        # expand("results/var_calls/{sample}.vcf", sample=SAMPLES),
        expand("results/stats/{sample}.var.stats", sample=SAMPLES),
        # expand("results/var_filtered/{sample}_QUAL_fltr.vcf", sample=SAMPLES),
        # expand("results/var_filtered/{sample}_DP_fltr.vcf", sample=SAMPLES),
        # "reference/vep/cache",
        "reference/vep/plugins",
        expand("results/var_vep_annotated/{sample}.vep.vcf", sample=SAMPLES),
        expand("results/var_vep_annotated/{sample}.vep.html", sample=SAMPLES)



###########################################################################################################################################################

## INPUT PREPARATION ##


# quality check reads
rule fastqc:
    input:
        "input/{sample}.R{ext}.paired.fq.gz"
    output:
        html="results/stats/{sample}.R{ext}.html",
        zip="results/stats/{sample}.R{ext}_fastqc.zip" 
    log:
        "logs/fastqc/{sample}.R{ext}_fastqc.log"
    wrapper:
        "v1.0.0/bio/fastqc"



# aggregate fastqc reports
rule multiqc:
    input:
        lambda wildcards: expand("results/stats/{sample}.R{ext}_fastqc.zip", sample=wildcards.sample, ext=EXT)
    output:
        "results/stats/{sample}_multiqc.html"
    log:
        "logs/multiqc/{sample}_multiqc.log"
    wrapper:
        "v1.1.0/bio/multiqc"



# QC and adapter trimmming 
rule trim_PE_reads:
    input:
        reads=["input/{sample}.R1.paired.fq.gz", "input/{sample}.R2.paired.fq.gz"]
    output:
        "results/trimmed_reads/{sample}.1_val_1.fq.gz",
        "results/trimmed_reads/{sample}.1.fastq.gz_trimming_report.txt",
        "results/trimmed_reads/{sample}.2_val_2.fq.gz",
        "results/trimmed_reads/{sample}.2.fastq.gz_trimming_report.txt"
    params:
        extra="--gzip -q 20"
    log:
        "logs/trim_galore/{sample}.R1.log",
        "logs/trim_galore/{sample}.R2.log"
    wrapper:
        "v1.0.0/bio/trim_galore/pe"




###########################################################################################################################################################

## READ MAPPING ##


# index reference genome
rule bwa_index:
    input:
        "reference/{ref}.fasta"
    output:
        idx=multiext("reference/{ref}", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index/{ref}.log"
    params:
        algorithm="bwtsw"
    wrapper:
        "v1.0.0/bio/bwa/index"



rule samtools_index:
    input:
        "reference/{ref}.fasta"
    output:
        "reference/{ref}.fasta.fai"
    log:
        "logs/samtools_faidx/{ref}.log"
    wrapper:
        "v1.1.0/bio/samtools/faidx"



# map PE reads to ref. genome and sort
rule map_reads:
    input:
        reads=["results/trimmed_reads/{sample}.1_val_1.fq.gz", "results/trimmed_reads/{sample}.2_val_2.fq.gz"],
        idx=expand("reference/{ref}{ext}", ref=REF, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"])
    output:
        "results/mapped/{sample}_srt.bam"
    log:
        "logs/bwa_mem/{sample}_srt.log"
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'
        sort_order="coordinate"
    threads: 2
    wrapper:
        "v1.1.0/bio/bwa/mem"



# remove duplicates
rule dedup:
    input:
        "results/mapped/{sample}_srt.bam"
    output:
        bam="results/dedup/{sample}_dedup.bam",
        bai="results/dedup/{sample}_dedup.bam.bai",
        metrics="results/dedup/{sample}.dedup.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        extra="--REMOVE_DUPLICATES true --CREATE_INDEX true"
    resources:
        mem_mb=1024
    wrapper:
        "v1.1.0/bio/picard/markduplicates"



###########################################################################################################################################################

## VARIANT CALLING and FILTERING ##


# short variant calling 
rule var_call:
    input:
        ref=expand("reference/{ref}.fasta", ref=REF),
        samples="results/dedup/{sample}_dedup.bam",
        indexes="results/dedup/{sample}_dedup.bam.bai"
    output:
        "results/var_calls/{sample}.vcf"
    log:
        "logs/freebayes/{sample}.log"
    threads: 2
    wrapper:
        "v1.1.0/bio/freebayes"



# check min, max DP before filtering
rule var_stats:
    input:
        "results/var_calls/{sample}.vcf"
    output:
        "results/stats/{sample}.var.stats"
    log:
        "logs/bcftools/{sample}.var.stats.log"
    wrapper:
        "v1.1.0/bio/bcftools/stats"



# quality filter
rule QUAL_filter:
    input:
        "results/var_calls/{sample}.vcf"
    output:
        "results/var_filtered/{sample}_QUAL_fltr.vcf"
    log:
        "logs/bcftools/{sample}_QUAL_fltr.log"
    params:
        filter="-i 'QUAL >= 20'"
    wrapper:
        "v1.1.0/bio/bcftools/filter"



# depth filter
rule DP_filter:
    input:
        "results/var_filtered/{sample}_QUAL_fltr.vcf"
    output:
        "results/var_filtered/{sample}_DP_fltr.vcf"
    log:
        "logs/bcftools/{sample}_DP_fltr.log"
    params:
        filter="-i 'FMT/DP >= 10'"
    wrapper:
        "v1.1.0/bio/bcftools/filter"



###########################################################################################################################################################

## VARIANT ANNOTATION ##


# get annotation files
# rule vep_cache:
#     output:
#         directory("reference/vep/cache")
#     params:
#         species="sars_cov_2",
#         build="ASM985889v3",
#         release="101"
#     log:
#         "logs/vep/cache.log"
#     cache: True  
#     wrapper:
#         "0.70.0/bio/vep/cache"


# get vep plugins
rule vep_plugins:
    output:
        directory("reference/vep/plugins")
    params:
        release=101
    wrapper:
        "v1.1.0/bio/vep/plugins"


# annotate variants
rule vep_annotate_vars:
    input:
        calls="results/var_filtered/{sample}_DP_fltr.vcf", 
        cache="reference/vep/cache",  
        plugins="reference/vep/plugins"
    output:
        calls="results/var_vep_annotated/{sample}.vep.vcf", 
        stats="results/var_vep_annotated/{sample}.vep.html"
    params:
        plugins=["GO"],  # can't omit plugin, mandatory in wrapper script
        # extra="--everything"  
    log:
        "logs/vep/{sample}.var.vep.log"
    threads: 4
    wrapper:
        "0.71.1/bio/vep/annotate"