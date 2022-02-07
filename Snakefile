### VARIANT CALLING PIPELINE ###

# author: Lorena Derezanin
# date:  2/4/2022

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
EXT=[1, 2]
REF="MN908947.3"

rule all:
    input:
        # expand("results/stats/{sample}.R{ext}.html", sample=SAMPLES, ext=EXT),
        # expand("results/stats/{sample}.R{ext}_fastqc.zip", sample=SAMPLES, ext=EXT),
        # expand("results/trimmed_reads/{sample}.{ext}_val_{ext}.fq.gz", sample=SAMPLES, ext=EXT),
        # expand("results/trimmed_reads/{sample}.{ext}.fastq.gz_trimming_report.txt", sample=SAMPLES, ext=EXT),
        # expand("reference/{ref}{ext}", ref=REF, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
        expand("reference/{ref}.fasta.fai", ref=REF),
        # expand("results/mapped/{sample}_srt.bam", sample=SAMPLES),
        # expand("results/dedup/{sample}_dedup.bam", sample=SAMPLES),
        # expand("results/dedup/{sample}_dedup.bam.bai", sample=SAMPLES),
        # expand("results/dedup/{sample}.dedup.metrics.txt", sample=SAMPLES),
        expand("results/var_calls/{sample}.vcf", sample=SAMPLES),

      



###########################################################################################################################################################

## PREPARE INPUT SEQUENCES ##


# quality check reads
rule reads_QC:
    input:
        # ["input/{sample}.R1.paired.fq.gz", "input/{sample}.R2.paired.fq.gz"]
        expand("input/{sample}.R{ext}.paired.fq.gz", sample=SAMPLES, ext=EXT)
    output:
        html="results/stats/{sample}.R{ext}.html",
        zip="results/stats/{sample}.R{ext}_fastqc.zip" 
    log:
        "logs/fastqc/{sample}.R{ext}_fastqc.log"
    wrapper:
        "v1.0.0/bio/fastqc"



# QC and adapter trimmming 
rule trim_PE_reads:
    input:
        expand("input/{sample}.R{ext}.paired.fq.gz", sample=SAMPLES, ext=EXT)
    output:
        "results/trimmed_reads/{sample}.{ext}_val_{ext}.fq.gz",
        "results/trimmed_reads/{sample}.{ext}.fastq.gz_trimming_report.txt"
    params:
        extra="--gzip -q 20"
    log:
        "logs/trim_galore/{sample}.R{ext}.log"
    wrapper:
        "v1.0.0/bio/trim_galore/pe"




###########################################################################################################################################################

## MAP READS ##


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
        reads=expand("input/{sample}.R{ext}.paired.fq.gz", sample=SAMPLES, ext=EXT),
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


# mosdepth - coverage check



###########################################################################################################################################################

## VARIANT CALLING ##


# short variant calling 
rule var_call:
    input:
        ref=expand("reference/{ref}.fasta", ref=REF),
        samples="results/dedup/{sample}_dedup.bam",
        indexes="results/dedup/{sample}_dedup.bam.bai"
        #regions="path/to/region-file.bed"
    output:
        "results/var_calls/{sample}.vcf"
    log:
        "logs/freebayes/{sample}.log"
    threads: 2
    wrapper:
        "v1.1.0/bio/freebayes"



# quality filter
rule qual_filter:
    input:
        "results/var_calls/{sample}.vcf"
    output:
        "results/var_calls/{sample}.QUAL_fltr.vcf"
    log:
        "log/bcftools/{sample}.QUAL_fltr.log"
    params:
        filter="-i 'QUAL >= 20'"
    wrapper:
        "v1.1.0/bio/bcftools/filter"









