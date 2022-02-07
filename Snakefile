### VARIANT CALLING PIPELINE ###

# author: Lorena Derezanin
# date:  2/4/2022

# snakemake v.6.15.2.
# conda v.4.11.0
# mamba v.0.15.3 - more robust pkg manager than conda, handles snakemake releases better
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
        expand("reference/{ref}{ext}", ref=REF, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
        expand("results/mapped/{sample}_srt.bam", sample=SAMPLES),

      



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
        "results/trimmed_reads/{sample}.{ext}.fastq.gz_trimming_report.txt",
    params:
        extra="--gzip -q 20"
    log:
        "logs/trim_galore/{sample}.R{ext}.log"
    wrapper:
        "v1.0.0/bio/trim_galore/pe"




###########################################################################################################################################################

## MAP READS ##


# index reference genome
rule ref_index:
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


# map PE reads to ref. genome
rule map_reads:
    input:
        reads=expand("input/{sample}.R{ext}.paired.fq.gz", sample=SAMPLES, ext=EXT),
        # idx="reference/{ref}"
        idx=expand("reference/{ref}{ext}", ref=REF, ext=[".amb", ".ann", ".bwt", ".pac", ".sa"]),
    output:
        "results/mapped/{sample}_srt.bam",
    log:
        "logs/bwa_mem/{sample}_srt.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'
        sort_order="coordinate", 
        # sort_extra="",  
    threads: 2
    wrapper:
        "v1.1.0/bio/bwa/mem"



# remove duplicates
rule rm_dups:





###########################################################################################################################################################

## VARIANT CALLING ##


# rule var_call:

# FreeBayes


# rule QC_filter:








