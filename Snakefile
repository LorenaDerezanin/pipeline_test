### VARIANT CALLING PIPELINE ###

# author: Lorena Derezanin
# date:  2/4/2022

# snakemake v.6.15.2.
# conda v.4.11.0
# mamba v.0.15.3 - more robust pkg manager than conda, handles snakemake releases better
# run in mamba env snek2
# snakemake --use-conda


# add unit tests 
# add config with default params - 2 CPUs, required input params: -1 reads1.fq/fa -2 reads2.fa/fq -ref 
# test pipeline on mice reads mapped to MT - extract from dedup bam



###########################################################################################################################################################

configfile: "config.yaml"

rule all:
    input:
        expand("results/stats/{sample}.html", sample=config["samples"]),
        expand("results/stats/{sample}_fastqc.zip", sample=config["samples"])




###########################################################################################################################################################

## PREPARE INPUT SEQUENCES ##


# quality check reads

rule reads_QC:
    input:
        expand("input/{sample}.fq.gz", sample=config["samples"])
    output:
        html="results/stats/{sample}.html",
        zip="results/stats/{sample}_fastqc.zip"
    # params:
    #     ""  
    log:
        "logs/fastqc/{sample}_fastqc.log"
    wrapper:
        "v1.0.0/bio/fastqc"


# QC and adapter trimmming 

rule trim_PE_reads:
    input:
        expand("input/{sample}.fq.gz", sample=config["samples"])
    output:
        "trimmed_reads/{sample}.1_val_1.fq.gz",
        "trimmed_reads/{sample}.1.fastq.gz_trimming_report.txt",
        "trimmed_reads/{sample}.2_val_2.fq.gz",
        "trimmed_reads/{sample}.2.fastq.gz_trimming_report.txt",
    params:
        extra="-q 20",
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v1.0.0/bio/trim_galore/pe"








# rule map_reads:

# bwa mem


# rule sam_to_bam:

# samtools view



# rule sort_bam:

# samtools sort



# rule rm_dups:

# picard MarkDuplicates


###########################################################################################################################################################

## VARIANT CALLING ##


# rule var_call:

# FreeBayes


# rule QC_filter:








