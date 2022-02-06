### VARIANT CALLING PIPELINE ###

# author: Lorena Derezanin
# date:  2/4/2022

# snakemake v
# run in conda env


# add unit tests before rule chunks
# add config with default params - 2 CPUs, required input params: -1 reads1.fq/fa -2 reads2.fa/fq -ref 
# test pipeline on mice reads mapped to MT - extract from dedup bam



###########################################################################################################################################################


rule all:
    input:





###########################################################################################################################################################

## PREPARE SEQUENCES ##

rule trim_reads:

trim_galore




rule reads_QC:

fastQC/MultiQC




rule map_reads:

bwa mem


rule sam_to_bam:

samtools view



rule sort_bam:

samtools sort



rule rm_dups:

picard MarkDuplicates


###########################################################################################################################################################

## VARIANT CALLING ##


rule var_call:

FreeBayes


rule QC_filter:







