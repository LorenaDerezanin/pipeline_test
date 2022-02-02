# Pipeline_test

Variant calling pipeline and statistical summary of metrics

* trim_galore
* FastQC or MultiQC
* bwa mem
* samtools sam > bam > sorted bam
* Picard MarkDups - rmdups
* variant calling - FreeBayes
* QC filtering
* plot summary - Rmd > knit to html
** het/hom SNPs/INDELs
** SNP/INDEL per 1kb-window
** genome synteny - gggenomes (R pkg)
** chromomap - heatmap with variants

Run in conda env, include yml file to recreate the env
Run unit tests
Run on 2 cores by default (to avoid using up all local CPUs)

Test pipeline on mice reads mapped to MT - extract from dedup bam:
* 1 mouse from control line (high cov set)
* 1 mouse from "marathon" line (high cov set)
* check cov with mosdepth
* downsample if necessary (seqtk)
* bam to fastq read pairs
* map to ref. mouse mitogenome
* run whole pipeline
