
# Variant calling pipeline 


### Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### Recreate conda environment:
`conda env create -f envs/snek.yml`

### Activate environment:
`conda activate snek`

### Install snakemake with mamba: 
`mamba install snakemake`
Robust package manager, handles snakemake releases and dependencies better than conda.

### Run pipeline as:
`snakemake --use-conda --cores 4 --verbose -S Snakefile`

* plot summary - Rmd > knit to html
** het/hom SNPs/INDELs
** genome synteny - gggenomes (R pkg)
** chromomap - heatmap with variants

Run unit tests

Test pipeline on mice reads mapped to MT - extract from dedup bam:
* 1 mouse from control line (high cov set)
* 1 mouse from "marathon" line (high cov set)
* check cov with mosdepth
* downsample if necessary (seqtk)
* bam to fastq read pairs
* map to ref. mouse mitogenome
* run whole pipeline
