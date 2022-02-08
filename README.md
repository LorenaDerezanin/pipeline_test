
# Variant calling pipeline 

A Snakemake workflow for calling and annotation of short variants.

## Usage

`git clone git@github.com:LorenaDerezanin/pipeline_test.git` 

### Step 1: Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

Minimal conda installer for running pipeline in an isolated conda environment to avoid dependency hell and ensure reproducibility.

### Step 2 (Recommended): Install mamba - faster package manager 

`conda install mamba -n base -c conda-forge`

Recommended installation to speed up env setup. Mamba is a more robust and faster package manager (parallel download of data), and handles releases and dependencies better than conda. If continuing with `conda`, `mamba` should be replaced with `conda` in Step 3.

### Step 3: Recreate conda environment

`mamba env create -n snek -f envs/snek.yml`

### Step 4: Activate environment

`conda activate snek`


### Step 5: Run pipeline

`snakemake --use-conda --cores 4 --verbose`


## Troubleshooting

If conda fails to install `snakemake v.6.15`, install snakemake with mamba: `mamba install snakemake`.


## Pipeline content

Bioinformatic tools used in the Snakemake workflow:

* fastQC
* multiQC
* trim_galore
* bwa 
* samtools
* picard
* freebayes
* bcftools
* vep


