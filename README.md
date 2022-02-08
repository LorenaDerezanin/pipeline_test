
# Variant calling pipeline 

A Snakemake workflow for calling and annotation of short variants.<br/> 
Workflow takes paired-end Illumina short read data (fastq files) as input and outputs annotated variant calls in a vcf file as the final result.
Input directory contains PE Illumina reads from a publicly available SARS-CoV-2 dataset [SRA accession SRR15660643](https://www.ncbi.nlm.nih.gov/sra/?term=SRR15660643) downsampled to 16000 paired reads (sample.R1.paired.fq.gz and sample.R2.paired.fq.gz).<br/> 
A fasta file with the Wuhan-Hu-1 reference genome [Genbank accession MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3) is included in the<br/>
reference directory (MN908947.3.fasta), along with the VEP cache for successful annotation of genomic features.

## Usage

`git clone git@github.com:LorenaDerezanin/pipeline_test.git` 

### Step 1: Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

Minimal conda installer for running pipeline in an isolated conda environment to avoid dependency hell and ensure reproducibility.

### Step 2 (Recommended): Install mamba - faster package manager 

`conda install mamba -n base -c conda-forge`

Recommended installation to speed up env setup. Mamba is a more robust and faster package manager (parallel download of data), and handles releases and dependencies better than conda. If continuing with `conda`, `mamba` should be replaced with `conda` in Step 3.

### Step 3: Recreate conda environment

`cd pipeline_test/`
`mamba env create -n snek -f envs/snek.yml`

### Step 4: Activate environment

`conda activate snek`


### Step 5: Run pipeline

`snakemake --use-conda --cores 4 --verbose`


## Troubleshooting

If conda fails to install `snakemake v.6.15`, install snakemake with mamba: `mamba install snakemake`.


## Pipeline content

Bioinformatic tools used in the Snakemake workflow in the form of Snakemake Wrappers obtained from [The Snakemake Wrappers Repository](https://snakemake-wrappers.readthedocs.io/en/stable/):

* fastQC
* multiQC
* trim_galore
* bwa 
* samtools
* picard
* freebayes
* bcftools
* vep


