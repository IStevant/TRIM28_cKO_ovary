# Understanding TRIM28 molecular function involved in granulosa-to-Sertoli cell transdifferentiation

## Introduction

The following pipeline was created to analyse CUTandRUN and ATAC-seq data produced frm purified somatic cells from 8-week-old ovaries, in control and Trim28 cKO mice.

The analysis consists first in analysing the H3K9me3 histone mark in control and cKO ovaries to evaluate the effect of Trim28 deletion in heterochromatin.

Second, we profiled the changes in chromatin accessibility to identify whether TRIM28 acts on active cis-regulatory elements. We identified the transcription factor binding sites enriched in the regions gaining and losing accessibility.

Finally, we combined ChIP-seq data of TRIM28 and four ovarian transcription factors which motifs were found enriched in regions affected by TRIM28 loss and identified that TRIM28 might contribute to the stabilisation of ovarian transcription factor hubs.

## How to run the pipeline

### Raw data mapping
Prior to run the pipeline, the raw CUTandRUn (v3.2.2) and ATAC-seq data should be processed using the nf-core/cutandrun and nf-core/atacseq (v2.1.2) pipelines.

The raw data are available on GEO (GSE339481 and GSE339478).

### Install Conda

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
# Specify your installation directory
```

### Clone the Trim28 git project

The following command downloads the last version of the pipeline in the Trim28_cKO_ovary directory.

```bash
git clone git@github.com:IStevant/Trim28_cKO_ovary.git
```

### Install the Conda trim28 environment

This installs a conda environment with R, Snakemake, as well as all the dependencies necessary to run the pipeline. It can some time to install.

```bash
cd Trim28_cKO_ovary
conda env create -f conda_env/trim28_environment.yml
```

### Install the R packages

#### Install missing dependencies

```bash
snakemake --core 1 -f install_packages
```

#### Run the full pipeline

Note that the pipeline was developped to work on an HPC with slurm. It might be addapted to run in another type of system.

```bash
snakemake --profile=slurm/ 
```
