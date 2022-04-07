# StrainSeq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
<!-- [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) -->

## Introduction

**StrainSeq** is a Nextflow pipeline for performing bacterial RNA-Seq analysis without a reference genome.

## Pipeline summary

The pipeline will perform the following steps:

1. Trim adaptors from reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Subset to genes of interest (e.g. the core gene set; or by majority sequence type)
4. Pseudo-alignment to strain-specific gene sets ([`kallisto`](https://pachterlab.github.io/kallisto/))
5. Merging and length-scaling of counts for orthologous genes
6. Size-factor scaling of merged counts ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or [`edgeR`](http://bioconductor.org/packages/release/bioc/html/edgeR.html))
7. Visualisation of gene expression across strains ([`UMAP`](https://umap-learn.readthedocs.io/))

## Required input


* Metadata file. Must contain the below named columns (additional columns are optional):
    ```console
    
    ```
