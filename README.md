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


- Metadata file. Must contain the following named columns:
  - `RNA_sample_id`: transcriptome sample identifier
  - `DNA_sample_id`: genome sample identifier
  - `sample_name`: strain name
  - `fastq`: path to fastq file for RNA-Seq data
  - `fasta`: path to fasta file containing gene sequences
  Additional columns are optional. See below example:

    ```console
    RNA_sample_id	DNA_sample_id	sample_name	majority_ST	level7000	ST	patient	infection_type	fastq	fasta
    SRX5123744	SRR8737283	ZG302367	274	4	274	7ZG302367	ear infection	/path/to/fastq/SRX5123744_T1.fastq.gz	/path/to/fasta/SRR8737283.fna
    SRX5123743	SRR8737284	ZG302359	244	5	244	7ZG302359	pyrexia	/path/to/fastq/SRX5123743_T1.fastq.gz	/path/to/fasta/SRR8737284.fna
    SRX5123714	SRR8737286	PSAE1649	313	2	313	PSAE1649	wound infection	/path/to/fastq/SRX5123714_T1.fastq.gz	/path/to/fasta/SRR8737286.fna
    SRX5123695	SRR8737287	MHH17441	235	6	235	4MHH17441	urinary tract	/path/to/fastq/SRX5123695_T1.fastq.gz	/path/to/fasta/SRR8737287.fna
    SRX5123726	SRR8737288	PSAE1975	395	7	395	1PSAE1975	wound infection	/path/to/fastq/SRX5123726_T1.fastq.gz	/path/to/fasta/SRR8737288.fna
    SRX5123719	SRR8737290	PSAE1745	111	3	111	3PSAE1745	respiratory tract	/path/to/fastq/SRX5123719_T1.fastq.gz	/path/to/fasta/SRR8737290.fna
    ```
