# StrainSeq

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
<!-- [![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) -->

## Introduction

**StrainSeq** is a Nextflow pipeline for performing strain-specific bacterial RNA-Seq analysis without a reference genome.

## Pipeline summary

The pipeline requires the output from a pan-genome analysis with ([`Panaroo`](https://gtonkinhill.github.io/panaroo/)) and will perform the following steps:

1. Trim adaptors from reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Subset to genes of interest (e.g. the core gene set; or by majority sequence type)
4. Pseudo-alignment to strain-specific gene sets ([`kallisto`](https://pachterlab.github.io/kallisto/))
5. Merging and length-scaling of counts for orthologous genes
6. Size-factor scaling of merged counts ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) or [`edgeR`](http://bioconductor.org/packages/release/bioc/html/edgeR.html))
7. Visualisation of gene expression across strains ([`UMAP`](https://umap-learn.readthedocs.io/))

## Required input


- __Metadata file__: Tab-delimited file, which must contain at least the following named columns:
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

- __Gene presence-absence file__: CSV-format output produced by ([`Panaroo`](https://gtonkinhill.github.io/panaroo/)).
  See below example (truncated):

    ```console
    Gene,Non-unique Gene name,Annotation,Pseudomonas_aeruginosa_PAO1_107_converted,SRR8737281,SRR8737282
    group_25769,,Protein of unknown function (DUF2845),PGD112333,SRR8737281_00662,SRR8737282_00657
    group_25760,;hypothetical protein,hypothetical protein;topoisomerase IIProtein of unknown function (DUF2790);Protein of unknown function (DUF2790);acyl- n-acyltransferaseUncharacterized protein conserved in bacteriaProtein of unknown function DUF482;,PGD109272,SRR8737281_05540,SRR8737282_06520
    group_25759,;hypothetical protein,Protein of unknown function (DUF3309);,PGD109276,SRR8737281_05538,SRR8737282_03481
    rubA2~~~Rubredoxin2,rubA2;Rubredoxin 2,rubredoxinRubredoxin-2anaerobic nitric oxide reductase flavorubredoxinRubredoxinRubredoxin;,PGD113575,SRR8737281_00021,SRR8737282_05845
    bfr_2~~~bacterioferritin,bfr_2;bacterioferritin,bacterioferritinBacterioferritinbacterioferritinBacterioferritin (cytochrome b1)bacterioferritinFerritin-like domain;,PGD109882,SRR8737281_05955,SRR8737282_01475
    mtnB~~~probablesugaraldolase,mtnB;probable sugar aldolase,methylthioribulose-1-phosphate dehydrataseMethylthioribulose-1-phosphate dehydratasemethylthioribulose-1-phosphate dehydrataseRibulose-5-phosphate 4-epimerase and related epimerases and aldolasesmethylthioribulose-1-phosphate dehydrataseClass II Aldolase and Adducin N-terminal domain;,PGD106137,SRR8737281_02504,SRR8737282_00940
    epd_2~~~epd_1~~~epd,epd_2;epd_1;epd,D-erythrose-4-phosphate dehydrogenaseD-erythrose-4-phosphate dehydrogenaseerythrose 4-phosphate dehydrogenaseerythrose-4-phosphate dehydrogenaseGlyceraldehyde 3-phosphate dehydrogenase C-terminal domain;D-erythrose-4-phosphate dehydrogenaseD-erythrose-4-phosphate dehydrogenaseerythrose 4-phosphate dehydrogenaseTransketolaseerythrose-4-phosphate dehydrogenaseGlyceraldehyde 3-phosphate dehydrogenase C-terminal domain,PGD103840,SRR8737281_03107,SRR8737282_02473
    emrE,emrE,SMR multidrug efflux transporterMethyl viologen resistance protein Cmultidrug efflux proteinMembrane transporters of cations and cationic drugsphosphonate utilization associated putative membrane proteinSmall Multidrug Resistance protein;SMR multidrug efflux transporterMethyl viologen resistance protein Cmultidrug efflux proteinMembrane transporters of cations and cationic drugsSmall Multidrug Resistance protein,PGD112851,SRR8737281_00400,SRR8737282_00337
    purA_1~~~purA_2~~~purA,purA_1;purA_2;purA,adenylosuccinate synthetaseAdenylosuccinate synthetaseadenylosuccinate synthetaseadenylosuccinate synthaseAdenylosuccinate synthetase,PGD112747,SRR8737281_00453,SRR8737282_00390
    ```

## Output

1. __trim_galore__ directory containing adaptor-trimmed RNA-Seq files and FastQC results.
2. __gene_counts__ directory containing:
    1. `gene_set_ST.tsv`: the subset of genes included in analysis.
    2. `raw_counts.tsv`: merged read counts per gene, scaled to the median gene length across strains.
    3. `norm_counts.tsv`: size factor scaled counts (normalised for library size).
    4. `rpkm_counts.tsv`: size factor scaled and gene length-scaled counts, expressed as reads per kilobase per million mapped reads (RPKM) (normalised for library size and gene length).
3. __umap_samples__ directory containing UMAP visualisation of gene expression across strains.
