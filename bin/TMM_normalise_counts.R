#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
library(edgeR)
library(tidyverse)

option_list <- list(
    make_option(c("-c", "--counts"), type="character", default=NULL,
        help="table of read counts per gene", metavar="character"),
    make_option(c("-l", "--lengths"), type="character", default=NULL,
        help="table of effective lengths per gene", metavar="character"),
    make_option(c("-g", "--genes"), type="character", default=NULL,
        help="core gene subset to be used", metavar="character"),
    make_option(c("-p", "--perc"), type="character", default=NULL,
        help="filter to core genes?", metavar="character"),
    make_option(c("-t", "--log_transform"), type="character", default=NULL,
        help="log transform the counts? default = FALSE", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default=NULL,
        help="output directory for results", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

counts_f <- opt$counts
lengths_f <- opt$lengths
gene_f <- opt$genes
perc <- if(opt$perc == "TRUE") TRUE else FALSE
log <- if(opt$log_transform == "TRUE") TRUE else FALSE
outdir <- opt$outdir


## Read data
counts_tab <- suppressMessages(read_tsv(counts_f))
lengths_tab <- suppressMessages(read_tsv(lengths_f))
core_genome <- suppressMessages(read_tsv(gene_f))

lengths_tab <- lengths_tab[match(counts_tab$Gene, lengths_tab$Gene),]

gene_ids <- counts_tab$Gene
counts_tab <- data.frame(counts_tab[,2:ncol(counts_tab)])
lengths_tab <- data.frame(lengths_tab[,2:ncol(lengths_tab)])
rownames(counts_tab) <- rownames(lengths_tab) <- gene_ids

## scale counts to reads per median gene length
median_lens <- rowMedians(as.matrix(lengths_tab), na.rm=TRUE)
names(median_lens) <- gene_ids

counts_tab_scaled <- (counts_tab/lengths_tab)
counts_tab_scaled <- sweep(counts_tab_scaled, 1, median_lens, "*")


colData <- data.frame(sample_name = colnames(counts_tab_scaled))

## subset to core genes
counts_tab_sub <- subset(counts_tab_scaled, rownames(
    counts_tab_scaled) %in% core_genome$gene)
median_sub <- median_lens[rownames(counts_tab_sub)]

## replace missing values with 0
counts_tab_sub[is.na(counts_tab_sub)] <- 0



## get size factors per sample, based on core genes
y <- DGEList(
    counts = counts_tab_sub,
    genes = data.frame(gene.length = median_sub)
)
y <- calcNormFactors(y, method = "TMM")

size_factors <- effectiveLibSizes(y)

## NOTE: = effectiveLibSizes(y) = (y$samples$lib.size)*(y$samples$norm.factors)

## NOTE: the output of edgeR effectiveLibSizes() is equivalent to the output of 
## sizeFactors() in DESeq2; see: https://support.bioconductor.org/p/46779/

if(isTRUE(perc)){
    ## get size factor-scaled counts
    res_df <- as.data.frame(cpm(y, log = log))
    ## get RPKM
    rpkm_df <- as.data.frame(edgeR::rpkm(y, log = log)) 

} else {
    ## get size factor-scaled counts
    res_df <- sweep(counts_tab_scaled, 2, size_factors, `/`)
    ## get RPKM values
    sf <- colSums(counts_tab_scaled, na.rm=TRUE)/1e6
    rpkm_df <- sweep(counts_tab_scaled, 2, sf, `/`)

    if(isTRUE(log)){
        res_df <- log2(res_df+1)
        rpkm_df <- log2(rpkm_df+1) 
    }

}

## convert rownames to column
res_df <- tibble::rownames_to_column(as.data.frame(res_df), "feature_id")
rpkm_df <- tibble::rownames_to_column(as.data.frame(rpkm_df), "feature_id")
counts_tab_scaled <- tibble::rownames_to_column(as.data.frame(
    counts_tab_scaled), "feature_id")

write.table(
    counts_tab_scaled, file.path(outdir,"raw_counts.tsv"), 
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

write.table(
    res_df, file.path(outdir,"norm_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

write.table(
    rpkm_df, file.path(outdir,"rpkm_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)
