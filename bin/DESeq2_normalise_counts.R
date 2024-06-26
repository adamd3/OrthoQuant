#!/usr/bin/env Rscript

library(optparse)
library(DESeq2)
if (!require("tidyverse")) {
    lapply(c("readr", "tibble", "purrr", "tidyr", "stringr", "ggplot2"),
        library,
        character.only = TRUE
    )
} else {
    library(tidyverse)
}


option_list <- list(
    make_option(c("-c", "--counts"),
        type = "character", default = NULL,
        help = "table of read counts per gene", metavar = "character"
    ),
    make_option(c("-l", "--lengths"),
        type = "character", default = NULL,
        help = "table of effective lengths per gene", metavar = "character"
    ),
    make_option(c("-g", "--genes"),
        type = "character", default = NULL,
        help = "core gene subset to be used", metavar = "character"
    ),
    make_option(c("-p", "--perc"),
        type = "character", default = NULL,
        help = "filter to core genes?", metavar = "character"
    ),
    make_option(c("-t", "--log_transform"),
        type = "character", default = NULL,
        help = "log transform the counts? default = FALSE", metavar = "character"
    ),
    make_option(c("-o", "--outdir"),
        type = "character", default = NULL,
        help = "output directory for results", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

counts_f <- opt$counts
lengths_f <- opt$lengths
gene_f <- opt$genes
perc <- if (opt$perc == "TRUE") TRUE else FALSE
log <- if (opt$log_transform == "TRUE") TRUE else FALSE
outdir <- opt$outdir


## Read data
counts_tab <- suppressMessages(read_tsv(counts_f))
lengths_tab <- suppressMessages(read_tsv(lengths_f))
core_genome <- if (!is.null(gene_f)) {
    suppressMessages(read_tsv(gene_f))
} else {
    NULL
}

lengths_tab <- lengths_tab[match(counts_tab$Gene, lengths_tab$Gene), ]

gene_ids <- counts_tab$Gene
counts_tab <- data.frame(counts_tab[, 2:ncol(counts_tab)])
lengths_tab <- data.frame(lengths_tab[, 2:ncol(lengths_tab)])
rownames(counts_tab) <- rownames(lengths_tab) <- gene_ids

## scale counts to reads per median gene length
median_lens <- rowMedians(as.matrix(lengths_tab), na.rm = TRUE)
names(median_lens) <- gene_ids

counts_tab_scaled <- (counts_tab / lengths_tab)
counts_tab_scaled <- sweep(counts_tab_scaled, 1, median_lens, "*")


colData <- data.frame(sample_name = colnames(counts_tab_scaled))

# get scaling factors from core gene set
counts_tab_core <- na.omit(counts_tab_scaled)

dds <- suppressMessages(DESeqDataSetFromMatrix(
    countData = round(counts_tab_core),
    colData = colData, design = ~1
))

dds <- estimateSizeFactors(dds)
size_factors <- sizeFactors(dds)



if (isTRUE(perc)) {
    ## subset to `perc` gene set
    counts_tab_perc <- subset(counts_tab_scaled, rownames(
        counts_tab_scaled
    ) %in% core_genome$gene)
    counts_tab_perc[is.na(counts_tab_perc)] <- 0

    median_sub <- median_lens[rownames(counts_tab_perc)]

    ## get size factor-scaled counts
    res_df <- sweep(counts_tab_perc, 2, size_factors, `/`)

    ## get RPKM values
    lib_sf <- (colSums(counts_tab_perc, na.rm = TRUE) / 1e6)
    gene_sf <- sapply(median_sub, function(l) {
        l / 1e3 * lib_sf
    })
    gene_sf <- as.data.frame(t(gene_sf))
    rpkm_df <- counts_tab_perc / gene_sf
} else {
    ## get size factor-scaled counts
    res_df <- sweep(counts_tab_scaled, 2, size_factors, `/`)

    ## get RPKM values
    lib_sf <- (colSums(counts_tab_scaled, na.rm = TRUE) / 1e6)
    gene_sf <- sapply(median_lens, function(l) {
        l / 1e3 * lib_sf
    })
    gene_sf <- as.data.frame(t(gene_sf))
    rpkm_df <- counts_tab_scaled / gene_sf
}


if (isTRUE(log)) {
    res_df <- log2(res_df + 1)
    rpkm_df <- log2(rpkm_df + 1)
}

## convert rownames to column
res_df <- tibble::rownames_to_column(as.data.frame(res_df), "feature_id")
rpkm_df <- tibble::rownames_to_column(as.data.frame(rpkm_df), "feature_id")
counts_tab_scaled <- tibble::rownames_to_column(as.data.frame(
    counts_tab_scaled
), "feature_id")

write.table(
    counts_tab_scaled, file.path(outdir, "raw_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

write.table(
    res_df, file.path(outdir, "norm_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)

write.table(
    rpkm_df, file.path(outdir, "rpkm_counts.tsv"),
    col.names = TRUE, row.names = FALSE,
    sep = "\t", quote = FALSE
)
