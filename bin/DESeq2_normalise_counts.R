#!/usr/bin/env Rscript


if (!require("optparse")){
    install.packages("optparse")
}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("DESeq2")){
    BiocManager::install("DESeq2")
}

library(optparse)
library(DESeq2)


option_list <- list(
    make_option(c("-c", "--counts"), type="character", default=NULL,
        help="table of read counts per gene", metavar="character"),
    make_option(c("-l", "--lengths"), type="character", default=NULL,
        help="table of effective lengths per gene", metavar="character"),
    make_option(c("-g", "--genes"), type="character", default=NULL,
        help="core gene subset to be used", metavar="character"),
    make_option(c("-p", "--perc"), type="character", default=NULL,
        help="was filtering based on percentage presence of gene?", metavar="character"),
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
counts_tab <- read.csv(
    counts_f, header = TRUE, na.strings=c("","NA"), sep = "\t",
    stringsAsFactors = FALSE
)
lengths_tab <- read.csv(
    lengths_f, header = TRUE, na.strings=c("","NA"), sep = "\t",
    stringsAsFactors = FALSE
)
core_genome <- read.csv(
    gene_f, header = TRUE, na.strings=c("","NA"), sep = "\t",
    stringsAsFactors = FALSE
)


if(isFALSE(perc)){
    colnames(core_genome) <- gsub("X", "ST_", colnames(core_genome))
    na_count <- sapply(core_genome, function(y) sum(length(which(is.na(y)))))
    n <- 1
    least_na <- names(sort(na_count)[1:n])
    core_genes <- na.omit(core_genome[least_na])[,1]
} else {
    core_genes <- core_genome$gene
}

colnames(counts_tab)[1] <- colnames(lengths_tab)[1] <- "Gene"

counts_tab <- subset(counts_tab, Gene %in% core_genes)
lengths_tab <- subset(lengths_tab, Gene %in% core_genes)


rownames(counts_tab) <- rownames(lengths_tab) <- counts_tab$Gene
counts_tab$Gene <- lengths_tab$Gene <- NULL

lengths_tab <- lengths_tab[match(rownames(counts_tab), rownames(lengths_tab)),]


## scale counts to reads per median gene length
median_lens <- rowMedians(as.matrix(lengths_tab), na.rm=TRUE)
counts_tab_scaled <- (counts_tab/lengths_tab)
counts_tab_scaled <- sweep(counts_tab_scaled, 1, median_lens, "*")

## Replace missing values with 0
counts_tab_scaled[is.na(counts_tab_scaled)] <- 0


colData <- data.frame(sample_name = colnames(counts_tab_scaled))

dds <- suppressMessages(DESeqDataSetFromMatrix(
    countData = round(counts_tab_scaled),
    colData = colData,
    design = ~ 1
))

dds <- estimateSizeFactors(dds)
res_df <- counts(dds, normalized=TRUE)


## get RPKM per gene
mcols(dds)$basepairs <- median_lens
rpkm_df <- fpkm(dds, robust = TRUE)


if(isTRUE(log)){
    res_df <- log2(res_df+1)
    # rpkm_df <- log2(rpkm_df+1) ## update: don't log transform the RPKM vals
}

write.table(
    res_df, file.path(outdir,"norm_counts.tsv"),
    col.names = TRUE, row.names = TRUE,
    sep = "\t", quote = FALSE
)

write.table(
    rpkm_df, file.path(outdir,"rpkm_counts.tsv"),
    col.names = TRUE, row.names = TRUE,
    sep = "\t", quote = FALSE
)
