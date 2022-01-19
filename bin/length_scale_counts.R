#!/usr/bin/env Rscript

## load / install packages
if (!require("optparse")){
    install.packages("optparse")
}
if (!require("matrixStats")){
    install.packages("matrixStats")
}


################################################
## Arguments
option_list <- list(
    make_option(c("-c", "--counts"), type="character", default=NULL,
        help="table of read counts per gene", metavar="character"),
    make_option(c("-l", "--lengths"), type="character", default=NULL,
        help="table of effective lengths per gene", metavar="character"),
    make_option(c("-g", "--genes"), type="character", default=NULL,
        help="core gene subset to be used", metavar="character"),
    make_option(c("-i", "--impute"), type="character", default=NULL,
        help="return 0-imputed counts", metavar="character"),
    make_option(c("-p", "--perc"), type="character", default=NULL,
        help="was filtering based on percentage presence of gene?", metavar="character"),
    make_option(c("-o", "--outf"), type="character", default=NULL,
        help="output file for results", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

counts_f <- opt$counts
lengths_f <- opt$lengths
gene_f <- opt$genes
impute <- if(opt$impute == "TRUE") TRUE else FALS
perc <- if(opt$perc == "TRUE") TRUE else FALSE
outf <- opt$outf
################################################


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


counts_tab <- subset(counts_tab, Gene %in% core_genes)
lengths_tab <- subset(lengths_tab, Gene %in% core_genes)

lengths_tab <- lengths_tab[match(rownames(counts_tab), rownames(lengths_tab)),]


rownames(counts_tab) <- rownames(lengths_tab) <- counts_tab$Gene
counts_tab$Gene <- lengths_tab$Gene <- NULL

## scale counts to reads per median gene length
median_lens <- rowMedians(as.matrix(lengths_tab), na.rm=TRUE)
counts_tab_scaled <- (counts_tab/lengths_tab)
counts_tab_scaled <- sweep(counts_tab_scaled, 1, median_lens, "*")

if(isTRUE(impute)){
    ## Replace missing values with 0
    counts_tab_scaled[is.na(counts_tab_scaled)] <- 0
}

write.table(
    counts_tab_scaled, outf, col.names = TRUE, row.names = TRUE,
    sep = "\t", quote = FALSE
)