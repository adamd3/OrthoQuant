#!/usr/bin/env Rscript

library(umap)
library(ggplot2)
library(RColorBrewer)

option_list <- list(
    make_option(c("-n", "--norm_counts"), type="character", default=NULL,
        help="read counts normalised for both gene length and library size", metavar="character"),
    make_option(c("-m", "--metadata_merged"), type="character", default=NULL,
        help="merged metadata mapping strain ID to RNA ID", metavar="character"),
    make_option(c("-g", "--group"), type="character", default=NULL,
        help="metadata column for grouping strains", metavar="character"),
    make_option(c("-o", "--outdir"), type="character", default=NULL,
        help="directory for results", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

counts_f <- opt$norm_counts
meta_f <- opt$metadata_merged
group <- opt$group
outdir <- opt$outdir


large_disc_pal <- brewer.pal.info[brewer.pal.info$category == 'qual',]
colpal_large <- unlist(
    mapply(brewer.pal, large_disc_pal$maxcolors, rownames(large_disc_pal)))
colpal_large[c(5:8)] <- colpal_large[c(70:73)] ## replace to avoid colour clashes


##------------------------------------------------------------------------------
## Read and process data
##------------------------------------------------------------------------------
norm_counts <- read.table(counts_f)
norm_counts <- log2(norm_counts+1)

clone_meta <- read.table(
    meta_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE
)
clone_meta[[group]] <- as.factor(clone_meta[[group]])
# clone_meta$infection_type <- as.factor(clone_meta$infection_type)
# clone_meta$level7000 <- as.factor(clone_meta$level7000)
# clone_meta$ST <- as.factor(clone_meta$ST)

clone_meta_sub <- subset(clone_meta, sample_name %in% colnames(norm_counts))
clone_meta_sub <- clone_meta_sub[match(colnames(norm_counts),clone_meta_sub$sample_name),]


##------------------------------------------------------------------------------
## UMAP
##------------------------------------------------------------------------------
# custom.config <- umap.defaults
# custom.config$n_neighbors <- 20

umap_adj <- umap(t(norm_counts))#, config=custom.config)
umap_adj_dims <- data.frame(umap_adj$layout)
colnames(umap_adj_dims) <- paste0("UMAP_", c(1,2))
umap_adj_dims$sample_name <- rownames(umap_adj_dims)


## add metadata
umap_adj_dims <- merge(umap_adj_dims, clone_meta_sub, all.x=TRUE)

## save umap object
saveRDS(umap_adj, file.path(outdir,'umap.rds'))

p1 <- ggplot(umap_adj_dims, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(size = 4,  shape = 21, colour = "black",
        aes_string(fill = group)
    ) +
    scale_colour_manual(values = colpal_large, guide = "none") +
    scale_fill_manual(values = colpal_large) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_bw(base_size=12) +
    theme(
        legend.position = "right",
        axis.text.x = element_text(colour = "black", size=12),
        axis.text.y = element_text(colour = "black", size=12)
    )

ggsave(
    p1, file = file.path(outdir,paste0('umap_',group,'.png')),
    device = "png", units = "in",
    width = 9, height = 7, dpi = 300
)
