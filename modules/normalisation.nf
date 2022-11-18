process SUBSET_GENES {
    tag "$meta_merged"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    path gpa_file
    path meta_merged
    val perc

    output:
    path 'gene_set_ST.tsv', emit: gene_subset

    script:
    """
    subset_genes.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --perc=$perc \
        --ref_only=False \
        --outf=gene_set_ST.tsv
    """
}


// process LENGTH_SCALE_COUNTS {
//     tag "$merged_lens"
//     label 'process_medium'
//     publishDir "${params.outdir}/gene_counts", mode: 'copy'

//     input:
//     tuple path(merged_counts), path(merged_lens)
//     path gene_subset

//     output:
//     path 'raw_counts.tsv', emit: scaled_counts

//     script:
//     """
//     length_scale_counts.R \
//         -c $merged_counts \
//         -l $merged_lens \
//         -g $gene_subset \
//         -i TRUE -p TRUE \
//         -o raw_counts.tsv
//     """
// }

process DESEQ_NORMALISE_COUNTS {
    tag "$merged_counts"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    tuple path(merged_counts), path(merged_lens)
    path gene_subset

    output:
    path 'norm_counts.tsv', emit: norm_counts
    path 'rpkm_counts.tsv', emit: rpkm_counts
    path 'raw_counts.tsv', emit: scaled_counts

    script:
    """
    DESeq2_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -p TRUE -t TRUE \
        -o ./
    """
}

process TMM_NORMALISE_COUNTS {
    tag "$merged_counts"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    tuple path(merged_counts), path(merged_lens)
    path gene_subset

    output:
    path 'norm_counts.tsv', emit: norm_counts
    path 'rpkm_counts.tsv', emit: rpkm_counts
    path 'raw_counts.tsv', emit: scaled_counts

    script:
    """
    TMM_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -p TRUE -t TRUE \
        -o ./
    """
}
