process SUBSET_GENES {
    tag "$meta_merged"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    path gpa_file
    path meta_merged
    // path st_file
    // path strain_file
    val min_st_count
    val perc

    output:
    path 'gene_set_ST.tsv', emit: gene_subset

    script:
    """
    subset_genes.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --min_ST_count=$min_st_count \
        --perc=$perc \
        --rm_split=True --ref_only=False \
        --outf=gene_set_ST.tsv
    """
}


process LENGTH_SCALE_COUNTS {
    tag "$merged_lens"
    label 'process_medium'
    publishDir "${params.outdir}/gene_counts", mode: 'copy'

    input:
    tuple path(merged_counts), path(merged_lens)
    path gene_subset

    output:
    path 'kallisto_scaled_counts.tsv', emit: scaled_counts

    script:
    """
    length_scale_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -i TRUE -p TRUE \
        -o kallisto_scaled_counts.tsv
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
    path 'kallisto_tmm_counts.tsv', emit: tmm_counts

    script:
    """
    TMM_normalise_counts.R \
        -c $merged_counts \
        -l $merged_lens \
        -g $gene_subset \
        -r FALSE -p TRUE -t TRUE \
        -o kallisto_tmm_counts.tsv
    """
}
