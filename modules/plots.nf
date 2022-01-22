process UMAP_SAMPLES {
    tag "$tmm_counts"
    label 'process_medium'
    publishDir "${params.outdir}/umap_samples", mode: 'copy'

    input:
    path tmm_counts
    path meta_merged

    output:
    path '*.{rds,png}', emit: umap_out

    script:
    """
    umap.R -n $tmm_counts -m $meta_merged -o ./
    """
}
