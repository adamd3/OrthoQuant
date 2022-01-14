process MERGE_METADATA {
    tag "$metadata"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata

    output:
    path 'metadata_merged.tsv', emit: meta_merged

    script:
    """
    merge_metadata.py $metadata ${params.data_dir} metadata_merged.tsv
    """
}
