process MERGE_METADATA {
    tag "$metadata"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata
    path sample_ID

    output:
    path 'metadata_merged.tsv', emit: meta_merged
    path 'id_mappings.csv', emit: id_mappings

    script:
    """
    merge_metadata.py $metadata $id_mappings ${params.data_dir} metadata_merged.tsv
    """
}
