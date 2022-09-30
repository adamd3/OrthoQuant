process MERGE_METADATA {
    tag "$metadata"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata
    path id_mappings

    output:
    path 'metadata_merged.tsv', emit: meta_merged

    script:
    """
    merge_metadata.py $metadata $id_mappings metadata_merged.tsv
    """
}
