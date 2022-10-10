process CHECK_META_FILE {
    tag "$metadata"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata

    output:
    path 'metadata_final.tsv', emit: sample_metadata

    script:
    """
    check_metadata.py $metadata metadata_final.tsv
    """
}
