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

process SUBSET_GENES {
    tag "$st_file"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path gpa_file
    path meta_merged
    path st_file
    val perc

    output:
    path 'gene_set_ST.tsv', emit: meta_merged

    script:
    """
    find_core_genome.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --ST_file=$st_file --perc=$perc \
        --outf=gene_set_ST.tsv
    """
}
