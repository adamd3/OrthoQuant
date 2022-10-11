process KALLISTO_QUANT {
    tag "$meta.sample_id"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    path gpa
    // tuple val(name), path(clone_fasta)
    tuple val(meta), path(reads), path(fasta)

    output:
    path "kallisto_${meta.sample_id}", emit: kallisto_out_dirs

    script:
    def name = task.ext.prefix ?: "${meta.sample_id}"

    if (meta.paired_end) {
        // if trimming has not been performed, must symlink to match expected file names
        // kallisto params -l, -s are estimated from paired end data, but are required when using --single
        """
        [ ! -f  ${name}_1_val_1.fq.gz ] && ln -s ${reads[0]} ${name}_1_val_1.fq.gz
        [ ! -f  ${name}_2_val_2.fq.gz ] && ln -s ${reads[1]} ${name}_2_val_2.fq.gz
        kallisto index -i ${name}.kidx $fasta
        kallisto quant -t $task.cpus -i ${name}.kidx \
            --fr-stranded -o kallisto_${name} ${name}_1_val_1.fq.gz ${name}_2_val_2.fq.gz
        """
    } else {
        """
        [ ! -f  ${name}_trimmed.fq.gz ] && ln -s $reads ${name}_trimmed.fq.gz
        kallisto index -i ${name}.kidx $fasta
        kallisto quant -t $task.cpus -i ${name}.kidx \
            --fr-stranded --single -l 150 -s 20 -o kallisto_${name} ${name}_trimmed.fq.gz
        """
    }
}


process MERGE_COUNTS_AND_LENS {
    tag "merge_counts_and_lens"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    path gpa_file
    path kallisto_dirs
    path meta_merged

    output:
    tuple path('kallisto_merged_counts.tsv'), path('kallisto_merged_lens.tsv'), emit: kallisto_merged_out

    script:
    """
    merge_kallisto_counts.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --outf=kallisto_merged_counts.tsv

    merge_kallisto_lens.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --outf=kallisto_merged_lens.tsv
    """
}
