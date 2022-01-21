// process MAKE_KALLISTO_INDEX {
//     tag "$name"
//     label 'process_high'
//     publishDir "${params.outdir}/kallisto_idx", mode: 'copy'
//
//     input:
//     tuple val(name), path(clone_fasta)
//
//     output:
//     tuple val(name), path('*.kidx'), emit: kallisto_idx
//
//     script:
//     """
//     kallisto index -i ${name}.kidx ${clone_fasta}
//     """
// }
//
// process KALLISTO_QUANT {
//     tag "$name"
//     label 'process_high'
//     publishDir "${params.outdir}/kallisto_quant", mode: 'copy'
//
//     input:
//     tuple val(name), path(reads)
//     tuple val(name2), path(idx)
//
//     output:
//     path "kallisto_${name}", emit: kallisto_out_dirs, optional: true
//
//     script:
//     """
//     if [ "$name" == "$name2" ]; then
//         kallisto quant -t $task.cpus --single -i $idx \
//             --fr-stranded --single -l 150 -s 20 -o kallisto_${name} $reads
//     fi
//     """
// }

// create fasta + build index + quantify (single process)
process KALLISTO_QUANT {
    tag "$name"
    label 'process_high'
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    path gpa
    tuple val(name), path(clone_fasta)
    path trimmed_reads

    output:
    path "kallisto_${name}", emit: kallisto_out_dirs, optional: true

    script:
    """
    make_single_clone_fasta.py $clone_fasta $gpa $name
    kallisto index -i ${name}.kidx "${name}_ss.fna"
    kallisto quant -t $task.cpus --single -i ${name}.kidx \
        --fr-stranded --single -l 150 -s 20 -o kallisto_${name} ${name}_trimmed.fq.gz
    """
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
        --ST_file=${params.st_file} \
        --outf=kallisto_merged_counts.tsv

    merge_kallisto_lens.py \
        --gene_presence_absence=$gpa_file \
        --metadata_merged=$meta_merged \
        --ST_file=${params.st_file} \
        --outf=kallisto_merged_lens.tsv
    """
}
