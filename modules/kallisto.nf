process MAKE_KALLISTO_INDEX {
    tag "$name"
    publishDir "${params.outdir}/kallisto_idx", mode: 'copy'

    input:
    tuple val(name), path(clone_fasta)

    output:
    tuple val(name), path('*.kidx'), emit: kallisto_idx

    script:
    """
    kallisto index -i ${name}.kidx ${clone_fasta}
    """
}

process KALLISTO_QUANT {
    tag "$name"
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    tuple val(name), path(reads)
    tuple val(name2), path(idx)

    output:
    path "$name", emit: kallisto_out

    script:
    """
    kallisto quant -t $task.cpus --single -i $idx \
        --fr-stranded --single -l 150 -s 20 -o $name $reads
    """
}
