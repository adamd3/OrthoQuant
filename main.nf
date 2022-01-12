#!/usr/bin/env nextflow

/*
================================================================================
    strain_seq: strain-specific analysis of bacterial RNA-Seq data
================================================================================
    Github : [github.com/adamd3/strain_seq]
*/
def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run strain_seq --data_dir [dir] --meta_file [file] --multifasta_file [file] --gpa_file [gene_presence_absence.csv] -profile docker

    Mandatory arguments:
      --data_dir [file]               Path to directory containing FastQ files retrieved using the nf-core/fetchngs pipeline.
      --meta_file [file]              Path to file containing sample metadata.
      --multifasta_file [file]        Path to multi-fasta file containing all strain-specific gene sequences combined.
      --gpa_file [file]               Path to file containing gene presence/absence per strain (from Panaroo output).
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: conda, docker

    Other options:
      --outdir [file]                 The output directory where the results will be saved (Default: './results').
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}


/*
================================================================================
    Configuration
================================================================================
*/
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}


/*
================================================================================
    Validate inputs
================================================================================
*/

if (params.data_dir) {
    ch_data_dir = file(params.data_dir, checkIfExists: true)
} else { exit 1, 'Data directory not specified!' }

if (params.meta_file) {
    ch_metadata = file(params.meta_file, checkIfExists: true)
} else { exit 1, 'Metadata file not specified!' }

if (params.multifasta_file) {
    ch_multifasta_file = file(params.multifasta_file, checkIfExists: true)
} else { exit 1, 'Multi-fasta file not specified!' }

if (params.gpa_file) {
    ch_gpa_file = file(params.gpa_file, checkIfExists: true)
} else { exit 1, 'Gene presence/absence file not specified!' }




/*
================================================================================
    Pre-processing steps
================================================================================
*/

/*
 * Merge metadata from WGS and RNA-Seq
 */

process MERGE_METADATA {
    tag "$metadata"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata from ch_metadata
    path data_dir from ch_data_dir

    output:
    path 'metadata_merged.tsv' into ch_meta_merged

    script:
    """
    merge_metadata.py $metadata $data_dir metadata_merged.tsv
    """
}




/*
 *  Create strain-specific FASTA files
 */

process MAKE_CLONE_FASTA {
     tag "$multifasta"
     publishDir "${params.outdir}/clone_fasta", mode: 'copy'

     input:
     path multifasta from ch_multifasta_file
     path gpa from ch_gpa_file

     output:
     path '*.fna' into ch_clone_fasta

     script:
     """
     make_clone_fasta.py $multifasta $gpa ./
     """
 }


/*
 *  Create a Kallisto index for each strain
 */

process MAKE_KALLISTO_INDEX {
    tag "$clone_fasta"
    publishDir "${params.outdir}/kallisto_idx", mode: 'copy'

    input:
    path clone_fasta from ch_clone_fasta

    output:
    path '*.kidx' into ch_kallisto_idx

    script:
    """
    kallisto index -i ${clone_fasta.baseName}.kidx ${clone_fasta}
    """
}





/*
================================================================================
    Create channel for input FastQ files
================================================================================
*/

ch_meta_merged
    .splitCsv(header:true, sep:'\t')
    .map { row -> [ row.sample_id, [ file(row.fastq, checkIfExists: true) ] ] }
    .set { ch_raw_reads_trimgalore }



/*
================================================================================
    Trim reads
================================================================================
*/

if (params.skip_trimming) {
    ch_trimmed_reads = ch_raw_reads_trimgalore
    ch_trimgalore_results_mqc = Channel.empty()
    ch_trimgalore_fastqc_reports_mqc = Channel.empty()
} else {
    process TRIMGALORE {
        tag "$name"
        label 'process_high'
        publishDir "${params.outdir}/trim_galore", mode: 'copy',
            saveAs: { filename ->
                          if (filename.endsWith('.html')) "fastqc/$filename"
                          else if (filename.endsWith('.zip')) "fastqc/zips/$filename"
                          else if (filename.endsWith('trimming_report.txt')) "logs/$filename"
                          else params.save_trimmed ? filename : null
                    }

        input:
        tuple val(name), path(reads) from ch_raw_reads_trimgalore

        output:
        tuple val(name), path('*.fq.gz') into ch_trimmed_reads
        path '*.txt' into ch_trimgalore_results_mqc
        path '*.{zip,html}' into ch_trimgalore_fastqc_reports_mqc

        script:
        // Calculate number of --cores for TrimGalore based on value of task.cpus
        // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
        // See: https://github.com/nf-core/atacseq/pull/65
        // (Max no cores = 4, since there are diminishing returns beyond this)
        def cores = 1
        if (task.cpus) {
            cores = (task.cpus as int) - 3
            if (cores < 1) cores = 1
            if (cores > 4) cores = 4
        }

        // Add symlinks to original fastqs for consistent naming in MultiQC
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        trim_galore --cores $cores --fastqc --gzip ${name}.fastq.gz
        """
    }
}


/*
================================================================================
    Quantify strain-specific gene sets with kallisto
================================================================================
*/

process KALLISTO_QUANT {
    tag "$name"
    publishDir "${params.outdir}/kallisto_quant", mode: 'copy'

    input:
    tuple val(name), path(reads) from ch_trimmed_reads
    path index from ch_kallisto_idx

    output:
    path "$name" into ch_kallisto_out

    script:
    """
    kallisto quant -t $task.cpus --single -i $index \
        --fr-stranded --single -l 150 -s 20 -o $name $reads
    """
}
