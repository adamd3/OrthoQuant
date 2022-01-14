#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
================================================================================
    strain_seq: strain-specific analysis of bacterial RNA-Seq data
================================================================================
    Github : [github.com/adamd3/strain_seq]
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
    Validate inputs and create channels for files
================================================================================
*/

// if (params.data_dir) {
//     ch_data_dir = file(params.data_dir, checkIfExists: true)
// } else { exit 1, 'Data directory not specified!' }

if (params.meta_file) {
    ch_metadata = file(params.meta_file, checkIfExists: true)
} else { exit 1, 'Metadata file not specified!' }

// if (params.multifasta_file) {
//     ch_multifasta_file = file(params.multifasta_file, checkIfExists: true)
// } else { exit 1, 'Multi-fasta file not specified!' }

// if (params.faidx_file) {
//     ch_faidx_file = file(params.faidx_file, checkIfExists: true)
// } else { exit 1, 'Index for multi-fasta file not specified!' }

if (params.gpa_file) {
    ch_gpa_file = file(params.gpa_file, checkIfExists: true)
} else { exit 1, 'Gene presence/absence file not specified!' }




/*
================================================================================
    Modules
================================================================================
*/
include {MERGE_METADATA} from './modules/merge_metadata'
include {MAKE_CLONE_FASTA} from './modules/make_clone_fasta'
include {TRIMGALORE} from './modules/trim_reads'
include {MAKE_KALLISTO_INDEX; KALLISTO_QUANT} from './modules/kallisto'


/*
================================================================================
    Main workflow
================================================================================
*/

workflow {
    
    /*
     * Merge metadata from WGS and RNA-Seq
     */
    MERGE_METADATA (
        ch_metadata
    )
    ch_meta_merged = MERGE_METADATA.out.meta_merged

    /*
     *  Create channels for input files
     */
    ch_meta_merged
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.sample_id, [ file(row.fastq, checkIfExists: true) ] ] }
        .set { ch_raw_reads_trimgalore }

    ch_meta_merged
        .splitCsv(header:true, sep:'\t')
        .map { row -> [ row.sample_id, [ file(row.fasta, checkIfExists: false) ] ] }
        .set { ch_clone_fasta_init }

    // extract the sample IDs only:
    // ch_meta_merged
    //     .splitCsv(header: true, sep:'\t')
    //     .map { row -> row.sample_id }
    //     .set { ch_clone_ids }

    /*
     *  Create strain-specific FASTA files
     */
    MAKE_CLONE_FASTA (
        ch_gpa_file,
        ch_clone_fasta_init
    )
    ch_clone_fasta = MAKE_CLONE_FASTA.out.clone_fasta

    /*
     *  Create a Kallisto index for each strain
     */
    MAKE_KALLISTO_INDEX (
        ch_clone_fasta
    )
    ch_kallisto_idx = MAKE_KALLISTO_INDEX.out.clone_fasta

    /*
     *  Trim reads
     */
    if (params.skip_trimming) {
        ch_trimmed_reads = ch_raw_reads_trimgalore
        ch_trimgalore_results_mqc = Channel.empty()
        ch_trimgalore_fastqc_reports_mqc = Channel.empty()
    } else {
        TRIMGALORE (
            ch_raw_reads_trimgalore
        )
        ch_trimmed_reads = TRIMGALORE.out.trimmed_reads
        ch_trimgalore_results_mqc = TRIMGALORE.out.trimgalore_results_mqc
        ch_trimgalore_fastqc_reports_mqc = TRIMGALORE.out.trimgalore_fastqc_reports_mqc
    }

    /*
     *  Create a Kallisto index for each strain
     */
    KALLISTO_QUANT (
        ch_trimmed_reads,
        ch_kallisto_idx
    )
    // NOTE: this is a directory containing the kallisto results
    ch_kallisto_out = KALLISTO_QUANT.out.kallisto_out

}


/*
================================================================================
    COMPLETION EMAIL AND SUMMARY
================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
    }
    NfcoreTemplate.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
}

c_green = "\033[0;32m";
c_reset = "\033[0m"

workflow.onComplete {
    log.info"""
    Execution status: ${ workflow.success ? 'OK' : 'failed' }
    ${c_green}Results are reported here: $params.outdir${c_reset}
    ￼""".stripIndent()
}


def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
      nextflow run strain_seq --data_dir [dir] --meta_file [file] --gpa_file [gene_presence_absence.csv] -profile docker

    Mandatory arguments:
      --data_dir [file]               Path to directory containing FastQ files retrieved using the nf-core/fetchngs pipeline.
      --meta_file [file]              Path to file containing sample metadata.
      --gpa_file [file]               Path to file containing gene presence/absence per strain (from Panaroo output).
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: conda, docker

    Other options:
      --outdir [file]                 The output directory where the results will be saved (Default: './results').
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}
