#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
================================================================================
    OrthoQuant: strain-specific analysis of bacterial RNA-Seq data
================================================================================
    Github : [github.com/adamd3/OrthoQuant]
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

if (params.meta_file) {
    ch_samples = file(params.meta_file, checkIfExists: true)
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
include {CHECK_META_FILE} from './modules/metadata'
// include {MAKE_CLONE_FASTA} from './modules/make_clone_fasta'
include {TRIMGALORE} from './modules/trim_reads'
include {KALLISTO_QUANT; MERGE_COUNTS_AND_LENS} from './modules/kallisto'
include {SUBSET_GENES; LENGTH_SCALE_COUNTS; TMM_NORMALISE_COUNTS; DESEQ_NORMALISE_COUNTS} from './modules/normalisation'
include {UMAP_SAMPLES} from './modules/plots'



/*
================================================================================
    Functions
================================================================================
*/
// Function to get list of [ meta, [ fastq_1, fastq2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create sample metadata
    def meta = [:]
    meta.sample_id    = row.dna_sample_id
    meta.paired_end   = row.paired.toBoolean()

    // add path(s) of the fastq file(s) to the metadata
    def seq_meta = []
    if (!file(row.fastq1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq1}"
    }
    if (meta.paired_end) {
        if (!file(row.fastq2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq2}"
        }
        seq_meta = [ meta, [ file(row.fastq1), file(row.fastq2) ], file(row.fasta) ]
    } else {
        seq_meta = [ meta, file(row.fastq1), file(row.fasta) ]
    }
    return seq_meta
}


/*
================================================================================
    Main workflow
================================================================================
*/
workflow {

    /*
     *  Create channels for input files
     */

    CHECK_META_FILE (
        ch_samples
    )
    ch_metadata = CHECK_META_FILE.out.sample_metadata


    // Channel
    //     .fromPath(params.meta_file)
    //     .splitCsv(header:true, sep:'\t')
    //     // .map { row -> [ row.dna_sample_id, [ file(row.fastq, checkIfExists: true) ] ] }
    //     .map { create_fastq_channel(it) }
    //     .set { ch_raw_reads_trimgalore }

    // Channel
    //     .fromPath(params.meta_file)
    //     .splitCsv(header:true, sep:'\t')
    //     .map { row -> [ row.dna_sample_id, [ file(row.fasta, checkIfExists: true) ] ] }
    //     .set { ch_clone_fasta_init }

    ch_metadata
        .splitCsv(header: true, sep:'\t')
        .map { create_fastq_channel(it) }
        .set { ch_raw_reads_trimgalore }


    /*
     *  Get the subset of genes to be included in the analysis
     */
    SUBSET_GENES (
        ch_gpa_file,
        ch_metadata,
        params.perc
    )
    ch_gene_subset = SUBSET_GENES.out.gene_subset

    /*
     *  Trim + pseudo-align reads
     */
    if (params.skip_trimming) {

        KALLISTO_QUANT (
            ch_gpa_file,
            ch_raw_reads_trimgalore,
            params.strandedness
        )
        ch_kallisto_out_dirs = KALLISTO_QUANT.out.kallisto_out_dirs.collect()

    } else {

        TRIMGALORE (
            ch_raw_reads_trimgalore
        )
        // ch_trimmed_reads = TRIMGALORE.out.trimmed_reads.collect()
        ch_trimmed_reads = TRIMGALORE.out.trimmed_reads
        ch_trimgalore_results_mqc = TRIMGALORE.out.trimgalore_results_mqc
        ch_trimgalore_fastqc_reports_mqc = TRIMGALORE.out.trimgalore_fastqc_reports_mqc

        KALLISTO_QUANT (
            ch_gpa_file,
            ch_trimmed_reads,
            params.strandedness
        )
        ch_kallisto_out_dirs = KALLISTO_QUANT.out.kallisto_out_dirs.collect()
    }

    /*
     *  Merge counts
     */
    MERGE_COUNTS_AND_LENS (
        ch_gpa_file,
        ch_kallisto_out_dirs,
        ch_metadata
    )
    ch_kallisto_merged_out = MERGE_COUNTS_AND_LENS.out.kallisto_merged_out

    /*
     *  Scale counts to median gene length across strains
     */
    LENGTH_SCALE_COUNTS (
        ch_kallisto_merged_out,
        ch_gene_subset
    )
    ch_scaled_counts = LENGTH_SCALE_COUNTS.out.scaled_counts

    /*
     *  Get size-factor-scaled counts
     */
    if (params.norm_method == 'DESeq') {
        DESEQ_NORMALISE_COUNTS (
            ch_kallisto_merged_out,
            ch_gene_subset
        )
        ch_norm_counts = DESEQ_NORMALISE_COUNTS.out.norm_counts
    } else if (params.norm_method == 'TMM') {
        TMM_NORMALISE_COUNTS (
            ch_kallisto_merged_out,
            ch_gene_subset
        )
        ch_norm_counts = TMM_NORMALISE_COUNTS.out.norm_counts
    }
    // NB the scaled counts are log-transformed by default; the RPKM counts are not

    /*
     *  UMAP of samples
     */
    UMAP_SAMPLES (
        ch_norm_counts,
        ch_metadata,
        params.group
    )
    ch_umap_out = UMAP_SAMPLES.out.umap_out

}

/*
================================================================================
    Completion summary
================================================================================
*/

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
      nextflow run OrthoQuant --meta_file [file] --gpa_file [gene_presence_absence.csv] --group [str] --perc [str] --norm_method [str] -profile conda

    Mandatory arguments:
      --meta_file [file]              Path to file containing sample metadata.
      --gpa_file [file]               Path to file containing gene presence/absence per strain (from Panaroo output).
      --group [str]                   Metadata column to use for grouping strains. Default = majority_ST.
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: conda, docker

    Other options:
      --strandedness [str]            Strandedness of the reads. 0 = unstranded, 1 = forward-stranded, 2 = reverse-stranded. Default = 2.
      --fragment_len [str]            Estimated average fragment length for kallisto transcript quantification (only required for single-end reads). Default = 150.
      --fragment_sd [str]             Estimated standard deviation of fragment length for kallisto transcript quantification (only required for single-end reads). Default = 20.
      --perc [str]                    Minimum percent of strains containing a gene for defining the core gene set. Default = 99.
      --norm_method [str]             How to perform size-factor scaling of counts for normalisation. Available options: DESeq (default), TMM.
      --skip_trimming [bool]          Do not trim adaptors from FastQ files.
      --outdir [file]                 The output directory where the results will be saved (default: './results').
      --cachedir [str]                Defines the path where Singularity images are cached [default: $params.cachedir]
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    """.stripIndent()
}
