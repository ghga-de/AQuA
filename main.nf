#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_qc_pipeline'
params.fasta     = getGenomeAttribute('fasta')
params.fasta_fai = getGenomeAttribute('fasta_fai')
params.intervals = getGenomeAttribute('intervals')
params.gtf       = getGenomeAttribute('gtf')
params.bed12     = getGenomeAttribute('bed12')


include { QCMETRICS               } from './workflows/qcmetrics'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_qc_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_qc_pipeline'
include { METADATA_TO_SAMPLESHEET } from './modules/local/metadata_to_samplesheet'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GHGA_QCMETRICS {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    QCMETRICS (
        samplesheet
    )
    emit:
    multiqc_report = QCMETRICS.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    if (!params.input && !params.metadata) {
        error("Please provide a valid --input file ending in .csv or --metadata ending in .json")
    }

    def ch_input

    if (params.metadata) {
        METADATA_TO_SAMPLESHEET (
            Channel.fromPath(params.metadata, checkIfExists: true)
        )
        ch_input = METADATA_TO_SAMPLESHEET.out.samplesheet
    } else {
        ch_input = Channel.fromPath(params.input, checkIfExists: true)
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        ch_input
    )

    //
    // WORKFLOW: Run main workflow
    //
    GHGA_QCMETRICS (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        GHGA_QCMETRICS.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/