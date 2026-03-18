#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_qc_pipeline'


include { AQUA                    } from './workflows/aqua'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_qc_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_qc_pipeline'
include { METADATA_TO_SAMPLESHEET } from './modules/local/metadata_to_samplesheet'
params.fasta     = getGenomeAttribute('fasta')
params.fasta_fai = getGenomeAttribute('fasta_fai')
params.intervals = getGenomeAttribute('intervals')
params.gtf       = getGenomeAttribute('gtf')
params.bed12     = getGenomeAttribute('bed12')
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    if (!params.input && !params.metadata) {
        error("Please provide a valid --input file ending in .csv or --metadata ending in .json")
    }

    def ch_input

    if (params.metadata) {
        METADATA_TO_SAMPLESHEET(
            Channel.fromPath(params.metadata, checkIfExists: true)
        )
        ch_input = METADATA_TO_SAMPLESHEET.out.samplesheet
    }
    else {
        ch_input = Channel.fromPath(params.input, checkIfExists: true)
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        ch_input,
    )

    //
    // WORKFLOW: Run main workflow
    //
    GHGA_AQUA(
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        GHGA_AQUA.out.multiqc_report,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GHGA_AQUA {
    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    AQUA(
        samplesheet
    )

    emit:
    multiqc_report = AQUA.out.multiqc_report // channel: /path/to/multiqc_report.html
}
