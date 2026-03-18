/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FASTQ         } from '../modules/nf-core/samtools/fastq/main'
include { FASTPLONG              } from '../modules/nf-core/fastplong/main'
include { MULTIQCMAPPER          } from '../modules/local/multiqcmapper/main'

include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'

include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_qc_pipeline'
include { STEP1                  } from '../subworkflows/local/step1'
include { STEP2                  } from '../subworkflows/local/step2'
include { STEP3                  } from '../subworkflows/local/step3'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AQUA {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    // create reference channels ////

    ch_fasta = Channel
        .fromPath(params.fasta, checkIfExists: true)
        .map { fasta -> tuple([id: fasta.getSimpleName()], fasta) }
        .collect()

    ch_fai = params.fasta_fai
        ? Channel.fromPath(params.fasta_fai, checkIfExists: true).map { fai -> tuple([id: fai.getSimpleName()], fai) }.collect()
        : Channel.empty()
    ch_intervals = params.intervals
        ? Channel.fromPath(params.intervals, checkIfExists: true).map { intervals -> tuple([id: intervals.getSimpleName()], intervals) }.collect()
        : Channel.empty()
    ch_gtf = params.gtf
        ? Channel.fromPath(params.gtf, checkIfExists: true).map { gtf -> tuple([id: gtf.getSimpleName()], gtf) }.collect()
        : Channel.empty()
    ch_bed12 = params.bed12
        ? Channel.fromPath(params.bed12, checkIfExists: true).map { bed12 -> tuple([id: bed12.getSimpleName()], bed12) }.collect()
        : Channel.empty()
    ch_refvcf = params.refvcf
        ? Channel.fromPath(params.refvcf, checkIfExists: true).collect()
        : Channel.empty()

    // create empty channels for versions and multiqc files
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_samplesheet.view()
    ch_samplesheet
        .branch { it ->
            def meta = it[0]
            step1: meta.step == 1
            step2: meta.step == 2
            step3: meta.step == 3
            other: false
        }
        .set { samplesheet }

    // STEP 1 // FASTQ FILE ANALYSIS //

    ch_pacbio_step2 = samplesheet.step2.filter { meta, bam, bai -> meta.experiment_method?.toLowerCase() in ["pacbio"] }
    ch_pacbio_step1 = samplesheet.step1.filter { meta, reads -> meta.experiment_method?.toLowerCase() in ["pacbio"] }
    //
    // // Convert read BAMs to FASTQs
    def samtools_fastq_interleave = false
    SAMTOOLS_FASTQ(
        ch_pacbio_step2,
        samtools_fastq_interleave,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)
    //
    // // Run fastplong for pacbio samples
    if (!params.skip_tools?.contains('fastplong')) {
        FASTPLONG(
            SAMTOOLS_FASTQ.out.other.mix(ch_pacbio_step1),
            [],
            false,
            false,
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.json.collect { _meta, json -> json })
        ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.html.collect { _meta, html -> html })
        ch_versions = ch_versions.mix(FASTPLONG.out.versions)
    }

    STEP1(
        samplesheet.step1
    )
    ch_versions = ch_versions.mix(STEP1.out.ch_versions)
    ch_multiqc_files = ch_multiqc_files.mix(STEP1.out.ch_multiqc_files)

    STEP2(
        samplesheet.step2,
        ch_fasta,
        ch_fai,
        ch_intervals,
        ch_refvcf,
    )
    ch_versions = ch_versions.mix(STEP2.out.ch_versions)
    ch_multiqc_files = ch_multiqc_files.mix(STEP2.out.ch_multiqc_files)

    STEP3(
        samplesheet.step3,
        ch_fasta,
    )
    ch_versions = ch_versions.mix(STEP3.out.ch_versions)
    ch_multiqc_files = ch_multiqc_files.mix(STEP3.out.ch_multiqc_files)

    //TODO: ADD A TOOL TO PREDICT CONTAMINATIONS BTW SOMATIC AND NORMAL CELLS

    // such a tool https://github.com/ghga-de/nf-snvcalling/blob/main/bin/PurityReloaded.py can be a nice addition
    // but it is very old using python2.7 plus it is hard-coded towards dkfz pipelines. If we want to use the logic we 
    // have to rewrite it and make it more general. It can be added in the next version of the pipeline if there is a need for it.

    // another tool would be a nice addition is https://github.com/ghga-de/nf-platypusindelcalling/blob/main/modules/local/sample_swap.nf
    // can calculate sample swaps between tumor and normal samples based on vcf files. The tool is old and requires some work to


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'qc_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    if (params.concepts_yaml) {
        MULTIQCMAPPER(
            MULTIQC.out.data,
            file(params.concepts_yaml),
        )
        ch_versions = ch_versions.mix(MULTIQCMAPPER.out.versions)
    }

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
