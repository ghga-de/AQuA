//
// STEP1: ANALYSIS TOOLS FOR RAW FILES
//

include { FASTQC      } from '../../../modules/nf-core/fastqc/main'
include { FASTP       } from '../../../modules/nf-core/fastp/main'
include { NANOPLOT    } from '../../../modules/nf-core/nanoplot/main'
include { SEQFU_STATS } from '../../../modules/nf-core/seqfu/stats/main'

workflow STEP1 {
    take:
    samplesheet // channel: [val(meta), fastq1, fastq2]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_no_pacbio = samplesheet.filter { meta, reads ->
        !(meta.experiment_method?.toLowerCase() in ["pacbio", "rna", "methylseq"])
    }

    ch_nanoplot_in = samplesheet.filter { meta, reads ->
        meta.experiment_method?.toLowerCase() in ["nanopore"]
    }

    //  Runs FASTQC
    if (!params.skip_tools?.contains('fastqc')) {
        FASTQC(
            samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    //  Runs FASTP for all exept pacbio, rna and methylseq samples, as FASTP does not support these data types
    if (!params.skip_tools?.contains('fastp')) {
        save_trimmed_fail = false
        save_merged = false
        FASTP(
            ch_no_pacbio,
            [],
            false,
            save_trimmed_fail,
            save_merged,
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { _meta, json -> json })
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { _meta, html -> html })
        ch_versions = ch_versions.mix(FASTP.out.versions)
    }

    // Runs SEQFU_STATS for all non-PacBio samples, as SEQFU_STATS does not support PacBio data
    if (!params.skip_tools?.contains('seqfu')) {
        SEQFU_STATS(
            ch_no_pacbio
        )
        ch_multiqc_files = ch_multiqc_files.mix(SEQFU_STATS.out.multiqc.collect { _meta, file -> file })
        ch_versions = ch_versions.mix(SEQFU_STATS.out.versions)
    }

    // Runs NANOPLOT for Nanopore data
    if (!params.skip_tools?.contains('nanoplot')) {
        NANOPLOT(
            ch_nanoplot_in
        )
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.log.collect { it[1] })
        ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    }

    emit:
    ch_versions      = ch_versions
    ch_multiqc_files = ch_multiqc_files
}
