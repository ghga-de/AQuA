//
// STEP1: ANALYSIS TOOLS FOR RAW FILES
//

include { FASTQC      } from '../../../modules/nf-core/fastqc/main'
include { FASTP       } from '../../../modules/nf-core/fastp/main'
include { FASTPLONG   } from '../../../modules/nf-core/fastplong/main'
include { SEQFU_STATS } from '../../../modules/nf-core/seqfu/stats/main'
include { NANOPLOT as NANOPLOT_FASTQ } from '../../../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_SIGNAL} from '../../../modules/nf-core/nanoplot/main'

workflow STEP1 {
    take:
    samplesheet // channel: [val(meta), fastq1, fastq2]
    tools       // list of tools to run
    nanopore_summary // channel: [val(meta), sequencing_summary.txt]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    samplesheet.branch { meta, reads ->
        short_ch:   !(meta.experiment_method?.toLowerCase() in ['nanopore', 'pacbio'])
        long_ch:    true
    }.set { inputs }

    //  Runs FASTQC
    if (tools.contains('fastqc')) {
        FASTQC(
            samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    //  Runs FASTP for all except pacbio, rna and methylseq samples, as FASTP does not support these data types
    if (tools.contains('fastp')) {
        save_trimmed_fail = false
        save_merged = false

        FASTP(
            inputs.short_ch,
            [],
            false,
            save_trimmed_fail,
            save_merged,
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { _meta, json -> json })
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { _meta, html -> html })
        ch_versions = ch_versions.mix(FASTP.out.versions)

    }
    // filter out fast5 and pod5 for now
    inputs.long_ch.branch { meta, reads ->
        def file_list = reads instanceof Collection ? reads : [reads]
        signal: file_list.any { it.name.endsWith('.fast5') || it.name.endsWith('.pod5') }  
        fastq: true 
    }.set { split_long_reads }

    // Run fastplong for long samples
    if (tools.contains('fastplong')) {
        FASTPLONG(
            split_long_reads.fastq,
            [],
            false,
            false,
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.json.collect { _meta, json -> json })
        ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.html.collect { _meta, html -> html })
        ch_versions = ch_versions.mix(FASTPLONG.out.versions)
    }

    split_long_reads.fastq.filter{ meta, file ->
        !(meta.experiment_method?.toLowerCase() in ['pacbio'])
    }.set{ch_no_pacbio}
    
    // Runs SEQFU_STATS for all non-PacBio samples, as SEQFU_STATS does not support PacBio data
    if (tools.contains('seqfu')) {
        SEQFU_STATS(
            ch_no_pacbio
        )
        ch_multiqc_files = ch_multiqc_files.mix(SEQFU_STATS.out.multiqc.collect { _meta, file -> file })
        ch_versions = ch_versions.mix(SEQFU_STATS.out.versions)
    }

    // Runs NANOPLOT for Nanopore data
    if (tools.contains('nanoplot')) {
        NANOPLOT_FASTQ(
            split_long_reads.fastq
        )
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_FASTQ.out.txt.collect { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_FASTQ.out.log.collect { it[1] })
        ch_versions = ch_versions.mix(NANOPLOT_FASTQ.out.versions)
    }

    // nanoplot is the only option for signal files
    NANOPLOT_SIGNAL(
        nanopore_summary
    )
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_SIGNAL.out.txt.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_SIGNAL.out.log.collect { it[1] })
    ch_versions = ch_versions.mix(NANOPLOT_SIGNAL.out.versions)

    emit:
    ch_versions      = ch_versions
    ch_multiqc_files = ch_multiqc_files
}
