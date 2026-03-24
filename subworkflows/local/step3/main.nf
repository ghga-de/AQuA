//
// STEP3: ANALYSIS TOOLS FOR VCF/BCF FILES
//

include { TABIX_BGZIP    } from '../../../modules/nf-core/tabix/bgzip/main'
include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../../modules/nf-core/bcftools/stats/main'

workflow STEP3 {
    take:
    samplesheet // channel: [val(meta), vcf]
    ch_fasta

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    samplesheet.branch { meta, vcf ->
            compressed: vcf.getName().endsWith('.gz')
            uncompressed: true
        }.set { ch_branched_vcfs }

    TABIX_BGZIP(
        ch_branched_vcfs.uncompressed
    )
    ch_ready_vcfs = ch_branched_vcfs.compressed.mix(TABIX_BGZIP.out.output)

    TABIX_TABIX(
        ch_ready_vcfs
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    regions = [[], []]
    targets = [[], []]
    samples = [[], []]
    exons = [[], []]

    if (!params.skip_tools?.contains('bcftools_stats')) {
        BCFTOOLS_STATS(
            samplesheet.join(TABIX_TABIX.out.tbi),
            regions,
            targets,
            samples,
            exons,
            ch_fasta,
        )
        ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    }

    emit:
    ch_versions      = ch_versions
    ch_multiqc_files = ch_multiqc_files
}
