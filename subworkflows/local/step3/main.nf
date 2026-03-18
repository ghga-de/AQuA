//
// STEP3: ANALYSIS TOOLS FOR VCF/BCF FILES
//

include { TABIX_TABIX    } from '../../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_STATS } from '../../../modules/nf-core/bcftools/stats/main'

workflow STEP3 {
    take:
    samplesheet // channel: [val(meta), vcf]
    ch_fasta

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Runs BCFTOOLS_STATS 
    // TODO add support for vcf and bcf-bcf.gz
    TABIX_TABIX(
        samplesheet
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
