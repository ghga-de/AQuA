//
// STEP2: ANALYSIS TOOLS FOR BAM/CRAM FILES
//

include { SAMTOOLS_FAIDX                } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS                } from '../../../modules/nf-core/samtools/stats/main'
include { MOSDEPTH                      } from '../../../modules/nf-core/mosdepth/main'
include { RSEQC_BAMSTAT                 } from '../../../modules/nf-core/rseqc/bamstat/main'
include { PRESEQ_LCEXTRAP               } from '../../../modules/nf-core/preseq/lcextrap/main'
include { NGSBITS_SAMPLEGENDER          } from '../../../modules/nf-core/ngsbits/samplegender/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { VERIFYBAMID_VERIFYBAMID       } from '../../../modules/nf-core/verifybamid/verifybamid/main'

workflow STEP2 {
    take:
    samplesheet // channel: [val(meta), bam, cram]
    ch_fasta
    ch_fai
    ch_intervals
    ch_refvcf

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    if (!params.fasta_fai) {
        // this can be replaced if more than once tool will use fai
        SAMTOOLS_FAIDX(
            ch_fasta,
            [[], []],
            false,
        )
        ch_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    // Filter out PacBio samples for coverage and general metrics tools
    ch_bam_for_metrics = samplesheet.filter { meta, bam, index ->
        !(meta.experiment_method?.toLowerCase() in ["pacbio"])
    }

    // Branch samples: those with index and those without
    ch_bam_for_metrics
        .branch { meta, bam, index ->
            has_index: index
            no_index: !index
        }
        .set { ch_bam_indexed_split }

    // Run indexing for samples missing an index
    SAMTOOLS_INDEX(
        ch_bam_indexed_split.no_index.map { meta, bam, _index -> tuple(meta, bam) }
    )

    // Merge the newly created indexes with the BAM files
    ch_newly_indexed = ch_bam_indexed_split.no_index
        .join(SAMTOOLS_INDEX.out.bai.mix(SAMTOOLS_INDEX.out.csi), remainder: false)
        .map { meta, bam, old_empty, new_index -> [meta, bam, new_index] }

    // Combine everything back into one channel for downstream tools
    ch_final_bam_indexed = ch_bam_indexed_split.has_index.mix(ch_newly_indexed)

    // Runs MOSDEPTH
    if (!params.skip_tools?.contains('mosdepth')) {

        ch_final_bam_indexed
            .branch { meta, file, index ->
                targeted: meta.experiment_method?.toLowerCase() in ["wes", "tes", "atacseq", "chipseq"]
                other: true
            }
            .set { mosdepth_split }

        ch_mosdepth_targeted = mosdepth_split.targeted
            .combine(ch_intervals)
            .map { meta, file, index, _meta2, intervals -> tuple(meta, file, index, intervals) }

        ch_mosdepth_other = mosdepth_split.other.map { meta, file, index -> tuple(meta, file, index, []) }

        ch_mosdepth = ch_mosdepth_targeted.mix(ch_mosdepth_other)

        MOSDEPTH(
            ch_mosdepth,
            ch_fasta,
        )

        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.map { _meta, file -> file }.collect())
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    }

    // Runs SAMTOOLS_STATS
    if (!params.skip_tools?.contains('samtools_stats')) {
        SAMTOOLS_STATS(
            ch_final_bam_indexed,
            ch_fasta,
        )
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    }

    // Runs PICARD_COLLECTMULTIPLEMETRICS
    // Run Picard CollectInsertSizeMetrics
    // if rna run picard_collectrnaseqmetrics
    // if small rna run mirdeep2
    // if atac seq run ATACseqQC
    // if methlation run Bismark or Qualimap
    // if chip seq Phantompeakqualtools
    if (!params.skip_tools?.contains('picard_collectmultiplemetrics')) {
        PICARD_COLLECTMULTIPLEMETRICS(
            ch_final_bam_indexed,
            ch_fasta,
            ch_fai,
        )
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
    }

    if (params.run_contamination_estimation && ch_refvcf) {
        // Runs VERIFYBAMID_VERIFYBAMID if samplesheet has bam/cram files and if reference vcf is given
        if (!params.skip_tools?.contains('verifybamid')) {
            VERIFYBAMID_VERIFYBAMID(
                ch_final_bam_indexed,
                ch_refvcf,
            )
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.selfsm.map { _meta, file -> file }.collect())
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.depthsm.map { _meta, file -> file }.collect())
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.selfrg.map { _meta, file -> file }.collect())
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.depthrg.map { _meta, file -> file }.collect())
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.bestsm.map { _meta, file -> file }.collect())
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.log.map { _meta, file -> file }.collect())
            ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID.out.versions)
        }
    }


    // Runs PRESEQ_LCEXTRAP if samplesheet has bam/cram files
    if (!params.skip_tools?.contains('preseq')) {
        PRESEQ_LCEXTRAP(
            ch_final_bam_indexed
        )
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.map { _meta, file -> file }.collect())
    }

    // Runs RSEQC_BAMSTAT if samplesheet has bam/cram files and method is rna
    ch_rseqc_in = ch_final_bam_indexed.filter { meta, bam, bai -> meta.experiment_method?.toLowerCase()?.contains("rna") }
    if (!params.skip_tools?.contains('rseqc')) {
        RSEQC_BAMSTAT(
            ch_rseqc_in.map { meta, bam, _bai -> tuple(meta, bam) }
        )
        ch_multiqc_files = ch_multiqc_files.mix(RSEQC_BAMSTAT.out.txt.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)
    }

    // Predict sex of the samples if not given and if samplesheet has bam/cram files
    ch_gender_in = ch_final_bam_indexed.filter { meta, bam, bai -> meta.experiment_method?.toLowerCase()?.contains("wgs") }
    if (params.predict_sex) {
        NGSBITS_SAMPLEGENDER(
            ch_gender_in,
            ch_fasta,
            ch_fai,
            params.samplegender_method ?: 'xy',
        )
        ch_multiqc_files = ch_multiqc_files.mix(NGSBITS_SAMPLEGENDER.out.tsv.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions)
    }

    emit:
    ch_versions      = ch_versions
    ch_multiqc_files = ch_multiqc_files
}
