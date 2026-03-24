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
include { ATAQV_ATAQV                   } from '../../../modules/nf-core/ataqv/ataqv/main'
include { PHANTOMPEAKQUALTOOLS          } from '../../../modules/nf-core/phantompeakqualtools/main'
include { CRAMINO                       } from '../../../modules/nf-core/cramino/main'
include { METHYLDACKEL_MBIAS            } from '../../../modules/nf-core/methyldackel/mbias/main'

workflow STEP2 {
    take:
    samplesheet  
    ch_fasta
    ch_fai
    ch_intervals
    ch_refvcf

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (!params.fasta_fai){
        SAMTOOLS_FAIDX(
            ch_fasta,
            [[],[]],
            false
        )
        ch_fai = SAMTOOLS_FAIDX.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    samplesheet.branch { meta, bam, index ->
        has_index: index
        no_index: !index
    }.set { ch_bam_indexed_split }

    SAMTOOLS_INDEX ( 
        ch_bam_indexed_split.no_index.map { meta, bam, _index -> tuple(meta, bam) }, 
    )

    ch_newly_indexed = ch_bam_indexed_split.no_index
        .join(SAMTOOLS_INDEX.out.bai.mix(SAMTOOLS_INDEX.out.csi), remainder: false)
        .map { meta, bam, old_empty, new_index -> [meta, bam, new_index] }

    ch_final_bam_indexed = ch_bam_indexed_split.has_index.mix(ch_newly_indexed)

    // ========================================================
    // MASTER BRANCHING BY EXPERIMENT METHOD
    // ========================================================

    ch_final_bam_indexed.branch { meta, bam, bai ->
        wgs:    meta.experiment_method?.toLowerCase() in ['wgs']
        wes:    meta.experiment_method?.toLowerCase() in ['wxs', 'wcs', 'wes', 'tes']
        rna:    meta.experiment_method?.toLowerCase() in ['rna', 'total_rna', 'm_rna', 'nc_rna']
        smrna:  meta.experiment_method?.toLowerCase() in ['mi_rna', 'smrna']
        atac:   meta.experiment_method?.toLowerCase() in ['atac', 'atacseq']
        meth:   meta.experiment_method?.toLowerCase() in ['methylation', 'methylseq']
        chip:   meta.experiment_method?.toLowerCase() in ['chip_seq', 'chipseq', 'chip']
        cfdna:  meta.experiment_method?.toLowerCase() in ['cfdna', 'other']
        longread: meta.experiment_method?.toLowerCase() in ['nanopore', 'pacbio']
    }.set { ch_assay_split }


    if (!params.skip_tools?.contains('samtools_stats')) {
        SAMTOOLS_STATS(
            ch_final_bam_indexed,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    }

    // Picard Multiple Metrics: Best for standard DNA fragmentation
    ch_picard_multi_in = ch_assay_split.wgs.mix(ch_assay_split.wes, ch_assay_split.cfdna)
    if (!params.skip_tools?.contains('picard_collectmultiplemetrics')) {
        PICARD_COLLECTMULTIPLEMETRICS(
            ch_picard_multi_in,
            ch_fasta,
            ch_fai
        )
        ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTMULTIPLEMETRICS.out.metrics.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions)
    }

    // Mosdepth: WGS needs no intervals while targeted methods need intervals
    if (!params.skip_tools?.contains('mosdepth')) {
        
        // Prepare targeted samples (WES, ATAC, CHIP, Methylation) - combine with intervals if provided, if not add empty list as placeholder
        ch_mosdepth_targeted = ch_assay_split.wes.mix(ch_assay_split.atac, ch_assay_split.chip, ch_assay_split.meth)
            .combine(ch_intervals.ifEmpty([[:], []])) 
            .map { meta, file, index, meta_int, intervals -> 
                def final_intervals = intervals instanceof List ? [] : intervals
                return tuple(meta, file, index, final_intervals) 
            }
        
        // Prepare other samples (WGS, cfDNA)
        ch_mosdepth_other = ch_assay_split.wgs.mix(ch_assay_split.cfdna, ch_assay_split.longread)
            .map { meta, file, index -> tuple(meta, file, index, []) }

        // Combine both streams
        ch_mosdepth_in = ch_mosdepth_targeted.mix(ch_mosdepth_other)

        MOSDEPTH (
            ch_mosdepth_in,
            ch_fasta
        )
        
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.map { _meta, file -> file }.collect())
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    }

    ch_preseq_in = ch_assay_split.wgs.mix(ch_assay_split.wes, ch_assay_split.rna, ch_assay_split.atac, ch_assay_split.chip)
    if (!params.skip_tools?.contains('preseq')) {
         PRESEQ_LCEXTRAP(
             ch_preseq_in
         )
         ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.map { _meta, file -> file }.collect())
    }

    ch_verifybamid_in = ch_assay_split.wgs.mix(ch_assay_split.wes)
    
    // Filter out cram files so only bam files proceed
    ch_verifybamid_bam_only = ch_verifybamid_in.filter { meta, file, index -> 
        file.name.endsWith('.bam') 
    }
    // TODO, VerifyBamID might be usefull especially for cancer analyiis, find a way to extract it from metadata     
    if (ch_refvcf){
        if (!params.skip_tools?.contains('verifybamid')) {
            VERIFYBAMID_VERIFYBAMID(
                ch_verifybamid_bam_only,
                ch_refvcf     
            )
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.selfsm.map { _meta, file -> file }.collect())
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.log.map { _meta, file -> file }.collect())
            ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID.out.versions)
        }
    }
    if (!params.skip_tools?.contains('rseqc')) {
         RSEQC_BAMSTAT(
             ch_assay_split.rna.map { meta, bam, bai -> tuple(meta, bam) }
         )
         ch_multiqc_files = ch_multiqc_files.mix(RSEQC_BAMSTAT.out.txt.map { _meta, file -> file }.collect())
         ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)
    }

    if (!params.skip_tools?.contains('ataqv')) {
        ATAQV_ATAQV(
            ch_assay_split.atac.map{ meta, bam, bai -> tuple(meta, bam, bai, []) },
            "human",
            [],[],[],[]
        )
        ch_multiqc_files = ch_multiqc_files.mix(ATAQV_ATAQV.out.json.map { _meta, file -> file }.collect())
    }

        if (!params.skip_tools?.contains('methyldackel')) {
        METHYLDACKEL_MBIAS(
            ch_assay_split.meth,
            ch_fasta,
            ch_fai,
        )
        ch_multiqc_files = ch_multiqc_files.mix(METHYLDACKEL_MBIAS.out.plots.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(METHYLDACKEL_MBIAS.out.versions)
    }

    if (!params.skip_tools?.contains('phantompeakqualtools')) {
        PHANTOMPEAKQUALTOOLS(
            ch_assay_split.chip.map { meta, bam, bai -> tuple(meta, bam) }
        )
        ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.spp.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions)
    }

    // TODO: find a way to extract gender, run this only if it is missing
    if (!params.skip_tools?.contains('ngsbits_samplegender')) {
         NGSBITS_SAMPLEGENDER(
             ch_assay_split.wgs,
             ch_fasta,
             ch_fai,
             params.samplegender_method ?: 'xy'
         )
         ch_multiqc_files = ch_multiqc_files.mix(NGSBITS_SAMPLEGENDER.out.tsv.map { _meta, file -> file }.collect())
         ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions)
    }

    if (!params.skip_tools?.contains('cramino')) {
         CRAMINO(
             ch_assay_split.longread
         )
         ch_multiqc_files = ch_multiqc_files.mix(CRAMINO.out.stats.map { _meta, file -> file }.collect())
         ch_multiqc_files = ch_multiqc_files.mix(CRAMINO.out.arrow.map { _meta, file -> file }.collect())
     }

    emit:
    ch_versions          = ch_versions
    ch_multiqc_files     = ch_multiqc_files 
}