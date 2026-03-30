//
// STEP2: ANALYSIS TOOLS FOR BAM/CRAM FILES
//

include { SAMTOOLS_FAIDX                } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { MOSDEPTH                      } from '../../../modules/nf-core/mosdepth/main'
include { RSEQC_BAMSTAT                 } from '../../../modules/nf-core/rseqc/bamstat/main'
include { PRESEQ_LCEXTRAP               } from '../../../modules/nf-core/preseq/lcextrap/main'
include { NGSBITS_SAMPLEGENDER          } from '../../../modules/nf-core/ngsbits/samplegender/main'
include { SAMBAMBA_FLAGSTAT             } from '../../../modules/nf-core/sambamba/flagstat/main'
include { VERIFYBAMID_VERIFYBAMID       } from '../../../modules/nf-core/verifybamid/verifybamid/main'
include { ATAQV_ATAQV                   } from '../../../modules/nf-core/ataqv/ataqv/main'
include { PHANTOMPEAKQUALTOOLS          } from '../../../modules/nf-core/phantompeakqualtools/main'
include { CRAMINO                       } from '../../../modules/nf-core/cramino/main'
include { METHYLDACKEL_MBIAS            } from '../../../modules/nf-core/methyldackel/mbias/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_BAMS } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_STATS as SAMTOOLS_STATS_CRAMS} from '../../../modules/nf-core/samtools/stats/main'

workflow STEP2 {
    take:
    samplesheet   // channel: [val(meta), bam, bai]
    ch_fasta      // reference: [val(meta), .fasta]
    ch_fai        // reference: [val(meta), .fai]
    ch_intervals  // reference: [val(meta), .bed]
    ch_refvcf     // reference: [val(meta), .vcf]
    tools         // list of tools to run

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    def needs_fai = tools.contains('picard_collectmultiplemetrics') || 
                    tools.contains('methyldackel') || 
                    tools.contains('ngsbits_samplegender')

    if (!params.fasta_fai && needs_fai) {
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

    // SAMTOOLS STATS AND SAMTOOLS FLAGSTATS CAN RUN ALL ANALYSIS TYPES
    if (tools.contains('samtools_stats')) {
        SAMTOOLS_STATS(
            ch_final_bam_indexed,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    }

    if (tools.contains('sambamba_flagstat')) {
        // CRAM files crashing Sambamba
        ch_final_bam_indexed.branch { meta, file, index -> 
            bams:  file.name.endsWith('.bam')
            crams: file.name.endsWith('.cram')
        }.set { ch_flagstat_split }

        SAMBAMBA_FLAGSTAT(
            ch_flagstat_split.bams.map{meta, bam, _bai -> [meta, bam]}
        )
        ch_multiqc_files = ch_multiqc_files.mix(SAMBAMBA_FLAGSTAT.out.stats.map { _meta, file -> file }.collect())

        // only option for cram is to use samtools stats!
        SAMTOOLS_STATS_CRAMS(
            ch_flagstat_split.crams,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_CRAMS.out.stats.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_CRAMS.out.versions)
    }
    // SAMTOOLS STATS AND SAMBAMBA FLAGSTATS CAN RUN ALL ANALYSIS TYPES
    if (tools.contains('samtools_stats')) {
        SAMTOOLS_STATS_BAMS(
            ch_final_bam_indexed,
            ch_fasta
        )
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS_BAMS.out.stats.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(SAMTOOLS_STATS_BAMS.out.versions)
    }

    // Branch out by analysis type
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

    // Mosdepth: WGS needs no intervals while targeted methods need intervals
    if (tools.contains('mosdepth')) {
        
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
    if (tools.contains('preseq')) {
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
    if (ch_refvcf){
        if (tools.contains('verifybamid')) {
            VERIFYBAMID_VERIFYBAMID(
                status.tumor,
                ch_refvcf     
            )
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.selfsm.map { _meta, file -> file }.collect())
            ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID.out.log.map { _meta, file -> file }.collect())
            ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID.out.versions)
        }
    }
    if (tools.contains('rseqc')) {
         RSEQC_BAMSTAT(
             ch_assay_split.rna.map { meta, bam, bai -> tuple(meta, bam) }
         )
         ch_multiqc_files = ch_multiqc_files.mix(RSEQC_BAMSTAT.out.txt.map { _meta, file -> file }.collect())
         ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)
    }

    if (tools.contains('ataqv')) {
        ATAQV_ATAQV(
            ch_assay_split.atac.map{ meta, bam, bai -> tuple(meta, bam, bai, []) },
            "human",
            [],[],[],[]
        )
        ch_multiqc_files = ch_multiqc_files.mix(ATAQV_ATAQV.out.json.map { _meta, file -> file }.collect())
    }

    if (tools.contains('methyldackel')) {
        METHYLDACKEL_MBIAS(
            ch_assay_split.meth,
            ch_fasta,
            ch_fai,
        )
        ch_multiqc_files = ch_multiqc_files.mix(METHYLDACKEL_MBIAS.out.plots.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(METHYLDACKEL_MBIAS.out.versions)
    }

    if (tools.contains('phantompeakqualtools')) {
        PHANTOMPEAKQUALTOOLS(
            ch_assay_split.chip.map { meta, bam, bai -> tuple(meta, bam) }
        )
        ch_multiqc_files = ch_multiqc_files.mix(PHANTOMPEAKQUALTOOLS.out.spp.map { _meta, file -> file }.collect())
        ch_versions = ch_versions.mix(PHANTOMPEAKQUALTOOLS.out.versions)
    }

    if (tools.contains('ngsbits_samplegender')) {
        ch_ngsbits_in = ch_assay_split.wgs.filter { meta, bam, bai -> 
            def sex = meta.sex?.toString()?.trim()?.toUpperCase()
            return !sex || sex in ['NA', 'UNKNOWN', 'NULL', '']
        }
         NGSBITS_SAMPLEGENDER(
             ch_assay_split.wgs,
             ch_fasta,
             ch_fai,
             params.samplegender_method ?: 'xy'
         )
         ch_multiqc_files = ch_multiqc_files.mix(NGSBITS_SAMPLEGENDER.out.tsv.map { _meta, file -> file }.collect())
         ch_versions = ch_versions.mix(NGSBITS_SAMPLEGENDER.out.versions)
    }

    if (tools.contains('cramino')) {
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