include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/local/samtools/flagstat'
include { GATK4_VARIANTEVAL                                         } from '../../../modules/local/gatk4/varianteval'
include { BCFTOOLS_STATS                                            } from '../../../modules/nf-core/bcftools/stats'
include { BCFTOOLS_INDEX                                            } from '../../../modules/nf-core/bcftools/index'



workflow GENERATE_STATS {

    take:
    ch_sn_indel_vcf
    ch_sv_vcf
    ch_bam
    ch_bai
    ch_genome
    ch_genome_fai

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // MAPPING STATS
    // -----------------------------------------------------------------

    SAMTOOLS_FLAGSTAT(
        ch_bam.join ( ch_bai )
    )

    // -----------------------------------------------------------------
    // VARIANT STATS
    // -----------------------------------------------------------------

    ch_vcf_stat_input = ch_sn_indel_vcf.mix( ch_sv_vcf )

    BCFTOOLS_INDEX ( ch_vcf_stat_input )
    ch_vcf_index = BCFTOOLS_INDEX.out.tbi

    BCFTOOLS_STATS (
        ch_vcf_stat_input.join( ch_vcf_index ),
        [], [], [], [],
        ch_genome.collect()
    )

    GATK4_VARIANTEVAL(
        ch_vcf_stat_input
        ch_genome.join( ch_genome_fai ).collect()
    )

    ch_versions = ch_versions
                    .mix ( BCFTOOLS_STATS.out.versions )
                    .mix ( BCFTOOLS_INDEX.out.versions )


    emit:
    flagstat                = SAMTOOLS_FLAGSTAT.out.flagstat
    vcf_stats               = BCFTOOLS_STATS.out.stats
    variant_eval            = GATK4_VARIANTEVAL.out.output
    versions                = ch_versions
}
