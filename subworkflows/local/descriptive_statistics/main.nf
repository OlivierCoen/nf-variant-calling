//include { GATK4_VARIANTEVAL                                         } from '../../../modules/local/gatk4/varianteval'
include { BCFTOOLS_STATS                                            } from '../../../modules/local/bcftools/stats'



workflow DESCRIPTIVE_STATISTICS {

    take:
    ch_variants
    ch_genome_fai_dict

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // VARIANT STATS
    // -----------------------------------------------------------------

    BCFTOOLS_STATS (
        ch_variants,
        ch_genome_fai_dict.map { meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect()
    )

    /*
    GATK4_VARIANTEVAL(
        ch_variants,
        ch_genome_fai_dict.collect()
    )
    */

}
