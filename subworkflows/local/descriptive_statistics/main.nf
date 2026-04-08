include { BCFTOOLS_STATS                                            } from '../../../modules/local/bcftools/stats'



workflow DESCRIPTIVE_STATISTICS {

    take:
    ch_variants
    ch_genome_fai

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // VARIANT STATS
    // -----------------------------------------------------------------

    BCFTOOLS_STATS (
        ch_variants,
        ch_genome_fai.collect()
    )

}
