include { GET_ALLELE_COUNTS                                      } from '../../../modules/local/get_allele_counts'
include { COCHRAN_MANTEL_HAENSEL_TEST                            } from '../../../modules/local/cochran_mantel_haensel_test'



workflow STATISTICAL_TESTS {

    take:
    ch_variants
    ch_design_file

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // VARIANT STATS
    // -----------------------------------------------------------------

    GET_ALLELE_COUNTS (
        ch_variants.map{ meta, vcf, tbi -> [ meta, vcf ] }
    )

    COCHRAN_MANTEL_HAENSEL_TEST(
        GET_ALLELE_COUNTS.out.counts,
        ch_design_file
    )

}
