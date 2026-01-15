include { GET_ALLELE_COUNTS                   } from '../../../modules/local/get_allele_counts'
include { CMH_TEST                            } from '../../../modules/local/cmh_test'



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

    CMH_TEST(
        GET_ALLELE_COUNTS.out.counts,
        ch_design_file
    )

    emit:
    pvalues                             = CMH_TEST.out.pvalues

}
