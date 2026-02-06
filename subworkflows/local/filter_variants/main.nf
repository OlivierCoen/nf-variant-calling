include { BCFTOOLS_VIEW as BASE_FILTERING                } from '../../../modules/local/bcftools/view'
include { ADDITIONAL_FILTERING                           } from '../../../modules/local/additional_filtering'
include { BCFTOOLS_INDEX                                 } from '../../../modules/local/bcftools/index'



workflow FILTER_VARIANTS {

    take:
    ch_variants
    min_depth_quantile
    max_depth_quantile


    main:

    // -----------------------------------------------------------------
    // BASE FILTERING
    // -----------------------------------------------------------------

    BASE_FILTERING ( ch_variants )

    // -----------------------------------------------------------------
    // ADDITIONAL FILTERING
    // -----------------------------------------------------------------

    ADDITIONAL_FILTERING (
        BASE_FILTERING.out.vcf_tbi.map { meta, vcf, tbi -> [ meta, vcf ] },
        min_depth_quantile,
        max_depth_quantile
    )
    ch_filtered_vcf = ADDITIONAL_FILTERING.out.vcf

    // -----------------------------------------------------------------
    // INDEXING
    // -----------------------------------------------------------------

    BCFTOOLS_INDEX ( ch_filtered_vcf )

    emit:
    vcf_tbi = ch_filtered_vcf.join( BCFTOOLS_INDEX.out.tbi )

}
