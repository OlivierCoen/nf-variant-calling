include { BCFTOOLS_FILL_TAGS                             } from '../../../modules/local/bcftools/fill_tags'
include { BCFTOOLS_VIEW as BASE_FILTERING                } from '../../../modules/local/bcftools/view'
//include { BCFTOOLS_VIEW as ADVANCED_FILTERING            } from '../../../modules/local/bcftools/view'
include { ADDITIONAL_FILTERING                           } from '../../../modules/local/additional_filtering'
include { BCFTOOLS_INDEX                                 } from '../../../modules/local/bcftools/index'



workflow FILTER_VARIANTS {

    take:
    ch_variants
    min_depth_quantile
    max_depth_quantile


    main:

    ch_vcf = ch_variants.map { meta, vcf, tbi -> [ meta, vcf ] }
    ch_tbi = ch_variants.map { meta, vcf, tbi -> [ meta, tbi ] }

    // -----------------------------------------------------------------
    // IN CASE, ADD ALLMISSING TAGS TO VCF
    // -----------------------------------------------------------------

    BCFTOOLS_FILL_TAGS( ch_vcf )

    // -----------------------------------------------------------------
    // BASE FILTERING
    // -----------------------------------------------------------------

    BASE_FILTERING (
        BCFTOOLS_FILL_TAGS.out.bcf.join( ch_tbi )
    )

    // -----------------------------------------------------------------
    // ADVANCED FILTERING BASED ON ALLELE FREQUENCY AT THE SAMPLE LEVEL
    // -----------------------------------------------------------------

    //ADVANCED_FILTERING( BASE_FILTERING.out.vcf_tbi )

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
