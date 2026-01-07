include { BCFTOOLS_VIEW                                        } from '../../../modules/local/bcftools/view'



workflow FILTER_VARIANTS {

    take:
    ch_variants

    main:

    // -----------------------------------------------------------------
    // BASE FILTERING
    // -----------------------------------------------------------------

    BCFTOOLS_VIEW ( ch_variants )


    emit:
    variants = BCFTOOLS_VIEW.out.vcf_tbi

}
