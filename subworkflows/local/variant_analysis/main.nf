include { SEPARATE_VCF_DATA                             } from '../../../modules/local/separate_vcf_data'
include { STATISTICAL_TEST                              } from '../../../modules/local/statistical_test'
include { AGGREGATE_DATA                                } from '../../../modules/local/aggregate_data'



workflow VARIANT_ANALYSIS {

    take:
    ch_vcf
    ch_design_file
    statistical_test
    window_size

    main:

    // -----------------------------------------------------------------
    // EXTRACT VARIANT DESCRIPTORS (LIKE REF AND ALT COUNTS) AND SEPARATE THEM FROM BASIC VARIANT FEATURESs
    // -----------------------------------------------------------------

    SEPARATE_VCF_DATA ( ch_vcf )

    ch_variants   = SEPARATE_VCF_DATA.out.variants
    ch_ref_counts = SEPARATE_VCF_DATA.out.ref_counts
    ch_alt_counts = SEPARATE_VCF_DATA.out.alt_counts

    // -----------------------------------------------------------------
    // COMPUTE STATISTICAL TEST
    // -----------------------------------------------------------------

    STATISTICAL_TEST(
        ch_ref_counts.join( ch_alt_counts ),
        ch_design_file.collect(),
        statistical_test
    )

    // -----------------------------------------------------------------
    // MAKE SLIDING WINDOWS WHERE EACH WINDOW IS REPRESENTED BY THE 0.05 QUANTILE OF PVALUE
    // -----------------------------------------------------------------

    AGGREGATE_DATA(
        ch_vcf.join( STATISTICAL_TEST.out.pvalues ),
        window_size

    )

    emit:
    variants                   = AGGREGATE_DATA.out.variants
    grouped_variants           = AGGREGATE_DATA.out.grouped_variants

}
