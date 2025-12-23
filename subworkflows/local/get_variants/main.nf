include { CALL_VARIANTS                                 } from '../call_variants'
include { MERGE_VARIANTS                                } from '../merge_variants'
include { VARIANT_STATISTICS as STATISTICS            } from '../variant_statistics'



workflow GET_VARIANTS {

    take:
    ch_bam_bai
    ch_genome_fai_dict
    ch_genome_region_file
    variant_type

    main:

    // -----------------------------------------------------------------
    // CALL VARIANTS
    // -----------------------------------------------------------------

    CALL_VARIANTS(
        ch_bam_bai,
        ch_genome_fai_dict,
        ch_genome_region_file,
        variant_type
    )

    // -----------------------------------------------------------------
    // MERGE ALL DATA
    // -----------------------------------------------------------------

    MERGE_VARIANTS (
        CALL_VARIANTS.out.vcf,
        ch_genome_fai_dict.map { meta, fasta, fai, dict -> [ meta, fasta, fai ] }
    )
    ch_variants = MERGE_VARIANTS.out.variants

    // -----------------------------------------------------------------
    // STATS
    // -----------------------------------------------------------------

    STATISTICS (
        ch_variants,
        ch_genome_fai_dict
    )


    emit:
    variants            = ch_variants

}
