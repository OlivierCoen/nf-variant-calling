include { BCFTOOLS_SORT as BCFTOOLS_SORT_PER_CHUNK                       } from '../../../modules/local/bcftools/sort'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_ALL                             } from '../../../modules/local/bcftools/sort'
include { BCFTOOLS_CONCAT                                                } from '../../../modules/local/bcftools/concat'
include { BCFTOOLS_INDEX                                                 } from '../../../modules/local/bcftools/index'
include { SURVIVOR_MERGE                                                 } from '../../../modules/local/survivor/merge'


workflow MERGE_VARIANTS {

    take:
    ch_variants_per_region

    main:

    // -----------------------------------------------------------------
    // ENSURE THAT THE INPUT FILES ARE SORTED BY CHR AND POSITION
    // -----------------------------------------------------------------

    BCFTOOLS_SORT_PER_CHUNK ( ch_variants_per_region )

    // -----------------------------------------------------------------
    // CONCAT VERTICALLY ON ALL GENOME REGIONS
    // -----------------------------------------------------------------

    ch_chunk_bcf_sorted     = BCFTOOLS_SORT_PER_CHUNK.out.bcf.groupTuple()
    ch_chunk_indexes_sorted = BCFTOOLS_SORT_PER_CHUNK.out.csi.groupTuple()

    BCFTOOLS_CONCAT (
        ch_chunk_bcf_sorted.join( ch_chunk_indexes_sorted )
    )


    // -----------------------------------------------------------------
    // MERGE VARIANTS ACROSS CALLERS
    // -----------------------------------------------------------------
    /*
    ch_survivor_input = BCFTOOLS_CONCAT.out.bcf
                            .map{ meta, file -> [ [ id: meta.id ], file ] }
                            .groupTuple()

    SURVIVOR_MERGE( ch_survivor_input.collect() )
    */
    // -----------------------------------------------------------------
    // ENSURE THAT THE INPUT FILES ARE SORTED BY CHR AND POSITION
    // -----------------------------------------------------------------

    BCFTOOLS_SORT_ALL( BCFTOOLS_CONCAT.out.bcf )
    ch_concatenated_vcf = BCFTOOLS_SORT_ALL.out.vcf

    // -----------------------------------------------------------------
    // INDEX
    // -----------------------------------------------------------------

    BCFTOOLS_INDEX ( ch_concatenated_vcf )
    ch_vcf_tbi = ch_concatenated_vcf.join( BCFTOOLS_INDEX.out.tbi )

    emit:
    vcf_tbi             = ch_vcf_tbi

}
