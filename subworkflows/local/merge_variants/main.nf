include { BCFTOOLS_SORT as BCFTOOLS_SORT_CHUNK                           } from '../../../modules/local/bcftools/sort'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_ALL                             } from '../../../modules/local/bcftools/sort'
include { BCFTOOLS_CONCAT                                                } from '../../../modules/local/bcftools/concat'
include { BCFTOOLS_MERGE                                                 } from '../../../modules/local/bcftools/merge'




workflow MERGE_VARIANTS {

    take:
    ch_variants_per_region
    ch_genome_fai

    main:

    // -----------------------------------------------------------------
    // ENSURE THAT THE INPUT FILES ARE SORTED BY CHR AND POSITION
    // -----------------------------------------------------------------

    BCFTOOLS_SORT_CHUNK ( ch_variants_per_region )

    // -----------------------------------------------------------------
    // CONCAT VERTICALLY ON ALL GENOME REGIONS, FOR EACH SAMPLE SEPARATELY
    // -----------------------------------------------------------------

    ch_chunk_bcf_sorted = BCFTOOLS_SORT_CHUNK.out.bcf.groupTuple()
    ch_chunk_indexes_sorted = BCFTOOLS_SORT_CHUNK.out.csi.groupTuple()

    BCFTOOLS_CONCAT (
        ch_chunk_bcf_sorted.join( ch_chunk_indexes_sorted )
    )

    // -----------------------------------------------------------------
    // ENSURE THAT THE INPUT FILES ARE SORTED BY CHR AND POSITION
    // -----------------------------------------------------------------

    BCFTOOLS_SORT_ALL ( BCFTOOLS_CONCAT.out.bcf )

    // -----------------------------------------------------------------
    // MERGE HORIZONTALLY ON SAMPLES
    // -----------------------------------------------------------------

    BCFTOOLS_MERGE (
        BCFTOOLS_SORT_ALL.out.bcf.map { meta, file -> file }.collect(),
        BCFTOOLS_SORT_ALL.out.csi.map { meta, file -> file }.collect(),
        ch_genome_fai.collect()
    )

    ch_vcf_tbi = BCFTOOLS_MERGE.out.vcf.join( BCFTOOLS_MERGE.out.tbi )

    emit:
    variants             = ch_vcf_tbi

}
