include { BCFTOOLS_SORT                                         } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_CONCAT                                       } from '../../../modules/nf-core/bcftools/concat'
include { BCFTOOLS_MERGE                                        } from '../../../modules/nf-core/bcftools/merge'
//include { BCFTOOLS_INDEX                                        } from '../modules/nf-core/bcftools/index'




workflow MERGE_CALLS {

    take:
    ch_variants_per_region
    ch_genome
    ch_genome_fai

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // ENSURE THAT THE INPUT FILES ARE SORTED BY CHR AND POSITION
    // -----------------------------------------------------------------

    BCFTOOLS_SORT ( ch_variants_per_region )

    // -----------------------------------------------------------------
    // CONCAT VERTICALLY ON ALL GENOME REGIONS, FOR EACH SAMPLE SEPARATELY
    // -----------------------------------------------------------------

    ch_to_concat = BCFTOOLS_SORT.out.vcf
                    .groupTuple()
                    .map { meta, vcfs -> [ meta, vcfs, [] ] }

    BCFTOOLS_CONCAT ( ch_to_concat )
/*
    // -----------------------------------------------------------------
    // MERGE HORIZONTALLY ON SAMPLES
    // -----------------------------------------------------------------

    ch_to_merge = BCFTOOLS_CONCAT.out.vcf
                    .map { meta, file -> file }
                    .collect()
                    .map { files -> [ files, [] ] }

    BCFTOOLS_MERGE (
        ch_to_merge,
        ch_genome.collect(),
        ch_genome_fai.collect()
        [:]
    )



    ch_versions = ch_versions
                    .mix ( BCFTOOLS_SORT.out.versions )
                    .mix ( BCFTOOLS_CONCAT.out.versions )
                    .mix ( BCFTOOLS_MERGE.out.versions )


    emit:
    vcf = BCFTOOLS_MERGE.out.vcf
    versions = ch_versions
*/
}
