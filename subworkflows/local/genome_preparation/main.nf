include { SAMTOOLS_FAIDX                                        } from '../../../modules/local/samtools/faidx'
include { MAKE_GENOME_REGIONS_WITH_OVERLAPS                     } from '../../../modules/local/make_genome_regions_with_overlaps'
include { GATK4_CREATESEQUENCEDICTIONARY                        } from '../../../modules/local/gatk4/createsequencedictionary'

workflow GENOME_PREPARATION {

    take:
    ch_genome
    reference_chunksize
    ratio_overlap_to_chunksize

    main:

    // -----------------------------------------------------------------
    // INDEX GENOME
    // -----------------------------------------------------------------

    SAMTOOLS_FAIDX( ch_genome )

    // -----------------------------------------------------------------
    // MAKE GENOME REGION CHUNKS
    // -----------------------------------------------------------------

    MAKE_GENOME_REGIONS_WITH_OVERLAPS(
        ch_genome,
        reference_chunksize,
        ratio_overlap_to_chunksize
    )

    // -----------------------------------------------------------------
    // MAKE GENOME REGION CHUNKS
    // -----------------------------------------------------------------

    GATK4_CREATESEQUENCEDICTIONARY ( ch_genome )



    emit:
    fai                 = SAMTOOLS_FAIDX.out.fai
    region_file         = MAKE_GENOME_REGIONS_WITH_OVERLAPS.out.regions
    dict                = GATK4_CREATESEQUENCEDICTIONARY.out.dict

}
