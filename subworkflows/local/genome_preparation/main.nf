include { SAMTOOLS_FAIDX                                        } from '../../../modules/local/samtools/faidx'
include { MAKE_GENOME_REGIONS                                   } from '../../../modules/local/make_genome_regions'
include { GATK4_CREATESEQUENCEDICTIONARY                        } from '../../../modules/local/gatk4/createsequencedictionary'

workflow GENOME_PREPARATION {

    take:
    ch_genome
    reference_chunksize

    main:

    // -----------------------------------------------------------------
    // INDEX GENOME
    // -----------------------------------------------------------------

    SAMTOOLS_FAIDX( ch_genome )

    // -----------------------------------------------------------------
    // MAKE GENOME REGION CHUNKS
    // -----------------------------------------------------------------

    MAKE_GENOME_REGIONS(
        ch_genome,
        reference_chunksize
    )

    // -----------------------------------------------------------------
    // MAKE GENOME REGION CHUNKS
    // -----------------------------------------------------------------

    GATK4_CREATESEQUENCEDICTIONARY ( ch_genome )



    emit:
    fai                 = SAMTOOLS_FAIDX.out.fai
    region_file         = MAKE_GENOME_REGIONS.out.regions
    dict                = GATK4_CREATESEQUENCEDICTIONARY.out.dict

}
