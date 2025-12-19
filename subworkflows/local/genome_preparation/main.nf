include { SAMTOOLS_FAIDX                                        } from '../../../modules/local/samtools/faidx'
include { MAKE_GENOME_REGIONS                                   } from '../../../modules/local/make_genome_regions'
include { GATK4_CREATESEQUENCEDICTIONARY                        } from '../../../modules/nf-core/gatk4/createsequencedictionary'

workflow GENOME_PREPARATION {

    take:
    ch_genome
    reference_chunksize

    main:

    ch_versions = Channel.empty()

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



    ch_versions = ch_versions
                    .mix ( GATK4_CREATESEQUENCEDICTIONARY.out.versions )


    emit:
    fai                 = SAMTOOLS_FAIDX.out.fai
    region_file         = MAKE_GENOME_REGIONS.out.regions
    dict                = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    versions            = ch_versions
}
