include { BWAMEM2_INDEX                                             } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM2_MEM                                               } from '../../../modules/local/bwamem2/mem'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/local/samtools/flagstat'
include { SAMTOOLS_MARKDUP                                          } from '../../../modules/local/samtools/markdup'



workflow MAPPING_MARK_DUPLICATES {

    take:
    ch_reads
    ch_genome

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // MAPPING
    // -------------------------- ---------------------------------------

    BWAMEM2_INDEX ( ch_genome )

    BWAMEM2_MEM (
        ch_reads,
        BWAMEM2_INDEX.out.index
    )

    // -----------------------------------------------------------------
    // MARK DUPLICATES
    // -----------------------------------------------------------------

    SAMTOOLS_MARKDUP ( BWAMEM2_MEM.out.bam )
    ch_bam = SAMTOOLS_MARKDUP.out.bam

    // -----------------------------------------------------------------
    // STATS
    // -----------------------------------------------------------------

    SAMTOOLS_INDEX ( ch_bam )
    ch_bai = SAMTOOLS_INDEX.out.bai

    SAMTOOLS_FLAGSTAT(
        ch_bam.join ( ch_bai )
    )


    ch_versions = ch_versions
                    .mix ( BWAMEM2_INDEX.out.versions )
                    .mix ( SAMTOOLS_INDEX.out.versions )


    emit:
    bam                     = ch_bam
    bai                     = ch_bai
    versions                = ch_versions
}
