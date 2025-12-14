include { BWAMEM2_INDEX                                             } from '../../../modules/local/bwamem/index'
include { BWAMEM2_MEM                                               } from '../../../modules/nf-core/bwamem/mem'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/nf-core/samtools/flagstat'
include { SAMTOOLS_MARKDUP                                          } from '../../../modules/local/samtools/markdup'



workflow MAPPING_MARK_DUPLICATES {

    take:
    ch_reads
    ch_fasta

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // MAPPING
    // -----------------------------------------------------------------

    BWAMEM2_INDEX ( ch_fasta )

    BWAMEM2_MEM (
        ch_reads,
        BWAMEM2_INDEX.out.index,
        ch_fasta,
        true
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

    SAMTOOLS_FLAGSTAT(
        ch_bam.join ( SAMTOOLS_INDEX.out.bai )
    )



    ch_versions = ch_versions
                    .mix ( BWAMEM2_INDEX.out.versions )
                    .mix ( BWAMEM2_MEM.out.versions )
                    .mix ( SAMTOOLS_INDEX.out.versions )
                    .mix ( SAMTOOLS_FLAGSTAT.out.versions )


    emit:
    bam                     = ch_bam
    flagstat                = SAMTOOLS_FLAGSTAT.out.flagstat
    versions                = ch_versions
}
