include { BWAMEM2_INDEX                                             } from '../../../modules/local/bwamem2/index'
include { BWAMEM2_MEM                                               } from '../../../modules/local/bwamem2/mem'
include { SAMTOOLS_MARKDUP                                          } from '../../../modules/local/samtools/markdup'
include { SAMTOOLS_MERGE                                            } from '../../../modules/local/samtools/merge'
include { PICARD_ADDORREPLACEREADGROUPS                             } from '../../../modules/nf-core/picard/addorreplacereadgroups'
include { SAMTOOLS_INDEX                                            } from '../../../modules/local/samtools/index'
include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/local/samtools/flagstat'



workflow MAPPING_MARK_DUPLICATES {

    take:
    ch_reads
    ch_bam
    ch_genome_fai

    main:

    // -----------------------------------------------------------------
    // MAPPING
    // -----------------------------------------------------------------

    BWAMEM2_INDEX (
        ch_genome_fai.map { meta, genome, fai -> [ meta, genome ] }
    )

    BWAMEM2_MEM (
        ch_reads,
        BWAMEM2_INDEX.out.index.collect()
    )

    ch_all_bam = ch_bam.mix( BWAMEM2_MEM.out.bam )

    // -----------------------------------------------------------------
    // ADD READ GROUP IN HEADER
    // -----------------------------------------------------------------

    PICARD_ADDORREPLACEREADGROUPS (
        ch_all_bam,
        ch_genome_fai.collect()
    )

    // -----------------------------------------------------------------
    // MERGE BY SAMPLE
    // -----------------------------------------------------------------

    ch_grouped_bams = PICARD_ADDORREPLACEREADGROUPS.out.bam
                        .map { meta, file -> [ [id: meta.id], file ] } // removing lane info
                        .groupTuple()

    SAMTOOLS_MERGE ( ch_grouped_bams )
    ch_merged_bam = SAMTOOLS_MERGE.out.bam

    // -----------------------------------------------------------------
    // MARK DUPLICATES
    // -----------------------------------------------------------------

    SAMTOOLS_MARKDUP ( ch_merged_bam )
    ch_markdupped_bam = SAMTOOLS_MARKDUP.out.bam

    // -----------------------------------------------------------------
    // INDEXING BAM
    // -----------------------------------------------------------------

    SAMTOOLS_INDEX ( ch_markdupped_bam )
    ch_markdupped_bai = SAMTOOLS_INDEX.out.bai

    // -----------------------------------------------------------------
    // MAPPING STATS
    // -----------------------------------------------------------------

    SAMTOOLS_FLAGSTAT(
        ch_markdupped_bam.join( ch_markdupped_bai )
    )


    emit:
    bam                = ch_markdupped_bam
    bai                = ch_markdupped_bai

}
