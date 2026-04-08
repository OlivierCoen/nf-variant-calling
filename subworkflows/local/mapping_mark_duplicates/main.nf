include { BWAMEM2_INDEX                                             } from '../../../modules/local/bwamem2/index'
include { BWAMEM2_MEM_MARKDUP                                       } from '../../../modules/local/bwamem2/mem_markdup'
include { SAMTOOLS_MERGE                                            } from '../../../modules/local/samtools/merge'
include { PICARD_ADDORREPLACEREADGROUPS                             } from '../../../modules/nf-core/picard/addorreplacereadgroups'
include { SAMTOOLS_INDEX                                            } from '../../../modules/local/samtools/index'
include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/local/samtools/flagstat'



workflow MAPPING_MARK_DUPLICATES {

    take:
    ch_reads
    ch_bam
    ch_genome_fai_dict

    main:

    // -----------------------------------------------------------------
    // MAPPING
    // -----------------------------------------------------------------

    BWAMEM2_INDEX (
        ch_genome_fai_dict.map { meta, genome, fai, dict -> [ meta, genome ] }
    )

    BWAMEM2_MEM (
        ch_reads,
        BWAMEM2_INDEX.out.index.collect()
    )

    ch_all_bam = ch_bam.mix( BWAMEM2_MEM.out.bam )

    // -----------------------------------------------------------------
    // ADD READ GROUP IN HEADER
    // -----------------------------------------------------------------

    GATK4_ADDORREPLACEREADGROUPS (
        ch_all_bam,
        ch_genome_fai_dict.map { meta, genome, fai, dict -> [ meta, genome, fai ] }.collect()
    )

    // -----------------------------------------------------------------
    // MERGE BY SAMPLE
    // -----------------------------------------------------------------

    ch_grouped_bams = PICARD_ADDORREPLACEREADGROUPS.out.bam
                        .map { meta, file -> [ [id: meta.id], file ] } // removing lane info
                        .groupTuple()

    SAMTOOLS_MERGE ( ch_grouped_bams )
    ch_bam = SAMTOOLS_MERGE.out.bam

    // -----------------------------------------------------------------
    // INDEXING BAM
    // -----------------------------------------------------------------

    SAMTOOLS_INDEX ( ch_bam )
    ch_markdupped_bai = SAMTOOLS_INDEX.out.bai

    // -----------------------------------------------------------------
    // MAPPING STATS
    // -----------------------------------------------------------------

    SAMTOOLS_FLAGSTAT(
        ch_bam.join( ch_markdupped_bai )
    )


    emit:
    bam                = ch_bam
    bai                = ch_markdupped_bai

}
