include { BWAMEM2_INDEX                                             } from '../../../modules/local/bwamem2/index'
include { BWAMEM2_MEM                                               } from '../../../modules/local/bwamem2/mem'
include { SAMTOOLS_MARKDUP                                          } from '../../../modules/local/samtools/markdup'
include { GATK4_MARKDUPLICATES                                      } from '../../../modules/local/gatk4/markduplicates'
include { SAMTOOLS_MERGE                                            } from '../../../modules/local/samtools/merge'
include { GATK4_ADDORREPLACEREADGROUPS                              } from '../../../modules/local/gatk4/addorreplacereadgroups'
include { SAMTOOLS_INDEX                                            } from '../../../modules/local/samtools/index'
include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/local/samtools/flagstat'



workflow MAPPING_MARK_DUPLICATES {

    take:
    ch_reads
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

    // -----------------------------------------------------------------
    // ADD READ GROUP IN HEADER
    // -----------------------------------------------------------------

    GATK4_ADDORREPLACEREADGROUPS (
        BWAMEM2_MEM.out.bam,
        ch_genome_fai_dict.map { meta, genome, fai, dict -> [ meta, genome, fai ] }.collect()
    )

    // -----------------------------------------------------------------
    // MERGE BY SAMPLE
    // -----------------------------------------------------------------

    ch_grouped_bams = GATK4_ADDORREPLACEREADGROUPS.out.bam
                        .map { meta, file -> [ [id: meta.id], file ] } // removing lane info
                        .groupTuple()

    SAMTOOLS_MERGE ( ch_grouped_bams )
    ch_bam = SAMTOOLS_MERGE.out.bam

    // -----------------------------------------------------------------
    // MARK DUPLICATES
    // -----------------------------------------------------------------

    if ( params.markdup_method == "gatk" ) {

        GATK4_MARKDUPLICATES (
            ch_bam,
            ch_genome_fai_dict.collect()
        )
        ch_markdupped_bam = GATK4_MARKDUPLICATES.out.output

    } else {

        SAMTOOLS_MARKDUP ( ch_bam )
        ch_markdupped_bam = SAMTOOLS_MARKDUP.out.bam

    }

    // -----------------------------------------------------------------
    // INDEXING BAM
    // -----------------------------------------------------------------

    SAMTOOLS_INDEX ( ch_markdupped_bam )
    ch_bam_bai = ch_markdupped_bam.join( SAMTOOLS_INDEX.out.bai )

    // -----------------------------------------------------------------
    // MAPPING STATS
    // -----------------------------------------------------------------

    SAMTOOLS_FLAGSTAT( ch_bam_bai )


    emit:
    bam_bai                 = ch_bam_bai

}
