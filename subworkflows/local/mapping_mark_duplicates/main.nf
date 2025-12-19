include { BWAMEM2_INDEX                                             } from '../../../modules/nf-core/bwamem2/index'
include { BWAMEM2_MEM                                               } from '../../../modules/local/bwamem2/mem'
include { SAMTOOLS_MARKDUP                                          } from '../../../modules/local/samtools/markdup'
include { GATK4SPARK_MARKDUPLICATES                                 } from '../../../modules/nf-core/gatk4spark/markduplicates'
include { SAMTOOLS_MERGE                                            } from '../../../modules/local/samtools/merge'
include { SAMTOOLS_INDEX                                            } from '../../../modules/nf-core/samtools/index'
include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/local/samtools/flagstat'



workflow MAPPING_MARK_DUPLICATES {

    take:
    ch_reads
    ch_genome
    ch_genome_fai
    ch_genome_dict

    main:

    ch_versions = Channel.empty()

    // -----------------------------------------------------------------
    // MAPPING
    // -----------------------------------------------------------------

    BWAMEM2_INDEX ( ch_genome )

    BWAMEM2_MEM (
        ch_reads,
        BWAMEM2_INDEX.out.index.collect()
    )

    // -----------------------------------------------------------------
    // MARK DUPLICATES
    // -----------------------------------------------------------------

    if ( params.markdup_method == "gatk" ) {

        GATK4SPARK_MARKDUPLICATES (
            BWAMEM2_MEM.out.bam,
            ch_genome.map { meta, file -> file },
            ch_genome_fai.map { meta, file -> file },
            ch_genome_dict.map { meta, file -> file },
        )
        ch_markdupped_bam = GATK4SPARK_MARKDUPLICATES.out.output
        ch_versions = ch_versions.mix(GATK4SPARK_MARKDUPLICATES.out.versions)

    } else {

        SAMTOOLS_MARKDUP ( BWAMEM2_MEM.out.bam )
        ch_markdupped_bam = SAMTOOLS_MARKDUP.out.bam

    }

    // -----------------------------------------------------------------
    // MERGE BY SAMPLE
    // -----------------------------------------------------------------

    ch_grouped_bams = ch_markdupped_bam
                        .map { meta, file -> [ [id: meta.id], file ] } // removing lane info
                        .groupTuple()

    SAMTOOLS_MERGE ( ch_grouped_bams )
    ch_bam = SAMTOOLS_MERGE.out.bam

    SAMTOOLS_INDEX ( ch_bam )



    ch_versions = ch_versions
                    .mix ( BWAMEM2_INDEX.out.versions )
                    .mix ( SAMTOOLS_INDEX.out.versions )


    emit:
    bam                     = ch_bam
    bai                     = SAMTOOLS_INDEX.out.bai
    versions                = ch_versions
}
