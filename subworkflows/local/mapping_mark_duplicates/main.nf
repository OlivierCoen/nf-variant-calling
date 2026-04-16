include { SEQKIT_SPLIT2 as SPLIT_FASTQ_INTO_CHUNKS                  } from '../../../modules/local/seqkit/split2'
include { BWAMEM2_INDEX                                             } from '../../../modules/local/bwamem2/index'
include { BWAMEM2_MEM                                               } from '../../../modules/local/bwamem2/mem'
include { SAMTOOLS_MARKDUP                                          } from '../../../modules/local/samtools/markdup'
include { SAMTOOLS_MERGE as MERGE_BAMS_FROM_DIFFERENT_CHUNKS        } from '../../../modules/local/samtools/merge'
include { SAMTOOLS_MERGE as MERGE_BAMS_FROM_DIFFERENT_LANES         } from '../../../modules/local/samtools/merge'
include { PICARD_ADDORREPLACEREADGROUPS                             } from '../../../modules/nf-core/picard/addorreplacereadgroups'
include { SAMTOOLS_INDEX                                            } from '../../../modules/local/samtools/index'
include { SAMTOOLS_FLAGSTAT                                         } from '../../../modules/local/samtools/flagstat'



workflow MAPPING_MARK_DUPLICATES {

    take:
    ch_reads
    ch_bam
    ch_genome_fai
    nb_seqs_fastq_chunks

    main:

    // -----------------------------------------------------------------
    // SPLIT READS FOR BETTER SCALABILITY
    // -----------------------------------------------------------------

    SPLIT_FASTQ_INTO_CHUNKS(
        ch_reads,
        nb_seqs_fastq_chunks
    )

    // group paires of reads (or reads alone if single end) per chunk
    ch_reads_chunks = SPLIT_FASTQ_INTO_CHUNKS.out.reads
                        .transpose() // split list of reads
                        .map {
                            meta, file ->
                                def part = file.name.tokenize('.')[-3]
                                new_meta = meta + [ part: part ]
                                [ new_meta, file ]
                        }
                        .groupTuple()

    // -----------------------------------------------------------------
    // MAPPING
    // -----------------------------------------------------------------

    BWAMEM2_INDEX (
        ch_genome_fai.map { meta, genome, fai -> [ meta, genome ] }
    )

    BWAMEM2_MEM (
        ch_reads_chunks,
        BWAMEM2_INDEX.out.index.collect()
    )

    // -----------------------------------------------------------------
    // MERGE CHUNKS
    // -----------------------------------------------------------------
    ch_merge_chunk_bam_input = BWAMEM2_MEM.out.bam
                                .map{ meta, bam -> [ [id: meta.id, lane: meta.lane], bam ] }
                                .groupTuple()

    MERGE_BAMS_FROM_DIFFERENT_CHUNKS ( ch_merge_chunk_bam_input )

    // -----------------------------------------------------------------
    // ADD READ GROUP IN HEADER
    // -----------------------------------------------------------------

    ch_all_bam = ch_bam.mix( MERGE_BAMS_FROM_DIFFERENT_CHUNKS.out.bam )

    PICARD_ADDORREPLACEREADGROUPS (
        ch_all_bam,
        ch_genome_fai.collect()
    )

    // -----------------------------------------------------------------
    // MERGE BY SAMPLE
    // -----------------------------------------------------------------

    ch_grouped_bams = PICARD_ADDORREPLACEREADGROUPS.out.bam
                        .map { meta, file -> [ [id: meta.id], file ] } // removing lane info: mergine by sample
                        .groupTuple()

    MERGE_BAMS_FROM_DIFFERENT_LANES ( ch_grouped_bams )
    ch_merged_bam = MERGE_BAMS_FROM_DIFFERENT_LANES.out.bam

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
