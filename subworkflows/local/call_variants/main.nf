include { FREEBAYES                                             } from '../../../modules/local/freebayes'
include { DELLY_CALL                                            } from '../../../modules/local/delly/call'
include { MANTA_GERMLINE                                        } from '../../../modules/local/manta/germline'


workflow CALL_VARIANTS {

    take:
    ch_bam
    ch_bai
    ch_genome_fai_dict
    ch_genome_region_file
    callers

    main:

    // -----------------------------------------------------------------
    // PREPARE INPUT FILES
    // -----------------------------------------------------------------

    ch_regions = ch_genome_region_file
                    .splitCsv(sep: "\t")
                    .map { row -> [chrom: row[0], start: row[1], end: row[2]] }


    ch_all_bams = ch_bam
                    .map { meta, bam -> bam }
                    .toSortedList()
                    .map { bam_files -> [ bam_files ] }

    ch_all_bais = ch_bai
                    .map { meta, bai -> bai }
                    .toSortedList()
                    .map { bai_files -> [ bai_files ] }

    ch_all_bam_bai_with_region = ch_all_bams
                                  .combine( ch_all_bais )
                                  .combine( ch_regions )

    ch_vcf = channel.empty()

    caller_list = callers.tokenize(',')

    // -----------------------------------------------------------------
    // FREEBAYES
    // -----------------------------------------------------------------

    if ( "freebayes" in caller_list ) {

        FREEBAYES (
            ch_all_bam_bai_with_region,
            ch_genome_fai_dict.map{ meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect()
        )
        ch_freebayes_vcf = FREEBAYES.out.vcf.map{ meta, file -> [ meta + [ caller: "freebayes", type: "snp_indel"], file ] }
        ch_vcf = ch_vcf.mix( ch_freebayes_vcf )

    }

    // -----------------------------------------------------------------
    // DELLY
    // -----------------------------------------------------------------

    if ( "delly" in caller_list ) {

        DELLY_CALL(
            ch_all_bam_bai_with_region,
            ch_genome_fai_dict.map{ meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect(),
            ch_genome_region_file.collect()
        )
        ch_delly_vcf = DELLY_CALL.out.vcf
                        .map{ meta, file -> [ meta + [ caller: "delly", type: "sv"], file ] }
        ch_vcf = ch_vcf.mix( ch_delly_vcf )

    }

    // -----------------------------------------------------------------
    // MANTA
    // -----------------------------------------------------------------

    if ( "manta" in caller_list ) {

        MANTA_GERMLINE(
            ch_all_bam_bai_with_region,
            ch_genome_fai_dict.map{ meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect(),
            []
        )
        //ch_manta_small_indel_vcf = MANTA_GERMLINE.out.candidate_small_indels_vcf.map{ meta, file -> [ meta + [ caller: "manta", type: "sv"], file ] }
        ch_manta_sv_vcf          = MANTA_GERMLINE.out.candidate_sv_vcf.map{ meta, file -> [ meta + [ caller: "manta", type: "sv"], file ] }
        //ch_manta_diploid_sv_vcf  = MANTA_GERMLINE.out.diploid_sv_vcf.map{ meta, file -> [ meta + [ caller: "manta", type: "sv"], file ] }

        ch_vcf = ch_vcf
                    //.mix( ch_manta_small_indel_vcf )
                    .mix( ch_manta_sv_vcf )
                    //.mix( ch_manta_diploid_sv_vcf )

    }


    emit:
    vcf                 = ch_vcf

}
