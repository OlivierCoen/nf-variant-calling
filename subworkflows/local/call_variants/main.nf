include { FREEBAYES                                             } from '../../../modules/local/freebayes'
include { DELLY_CALL                                            } from '../../../modules/local/delly/call'


workflow CALL_VARIANTS {

    take:
    ch_bam
    ch_bai
    ch_genome_fai_dict
    ch_genome_region_file
    skip_call_snps_indels
    skip_call_svs

    main:

    // -----------------------------------------------------------------
    // PREPARE INPUT FILES
    // -----------------------------------------------------------------

    ch_regions = ch_genome_region_file
                    .splitCsv(sep: "\t")
                    .map { row -> [chrom: row[0], start: row[1], end: row[2]] }


    ch_all_bams = ch_bam
                    .map { meta, bam -> bam }
                    .collect()
                    .map { bam_files -> [ bam_files ] }

    ch_all_bais = ch_bai
                    .map { meta, bai -> bai }
                    .collect()
                    .map { bai_files -> [ bai_files ] }

    ch_all_bam_bai_with_region = ch_all_bams
                                  .combine( ch_all_bais )
                                  .combine( ch_regions )

    ch_vcf = channel.empty()

    // -----------------------------------------------------------------
    // SNPS AND INDELS
    // -----------------------------------------------------------------

    if ( !skip_call_snps_indels ) {

        FREEBAYES (
            ch_all_bam_bai_with_region,
            ch_genome_fai_dict.map{ meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect()
        )

        ch_sns_indels = FREEBAYES.out.vcf
                           .map {
                                meta, vcf ->
                                    [ meta + [ variant_type: 'snp_indel' ], vcf ]
                            }
        ch_vcf = ch_vcf.mix( ch_sns_indels )

    }

    // -----------------------------------------------------------------
    // SVs
    // -----------------------------------------------------------------

    if ( !skip_call_svs ) {

        DELLY_CALL(
            ch_all_bam_bai_with_region,
            ch_genome_fai_dict.map{ meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect(),
            ch_genome_region_file.collect()
        )

        ch_svs = DELLY_CALL.out.vcf
                           .map {
                                meta, vcf ->
                                    [ meta + [ variant_type: 'structural_variant' ], vcf ]
                            }
        ch_vcf = ch_vcf.mix( ch_svs )

    }

    emit:
    vcf                 = ch_vcf

}
