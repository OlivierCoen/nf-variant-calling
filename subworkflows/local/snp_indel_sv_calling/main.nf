include { FREEBAYES                                             } from '../../../modules/local/freebayes'
include { DELLY_CALL                                            } from '../../../modules/local/delly/call'



workflow SNP_INDEL_SV_CALLING {

    take:
    ch_bam
    ch_bai
    ch_genome
    ch_genome_fai
    ch_genome_region_file
    ch_genome_dict

    main:

    ch_regions = ch_genome_region_file
                    .splitCsv(sep: "\t")
                    .map { row -> [chrom: row[0], start: row[1], end: row[2]] }


    // -----------------------------------------------------------------
    // SNPS AND INDELS
    // -----------------------------------------------------------------

    FREEBAYES (
        ch_bam.join( ch_bai ),
        ch_genome.join( ch_genome_fai ).collect(),
        ch_regions
    )
    ch_snp_indel_vcf = FREEBAYES.out.vcf


    // -----------------------------------------------------------------
    // SVs
    // -----------------------------------------------------------------

    DELLY_CALL(
        ch_bam.join( ch_bai ),
        ch_genome.join( ch_genome_fai ).collect(),
        ch_genome_region_file.combine( ch_regions )
    )
    ch_sv_vcf = DELLY_CALL.out.vcf


    emit:
    snp_indel_per_region                 = ch_snp_indel_vcf
    sv_per_region                        = ch_sv_vcf
}
