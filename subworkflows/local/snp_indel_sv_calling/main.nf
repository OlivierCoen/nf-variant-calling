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

    ch_bam_bai_regions = ch_bam
                            .join( ch_bai )
                            .combine( ch_regions )

    // -----------------------------------------------------------------
    // SNPS AND INDELS
    // -----------------------------------------------------------------

    FREEBAYES (
        ch_bam_bai_regions,
        ch_genome.join( ch_genome_fai ).collect()
    )

    // -----------------------------------------------------------------
    // SVs
    // -----------------------------------------------------------------

    DELLY_CALL(
        ch_bam_bai_regions,
        ch_genome.join( ch_genome_fai ).collect(),
        ch_genome_region_file.collect()
    )

    emit:
    snp_indel_per_region                 = FREEBAYES.out.vcf
    sv_per_region                        = DELLY_CALL.out.vcf
}
