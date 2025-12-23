include { FREEBAYES                                             } from '../../../modules/local/freebayes'



workflow CALL_SNPS_INDELS {

    take:
    ch_bam
    ch_bai
    ch_genome
    ch_genome_fai
    ch_genome_region_file
    ch_genome_dict

    main:

    ch_versions = channel.empty()

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

    emit:
    vcf                 = FREEBAYES.out.vcf
    versions            = ch_versions

}
