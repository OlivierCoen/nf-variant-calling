include { DELLY_CALL                                            } from '../../../modules/local/delly/call'


workflow CALL_STRUCTURAL_VARIANTS {

    take:
    ch_bam
    ch_bai
    ch_genome
    ch_genome_fai
    ch_genome_region_file

    main:

    ch_versions = channel.empty()

    ch_regions = ch_genome_region_file
                    .splitCsv(sep: "\t")
                    .map { row -> [chrom: row[0], start: row[1], end: row[2]] }

    ch_bam_bai_regions = ch_bam
                            .join( ch_bai )
                            .combine( ch_regions )

    // -----------------------------------------------------------------
    // SVs
    // -----------------------------------------------------------------

    DELLY_CALL(
        ch_bam_bai_regions,
        ch_genome.join( ch_genome_fai ).collect(),
        ch_genome_region_file.collect()
    )

    emit:
    vcf                        = DELLY_CALL.out.vcf
    versions                   = ch_versions
}
