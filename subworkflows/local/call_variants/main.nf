include { FREEBAYES                                             } from '../../../modules/local/freebayes'
include { DELLY_CALL                                            } from '../../../modules/local/delly/call'


workflow CALL_VARIANTS {

    take:
    ch_bam_bai
    ch_genome_fai_dict
    ch_genome_region_file
    variant_type

    main:

    ch_regions = ch_genome_region_file
                    .splitCsv(sep: "\t")
                    .map { row -> [chrom: row[0], start: row[1], end: row[2]] }

    ch_bam_bai_regions = ch_bam_bai.combine( ch_regions )

    if ( variant_type == "snps_indels" ) {

        // -----------------------------------------------------------------
        // SNPS AND INDELS
        // -----------------------------------------------------------------

        FREEBAYES (
            ch_bam_bai_regions,
            ch_genome_fai_dict.map{ meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect()
        )
        ch_vcf = FREEBAYES.out.vcf

    } else if ( variant_type == "svs" ) {

        // -----------------------------------------------------------------
        // SVs
        // -----------------------------------------------------------------

        DELLY_CALL(
            ch_bam_bai_regions,
            ch_genome_fai_dict.map{ meta, fasta, fai, dict -> [ meta, fasta, fai ] }.collect(),
            ch_genome_region_file.collect()
        )
        ch_vcf = DELLY_CALL.out.vcf

    } else {
        error("Unsupported variant type: $variant_type")
    }

    emit:
    vcf                 = ch_vcf

}
