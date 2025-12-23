include { CALL_SNPS_INDELS as CALL                                     } from '../call_snps_indels'
include { MERGE_VARIANTS   as MERGE                                    } from '../merge_variants'



workflow SNPS_INDELS {

    take:
    ch_bam
    ch_bai
    ch_genome
    ch_genome_fai
    ch_genome_region_file
    ch_genome_dict

    main:

    ch_versions = channel.empty()

    CALL(
        ch_bam,
        ch_bai,
        ch_genome,
        ch_genome_fai,
        ch_genome_region_file,
        ch_genome_dict
    )

    MERGE (
        CALL.out.vcf,
        ch_genome,
        ch_genome_fai
    )


    ch_versions = ch_versions
                    .mix( CALL.out.versions )

    emit:
    vcf                 = MERGE.out.vcf
    index               = MERGE.out.index
    versions            = ch_versions
}
