include { CALL_STRUCTURAL_VARIANTS as CALL                                     } from '../call_structural_variants'
include { MERGE_VARIANTS           as MERGE                                    } from '../merge_variants'



workflow STRUCTURAL_VARIANTS {

    take:
    ch_bam
    ch_bai
    ch_genome
    ch_genome_fai
    ch_genome_region_file

    main:

    ch_versions = channel.empty()

    CALL(
        ch_bam,
        ch_bai,
        ch_genome,
        ch_genome_fai,
        ch_genome_region_file
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
