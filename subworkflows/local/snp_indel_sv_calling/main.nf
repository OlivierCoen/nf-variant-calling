include { FREEBAYES                                             } from '../../../modules/local/freebayes'
include { DELLY_CALL                                            } from '../../../modules/local/delly/call'
include { SAMTOOLS_FAIDX                                        } from '../../../modules/local/samtools/faidx'



workflow SNP_INDEL_SV_CALLING {

    take:
    ch_bam
    ch_bai
    ch_genome

    main:

    // -----------------------------------------------------------------
    // SNPS AND INDELS
    // -------------------------- ---------------------------------------

    FREEBAYES (
        ch_bam.join( ch_bai ),
        ch_genome
    )


    // -----------------------------------------------------------------
    // SVs
    // -------------------------- ---------------------------------------

    SAMTOOLS_FAIDX( ch_genome )

    DELLY_CALL(
        ch_bam.join( ch_bai ),
        ch_genome.join( SAMTOOLS_FAIDX.out.fai )
    )

    ch_vcf = FREEBAYES.out.vcf
                .mix( DELLY_CALL.out.vcf )


    emit:
    vcf                     = ch_vcf
}
