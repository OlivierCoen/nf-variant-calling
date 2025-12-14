include { FASTP                                              } from '../../../modules/local/fastp'
include { FASTQC                                             } from '../../../modules/local/fastqc'


workflow READS_PREPARATION {

    take:
    ch_reads

    main:

    // ---------------------------------------------------------------------
    // Quality control on raw reads
    // ---------------------------------------------------------------------

    if ( !params.skip_fastqc ) {
        FASTQC ( ch_reads )
    }

    // ---------------------------------------------------------------------
    // Trimming / Filtering
    // ---------------------------------------------------------------------

    if ( !params.skip_fastp ) {

        FASTP (
            ch_reads,
            [], false, false, false
        )
        ch_reads      = FASTP.out.reads
    }

    emit:
    reads           = ch_reads
}
