include { FASTP                                              } from '../../../modules/nf-core/fastp'
include { FASTQC                                             } from '../../../modules/local/fastqc'


workflow READS_PREPARATION {

    take:
    ch_reads

    main:

    ch_versions = Channel.empty()

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
            [], false, false, true
        )
        ch_reads = FASTP.out.reads
        ch_versions = ch_versions.mix ( FASTP.out.versions )
    }

    emit:
    reads           = ch_reads
    versions        = ch_versions
}
