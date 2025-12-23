/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_PREPARATION                                    } from '../subworkflows/local/genome_preparation'
include { READS_PREPARATION                                     } from '../subworkflows/local/reads_preparation'
include { MAPPING_MARK_DUPLICATES                               } from '../subworkflows/local/mapping_mark_duplicates'
include { GET_VARIANTS as SNPS_INDELS                           } from '../subworkflows/local/get_variants'
include { GET_VARIANTS as STRUCTURAL_VARIANTS                   } from '../subworkflows/local/get_variants'

include { MULTIQC_WORKFLOW                                      } from '../subworkflows/local/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_CALLING {

    take:
    ch_reads // channel: samplesheet read in from --input
    ch_genome

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    GENOME_PREPARATION(
        ch_genome,
        params.reference_chunksize
    )

    ch_genome_fai_dict = ch_genome
                            .join( GENOME_PREPARATION.out.fai )
                            .join( GENOME_PREPARATION.out.dict )

    ch_genome_region_file = GENOME_PREPARATION.out.region_file

    // -----------------------------------------------------------------
    // PREPARE READS
    // -----------------------------------------------------------------

    READS_PREPARATION( ch_reads )
    ch_reads = READS_PREPARATION.out.reads

    // -----------------------------------------------------------------
    // MAP READS TO REFERENCE GENOME
    // -----------------------------------------------------------------

    MAPPING_MARK_DUPLICATES(
        ch_reads,
        ch_genome_fai_dict
    )
    ch_bam_bai = MAPPING_MARK_DUPLICATES.out.bam_bai

    // -----------------------------------------------------------------
    // VARIANT CALLING / MERGING / STATS
    // -----------------------------------------------------------------

    SNPS_INDELS(
        ch_bam_bai,
        ch_genome_fai_dict,
        ch_genome_region_file,
        "snps_indels"
    )

    STRUCTURAL_VARIANTS(
        ch_bam_bai,
        ch_genome_fai_dict,
        ch_genome_region_file,
        "svs"
    )

    // -----------------------------------------------------------------
    // MULTIQC
    // -----------------------------------------------------------------

    ch_versions = ch_versions
                    .mix ( GENOME_PREPARATION.out.versions )
                    .mix ( MAPPING_MARK_DUPLICATES.out.versions )

    MULTIQC_WORKFLOW(
        ch_multiqc_files,
        ch_versions,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
    )


    emit:
    multiqc_report            = MULTIQC_WORKFLOW.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
