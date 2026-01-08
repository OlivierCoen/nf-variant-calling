/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_PREPARATION                                    } from '../subworkflows/local/genome_preparation'
include { READS_PREPARATION                                     } from '../subworkflows/local/reads_preparation'
include { MAPPING_MARK_DUPLICATES                               } from '../subworkflows/local/mapping_mark_duplicates'
include { CALL_VARIANTS                                         } from '../subworkflows/local/call_variants'
include { MERGE_VARIANTS                                        } from '../subworkflows/local/merge_variants'
include { FILTER_VARIANTS                                       } from '../subworkflows/local/filter_variants'
include { DESCRIPTIVE_STATISTICS                                } from '../subworkflows/local/descriptive_statistics'
include { STATISTICAL_TESTS                                     } from '../subworkflows/local/statistical_tests'

include { MULTIQC_WORKFLOW                                      } from '../subworkflows/local/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_CALLING {

    take:
    ch_reads
    ch_design_file
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

    /// -----------------------------------------------------------------
    // CALL VARIANTS
    // -----------------------------------------------------------------

    CALL_VARIANTS(
        ch_bam_bai,
        ch_genome_fai_dict,
        ch_genome_region_file,
        params.skip_call_snps_indels,
        params.skip_call_svs
    )

    // -----------------------------------------------------------------
    // MERGE ALL DATA
    // -----------------------------------------------------------------

    MERGE_VARIANTS ( CALL_VARIANTS.out.vcf )

    // -----------------------------------------------------------------
    // FILTERING
    // -----------------------------------------------------------------

    FILTER_VARIANTS(
        MERGE_VARIANTS.out.vcf_tbi,
        ch_genome_fai_dict
    )
    ch_filtered_variants = FILTER_VARIANTS.out.vcf_tbi

    // -----------------------------------------------------------------
    // DESCRIPTIVE STATS
    // -----------------------------------------------------------------

    DESCRIPTIVE_STATISTICS (
        ch_filtered_variants,
        ch_genome_fai_dict
    )

    // -----------------------------------------------------------------
    // STATISTICAL TESTS
    // -----------------------------------------------------------------

    STATISTICAL_TESTS (
        ch_filtered_variants,
        ch_design_file
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
