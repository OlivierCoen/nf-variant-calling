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
include { VARIANT_ANALYSIS                                      } from '../subworkflows/local/variant_analysis'
include { REPORTING                                             } from '../subworkflows/local/reporting'

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

    GENOME_PREPARATION(
        ch_genome,
        params.reference_chunksize,
        params.ratio_overlap_to_chunksize
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

    /// -----------------------------------------------------------------
    // CALL VARIANTS
    // -----------------------------------------------------------------

    CALL_VARIANTS(
        MAPPING_MARK_DUPLICATES.out.bam,
        MAPPING_MARK_DUPLICATES.out.bai,
        ch_genome_fai_dict,
        ch_genome_region_file,
        params.skip_call_snps_indels,
        params.skip_call_svs
    )

    // -----------------------------------------------------------------
    // MERGE ALL DATA
    // -----------------------------------------------------------------

    MERGE_VARIANTS ( CALL_VARIANTS.out.vcf )

    ch_vcf_tbi = MERGE_VARIANTS.out.vcf_tbi

    // -----------------------------------------------------------------
    // FILTERING
    // -----------------------------------------------------------------

    FILTER_VARIANTS(
        ch_vcf_tbi,
        params.min_depth_quantile,
        params.max_depth_quantile
    )

    ch_filtered_vcf_tbi = FILTER_VARIANTS.out.vcf_tbi

    // -----------------------------------------------------------------
    // DESCRIPTIVE STATS
    // -----------------------------------------------------------------

    DESCRIPTIVE_STATISTICS (
        ch_filtered_vcf_tbi,
        ch_genome_fai_dict
    )

    // -----------------------------------------------------------------
    // STATISTICAL TESTS
    // -----------------------------------------------------------------

    VARIANT_ANALYSIS (
        ch_filtered_vcf_tbi.map{ meta, vcf, tbi -> [ meta, vcf ] },
        ch_design_file,
        params.window_size
    )

    ch_grouped_variants = VARIANT_ANALYSIS.out.grouped_variants

    // -----------------------------------------------------------------
    // REPORTING (DASH APP, MULTIQC, ...)
    // -----------------------------------------------------------------


    REPORTING(
        ch_filtered_vcf_tbi,
        ch_vcf_tbi,
        ch_grouped_variants,
        ch_genome_fai_dict,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
    )


    emit:
    multiqc_report            = REPORTING.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
