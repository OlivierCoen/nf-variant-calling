/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOME_PREPARATION                                    } from '../subworkflows/local/genome_preparation'
include { READS_PREPARATION                                     } from '../subworkflows/local/reads_preparation'
include { MAPPING_MARK_DUPLICATES                               } from '../subworkflows/local/mapping_mark_duplicates'
include { SNP_INDEL_SV_CALLING                                  } from '../subworkflows/local/snp_indel_sv_calling'
include { MERGE_CALLS as MERGE_SNPS_INDELS                      } from '../subworkflows/local/merge_calls'
include { MERGE_CALLS as MERGE_SVS                              } from '../subworkflows/local/merge_calls'
include { GENERATE_STATS                                        } from '../subworkflows/local/generate_stats'
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
    ch_genome_fai         = GENOME_PREPARATION.out.fai
    ch_genome_region_file = GENOME_PREPARATION.out.region_file
    ch_genome_dict        = GENOME_PREPARATION.out.dict

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
        ch_genome,
        ch_genome_fai,
        ch_genome_dict
    )

    // -----------------------------------------------------------------
    // VARIANT CALLING
    // -----------------------------------------------------------------

    SNP_INDEL_SV_CALLING(
        MAPPING_MARK_DUPLICATES.out.bam,
        MAPPING_MARK_DUPLICATES.out.bai,
        ch_genome,
        ch_genome_fai,
        ch_genome_region_file,
        ch_genome_dict
    )

    // -----------------------------------------------------------------
    // MERGE
    // -----------------------------------------------------------------


    MERGE_SNPS_INDELS (
        SNP_INDEL_SV_CALLING.out.snp_indel_per_region,
        ch_genome,
        ch_genome_fai
     )
/*
    MERGE_SVS (
        SNP_INDEL_SV_CALLING.out.sv_per_region,
        ch_genome,
        ch_genome_fai
     )
*/
    // -----------------------------------------------------------------
    // STATS
    // -----------------------------------------------------------------


    // -----------------------------------------------------------------
    // MULTIQC
    // -----------------------------------------------------------------

    ch_versions = ch_versions
                    .mix ( GENOME_PREPARATION.out.versions )
                    .mix ( MAPPING_MARK_DUPLICATES.out.versions )
                    //.mix ( MERGE_SNPS_INDELS.out.versions )


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
