include { EVALUATE_EFFECT_OF_FILTERS             } from '../../../modules/local/evaluate_effect_of_filters'
include { DASH_APP                               } from '../../../modules/local/dash_app'
include { MULTIQC                                } from '../../../modules/nf-core/multiqc'

include { methodsDescriptionText                 } from '../utils_nfcore_nf_variant_calling_pipeline'
include { paramsSummaryMultiqc                   } from '../../nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap                       } from 'plugin/nf-schema'

/*
========================================================================================
    SUBWORKFLOW TO DOWNLOAD EXPRESSIONATLAS ACCESSIONS AND DATASETS
========================================================================================
*/

workflow REPORTING {

    take:
    ch_filtered_vcf_tbi
    ch_vcf_tbi
    ch_grouped_variants
    ch_genome_fai_dict
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // -----------------------------------------------------------------
    // MAKE SCATTERPLOT OF FILTERED VARIANTS AGAINST VARIANTS
    // -----------------------------------------------------------------

    ch_input = ch_vcf_tbi
                 .join( ch_filtered_vcf_tbi )
                 .map{
                    meta, vcf, tbi, filtered_vcf, filtered_tbi ->
                        [ meta, vcf, filtered_vcf ]
                }

    EVALUATE_EFFECT_OF_FILTERS(
        ch_input,
        ch_genome_fai_dict.map { meta, genome, fai, dict -> [ meta, fai ] }.collect()
    )

    // -----------------------------------------------------------------
    // DASH APPLICATION
    // -----------------------------------------------------------------

    DASH_APP(
        ch_grouped_variants.map{ meta, file -> file }.collect()
    )

    // ------------------------------------------------------------------------------------
    // VERSIONS
    // ------------------------------------------------------------------------------------

    ch_versions = ch_versions
                    .mix ( DASH_APP.out.versions )

    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
                            .mix(topic_versions_string)
                            .collectFile(
                                storeDir: "${outdir}/pipeline_info",
                                name: 'nf_core_'  +  'variant_calling_software_'  + 'mqc_'  + 'versions.yml',
                                sort: true,
                                newLine: true
                            )

    // ------------------------------------------------------------------------------------
    // CONFIG
    // ------------------------------------------------------------------------------------

    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    ch_multiqc_custom_config = multiqc_config ?
        channel.fromPath(multiqc_config, checkIfExists: true) :
        channel.empty()

    ch_multiqc_logo          = multiqc_logo ?
        channel.fromPath(multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_files = ch_multiqc_files
        .mix( ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml') )

    ch_multiqc_custom_methods_description = multiqc_methods_description ?
        file(multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description     = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files
        .mix( ch_collated_versions )
        .mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

    ch_multiqc_files = ch_multiqc_files
                        .mix ( channel.topic("fastqc_zip") )
                        .mix ( channel.topic("fastp_json") )
                        .mix ( channel.topic("markdup_log") )
                        .mix ( channel.topic("flagstat") )
                        .mix ( channel.topic('bcftools_stats') )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    report = MULTIQC.out.report
}
