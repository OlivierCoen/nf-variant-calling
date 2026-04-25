process STATISTICAL_TEST {

    tag "${meta.id} - ${meta.type}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        '':
        '' }"

    input:
    tuple val(meta), path(reference_count_file), path(alternative_count_file)
    path design
    val statistical_test

    output:
    tuple val(meta), path("*.cmh_pvalues.txt"),                                                                    emit: pvalues
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'), topic: versions
    tuple val("${task.process}"), val('dplyr'), eval('Rscript -e "cat(as.character(packageVersion(\'dplyr\')))"'), topic: versions

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    compute_statistical_test.R \\
        --method $statistical_test \\
        --RO $reference_count_file \\
        --AO $alternative_count_file \\
        --design $design \\
        --out ${prefix}.cmh_pvalues.txt
    """

}
