process COCHRAN_MANTEL_HAENSEL_TEST {

    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4cb08d96e62942e7b6288abf2cfd30e813521a022459700e610325a3a7c0b1c8/data':
        'community.wave.seqera.io/library/bioconductor-geoquery_r-base_r-dplyr_r-optparse:fcd002470b7d6809' }"

    input:
    tuple val(meta), path(reference_count_file), path(alternative_count_file)
    path design

    output:
    //path("*.counts.csv"),                    optional: true,                                                    emit: counts
    tuple val("${task.process}"), val('R'),     eval('Rscript -e "cat(R.version.string)" | sed "s/R version //"'), topic: versions
    tuple val("${task.process}"), val('dplyr'), eval('Rscript -e "cat(as.character(packageVersion(\'dplyr\')))"'), topic: versions

    script:
    """
    compute_cmh_test.R \\
        --RO $reference_count_file \\
        --AO $alternative_count_file \\
        --design $design
    """

}
