process SAMTOOLS_FLAGSTAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path("*.flagstat"), topic: flagstat
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | awk "{print $2}"'),    topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools \\
        flagstat \\
        --threads ${task.cpus} \\
        $bam \\
        > ${prefix}.flagstat
    """

}
