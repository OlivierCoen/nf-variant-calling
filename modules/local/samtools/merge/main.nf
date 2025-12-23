process SAMTOOLS_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(bam_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${prefix}.bam") , optional:true, emit: bam
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | awk "{print $2}"'),    topic: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        merge \\
        --threads ${task.cpus-1} \\
        $args \\
        ${prefix}.bam \\
        $bam_files
    """
}
