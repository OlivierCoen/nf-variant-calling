process SAMTOOLS_MARKDUP {
    tag "${meta.id} - ${meta.lane}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(input_bam)

    output:
    tuple val(meta), path("*.bam"),     emit: bam
    path("*.markdup.log"),              topic: markdup_log
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | awk "{print $2}"'),    topic: versions


    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def args4 = task.ext.args4 ?: ''
    def prefix = "${input_bam.baseName}.markdup"
    """
    samtools collate ${args} -@ ${task.cpus-1} -O -u $input_bam \\
        | samtools fixmate ${args2} -@ ${task.cpus-1} -m -u - - \\
        | samtools sort ${args3} -@ ${task.cpus-1} -u - \\
        | samtools markdup ${args4} -@ ${task.cpus-1} -T $prefix -s - ${prefix}.bam \\
        2> ${prefix}.markdup.log
    """
}
