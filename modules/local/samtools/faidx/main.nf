process SAMTOOLS_FAIDX {
    tag "$fasta"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("*.{fa,fasta}") , emit: fa, optional: true
    tuple val(meta), path ("*.sizes")      , emit: sizes, optional: true
    tuple val(meta), path ("*.fai")        , emit: fai, optional: true
    tuple val(meta), path ("*.gzi")        , emit: gzi, optional: true
    tuple val("${task.process}"), val('samtools'), eval('samtools --version | head -1 | awk "{print $2}"'),    topic: versions

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        faidx \\
        $fasta \\
        $args
    """
}
