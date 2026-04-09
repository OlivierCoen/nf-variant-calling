process BWAMEM2_MEM {
    tag "${meta.id} - ${meta.lane}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c6/c622ec4279b71ea9f3edd09124efdf118431111b3129e10fb1faad05ffdde7c7/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samblaster_samtools:0b5cd60ff7adfcb8' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('bwamem2'),    eval("bwa-mem2 version 2>&1 | tail -1 | sed 's/.* //'"),                  topic: versions
    tuple val("${task.process}"), val('samtools'),   eval("samtools --version | sed '1!d; s/samtools //'"),                    topic: versions


    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.lane}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    bwa-mem2 \\
        mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools view -b -@ $task.cpus $args2 - \\
        | samtools sort $args3 -@ $task.cpus -o ${prefix}.bam -

    """

}
