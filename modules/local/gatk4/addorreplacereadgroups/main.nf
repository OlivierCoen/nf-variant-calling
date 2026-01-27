process GATK4_ADDORREPLACEREADGROUPS {
    tag "${meta.id} - ${meta.lane}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.bam"),  emit: bam,  optional: true
    tuple val(meta), path("*.bai"),  emit: bai,  optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version 2>&1 | grep 'The Genome Analysis Toolkit' | cut -d' ' -f6"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.markdup"
    def suffix = task.ext.suffix ?: "${bam.getExtension()}"
    def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
    def create_index = suffix == "bam" ? "--CREATE_INDEX" : ""
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK AddOrReplaceReadGroups] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    if ("${bam}" == "${prefix}.${suffix}") {
        error("Input and output names are the same, use \"task.ext.prefix\" to disambiguate!")
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        AddOrReplaceReadGroups \\
        ${args} \\
        ${reference} \\
        ${create_index} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.${suffix} \\
        --RGID ${meta.id} \\
        --RGLB Library1 \\
        --RGPL Illumina \\
        --RGSM ${meta.id} \\
        --RGPU Unknown
    """
}
