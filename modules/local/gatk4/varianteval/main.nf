process GATK4_VARIANTEVAL {
    tag "${fasta}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta), path(fai)

    output:
    path('output.eval.grp'), emit: output
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version 2>&1 | sed 's/^.*(GATK) v//; s/ .*\$//'"), topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${vcf.baseName}"
    def avail_mem = 6144
    if (!task.memory) {
        log.info('[GATK VariantEval] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantEval \\
        --REFERENCE ${fasta} \\
        --OUTPUT ${prefix}.eval.grp \\
        --eval set:${vcf} \\
        ${args}
    """
}
