process FREEBAYES {

    tag "${meta.id} on ${region.chrom}:${region.start}-${region.end}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freebayes:1.3.10--hbefcdb2_0'
        : 'biocontainers/freebayes:1.3.10--hbefcdb2_0'}"

    input:
    tuple val(meta), path(bam), path(bai), val(region)
    tuple val(meta2), path(fasta), path(fai)

    output:
    path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('freebayes'), eval("freebayes --version 2>&1 | sed 's/version:\s*v//g'"),    topic: versions


    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${region.chrom}_${region.start}_${region.end}"
    """
    freebayes \\
        --fasta-reference ${fasta} \\
        --gvcf \\
        ${args} \\
        --bam ${bam} \\
        --region ${region.chrom}:${region.start}-${region.end} \\
        --vcf ${prefix}.vcf

    bgzip ${prefix}.vcf
    """
}
