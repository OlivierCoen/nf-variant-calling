process FREEBAYES {

    tag "${meta.id} on ${chrom}:${start}-${end}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freebayes:1.3.10--hbefcdb2_0'
        : 'biocontainers/freebayes:1.3.10--hbefcdb2_0'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(fasta), path(fai)
    tuple val(chrom), val(start), val(end)

    output:
    path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('freebayes'), eval("freebayes --version 2>&1 | sed 's/version:\s*v//g'"),    topic: versions


    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${chrom}_${start}_${end}"
    """
    freebayes \\
        --fasta-reference ${fasta} \\
        --gvcf \\
        ${args} \\
        --bam ${bam} \\
        --region ${chrom}:${start}-${end} \\
        --vcf ${prefix}.vcf

    bgzip ${prefix}.vcf
    """
}
