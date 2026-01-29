process FREEBAYES {

    tag "${meta.id} on ${region.chrom}:${region.start}-${region.end}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/freebayes:1.3.10--hbefcdb2_0'
        : 'biocontainers/freebayes:1.3.10--hbefcdb2_0'}"

    input:
    tuple path(bam_files), path(bai_files), val(region)
    tuple val(meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('freebayes'), eval("freebayes --version 2>&1 | sed 's/version:\s*v//g'"), topic: versions


    script:
    def args    = task.ext.args   ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}_${region.chrom}_${region.start}_${region.end}"
    """
    # freebayes uses only 1 core
    # that's why we split the genome in multiple chunks and run freebayes in parallel on these chunks
    freebayes \\
        --fasta-reference ${fasta} \\
        ${args} \\
        --bam ${bam_files} \\
        --region ${region.chrom}:${region.start}-${region.end} \\
        --vcf ${prefix}.vcf

    bgzip ${prefix}.vcf
    """
}
