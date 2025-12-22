process DELLY_CALL {

    tag "${meta.id} on ${region.chrom}:${region.start}-${region.end}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/delly:1.3.3--h4d20210_0' :
        'biocontainers/delly:1.3.3--h4d20210_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(region)
    tuple val(meta2), path(fasta), path(fai)
    path(all_regions_file)

    output:
    tuple val(meta), path("*.vcf.gz")  , emit: vcf
    tuple val(meta), path("*.{csi,tbi}")     , emit: csi
    tuple val("${task.process}"), val('delly'), eval("delly --version 2>&1 | sed 's/^.*Delly version: v//; s/ using.*\$//'"), topic: versions


    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${region.chrom}_${region.start}_${region.end}"
    """
    # making file for regions to exclude
    egrep -v "^${region.chrom}\s+${region.start}\s+${region.end}\$" ${all_regions_file} > ${prefix}.regions.bed

    delly \\
        call \\
        ${args} \\
        --genome ${fasta} \\
        ${bam} \\
        | bgzip ${args2} --threads ${task.cpus} --stdout \\
        > ${prefix}.vcf.gz && tabix ${prefix}.vcf.gz
    """
}
