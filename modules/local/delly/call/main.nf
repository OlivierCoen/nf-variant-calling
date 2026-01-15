process DELLY_CALL {

    tag "${meta.id} on ${region.chrom}:${region.start}-${region.end}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d8/d871bbc11d83f49db1090b0613d4dc2244821864e280f092d4bc5c5121a8db1c/data' :
        'community.wave.seqera.io/library/delly:1.7.2--a8cf3dbcfe03c0dd' }"

    input:
    tuple val(meta), path(bam), path(bai), val(region)
    tuple val(meta2), path(fasta), path(fai)
    path(all_regions_file)

    output:
    tuple val(meta), path("*.vcf.gz"),    emit: vcf
    tuple val(meta), path("*.{csi,tbi}"), emit: csi
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
        --threads ${task.cpus} \\
        ${args} \\
        --genome ${fasta} \\
        ${bam} \\
        | bgzip ${args2} --threads ${task.cpus} --stdout \\
        > ${prefix}.vcf.gz && tabix ${prefix}.vcf.gz
    """
}
