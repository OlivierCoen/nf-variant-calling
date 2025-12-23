process BCFTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta),  path(vcf), path(tbi)
    tuple val(meta6), path(fasta), path(fai)

    output:
    tuple val(meta), path("*stats.txt"), topic: bcftools_stats
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_fasta = fasta ? "--fasta-ref ${fasta}" : ""
    """
    bcftools stats \\
        $args \\
        $reference_fasta \\
        $vcf > ${prefix}.bcftools_stats.txt
    """
}
