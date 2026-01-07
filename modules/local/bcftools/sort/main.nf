process BCFTOOLS_SORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${prefix}.bcf.gz"), emit: bcf
    tuple val(meta), path("*.csi")           , emit: csi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'"), topic: versions

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools \\
        sort \\
        --output-type b \\
        --output ${prefix}.bcf.gz \\
        --temp-dir . \\
        --write-index=csi \\
        ${args} \\
        ${input}
    """
}
