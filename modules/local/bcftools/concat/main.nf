process BCFTOOLS_CONCAT {
    tag "${meta.id} - ${meta.variant_type}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(inputs), path(indexes)

    output:
    tuple val(meta), path("${prefix}.bcf.gz"), emit: bcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'"), topic: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    bcftools concat \\
        --output-type b \\
        --allow-overlaps \\
        --output ${prefix}.bcf.gz \\
        $args \\
        --threads ${task.cpus} \\
        ${inputs}
    """
}
