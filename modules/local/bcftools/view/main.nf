process BCFTOOLS_VIEW {
    tag "${meta.id} - ${meta.type}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val biallelic_only
    val min_qual
    val min_overall_depth
    val extra_filters

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"),                                                       emit: vcf_tbi
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'"), topic: versions

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}.filtered"
    def biallelic_args = biallelic_only ? "-m 2 -M 2" : ""
    def qual_filter = min_qual ? "QUAL>=${min_qual}" : ""
    def depth_filter = min_overall_depth ? "INFO/DP>=${min_overall_depth}" : ""
    def filters = [qual_filter, depth_filter, extra_filters].findAll { it }.join(" & ")
    def filter_args = filters ? "-i '${filters}'" : ""
    """
    bcftools view \\
        --output-type z \\
        --output ${prefix}.vcf.gz \\
        --write-index=tbi \\
        --threads ${task.cpus} \\
        $biallelic_args \\
        $filter_args \\
        $args \\
        $vcf
    """
}
