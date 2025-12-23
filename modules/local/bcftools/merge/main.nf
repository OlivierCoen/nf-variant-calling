process BCFTOOLS_MERGE {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
    path(vcfs)
    path(tbis)
    tuple val(meta), path(fasta)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("*.{bcf,vcf}{,.gz}"), emit: vcf
    tuple val(meta), path("*.{csi,tbi}")      , emit: index, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "output"
    def input = (vcfs.collect().size() > 1) ? vcfs.sort{ it.name } : vcfs

    """
    bcftools merge \\
        $args \\
        --gvcf ${fasta} \\
        --threads $task.cpus \\
        --output ${prefix}.vcf.gz \\
        --merge none \\
        --write-index=csi \\
        $input
    """
}
