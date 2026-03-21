process SURVIVOR_MERGE {

    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/48c985f5c905e20a0298c5e2a937451b9388ac7d2f4ec860243acc7e96314678/data':
        'community.wave.seqera.io/library/survivor_tabix:65c4fc27f50bed38' }"

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.merged_matrix.txt"), emit: merged_matrix
    tuple val("${task.process}"), val('survivor'), eval("SURVIVOR 2>&1 | grep 'Version' | sed 's/Version: //'"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.survivor_merged"
    """
    SURVIVOR merge \\
        <(ls *.vcf) \\
        $args \\
        ${prefix}.vcf

    SURVIVOR genComp ${prefix}.vcf ${prefix}.merged_matrix.txt

    bgzip --threads ${task.cpus} ${prefix}.vcf.gz
    """
}
