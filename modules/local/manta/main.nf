process MANTA_GERMLINE {
    tag "${meta.id} on ${region.chrom}:${region.start}-${region.end}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f6/f696c93e6209e33ac0d15f1ecfa799bc67329eec07b0569e065ea8b220b53953/data' :
        'community.wave.seqera.io/library/manta_python:0eb71149179b3920' }"

    input:
    tuple path(bam_files), path(bai_files), val(region)
    tuple val(meta), path(fasta), path(fai)
    path(config)

    output:
    tuple val(meta), path("*candidate_small_indels.vcf.gz")    , emit: candidate_small_indels_vcf
    tuple val(meta), path("*candidate_small_indels.vcf.gz.tbi"), emit: candidate_small_indels_vcf_tbi
    tuple val(meta), path("*candidate_sv.vcf.gz")              , emit: candidate_sv_vcf
    tuple val(meta), path("*candidate_sv.vcf.gz.tbi")          , emit: candidate_sv_vcf_tbi
    tuple val(meta), path("*diploid_sv.vcf.gz")                , emit: diploid_sv_vcf
    tuple val(meta), path("*diploid_sv.vcf.gz.tbi")            , emit: diploid_sv_vcf_tbi
    tuple val("${task.process}"), val('manta'), eval("manta: \$( configManta.py --version )"), topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_files = bam_files.collect{"--bam ${it}"}.join(' ')
    def config_option = config ? "--config ${config}" : ""
    """
    echo "${region.chrom}\t${region.start}\t${region.end}" > region.bed

    configManta.py \\
        ${input_files} \\
        ${config_option} \\
        --reference $fasta \\
        --runDir manta \\
        --callRegions region.bed \\
        $args

    python manta/runWorkflow.py -m local -j $task.cpus

    mv manta/results/variants/candidateSmallIndels.vcf.gz \\
        ${prefix}.candidate_small_indels.vcf.gz
    mv manta/results/variants/candidateSmallIndels.vcf.gz.tbi \\
        ${prefix}.candidate_small_indels.vcf.gz.tbi
    mv manta/results/variants/candidateSV.vcf.gz \\
        ${prefix}.candidate_sv.vcf.gz
    mv manta/results/variants/candidateSV.vcf.gz.tbi \\
        ${prefix}.candidate_sv.vcf.gz.tbi
    mv manta/results/variants/diploidSV.vcf.gz \\
        ${prefix}.diploid_sv.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi \\
        ${prefix}.diploid_sv.vcf.gz.tbi
    """
}
