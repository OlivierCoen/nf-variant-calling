process SPLIT_FASTA_WITH_OVERLAPS {

    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/76/7604ee1c2b28c6c949061fa659fc5cb55e4dbbf17aec15d5f9d623ac654a077e/data':
        'community.wave.seqera.io/library/biopython:1.86--33dad2a816ac5c52' }"

    input:
    tuple val(meta), path(fasta)
    val chunksize

    output:
    tuple val(meta), path("chunks/*.fasta"),                                                                                  emit: chunks
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                         topic: versions
    tuple val("${task.process}"), val('biopython'),    eval('python3 -c "import biopython; print(biopython.__version__)"'),   topic: versions

    script:
    """
    split_fasta_with_overlaps.py \\
        --fasta $fasta \\
        --chunksize $chunksize \\
        --outdir chunks
    """

}
