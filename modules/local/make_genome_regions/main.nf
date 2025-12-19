process MAKE_GENOME_REGIONS {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/76/7604ee1c2b28c6c949061fa659fc5cb55e4dbbf17aec15d5f9d623ac654a077e/data':
        'community.wave.seqera.io/library/biopython:1.86--33dad2a816ac5c52' }"

    input:
    tuple val(meta), path(genome)
    val chunksize

    output:
    path 'genome_regions.tsv',                                                                                        emit: regions
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                     topic: versions
    tuple val("${task.process}"), val('Bio'),      eval('python3 -c "import Bio; print(Bio.__version__)"'),           topic: versions

    script:
    def chunksize_int = chunksize.toInteger()
    """
    make_genome_regions.py \\
        --genome $genome \\
        --chunksize $chunksize_int
    """

}
