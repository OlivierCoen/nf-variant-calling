process VCF_ADDTIONAL_FILTERS {

    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b5/b51efb1ddfde9f1ab70471cd39d363963434e61ea4f40afe51b57eb945412a66/data':
        'community.wave.seqera.io/library/tabix_matplotlib_polars_python:60e45a5643dbfbbc' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(genome_fai)
    val min_depth
    val max_depth_quantile


    output:
    tuple val(meta), path('*.filtered.vcf'),                                                                                emit: vcf
    tuple val("${task.process}"), val('python'),   eval("python3 --version | sed 's/Python //'"),                                topic: versions
    tuple val("${task.process}"), val('polars'),      eval('python3 -c "import polars; print(polars.__version__)"'),             topic: versions
    tuple val("${task.process}"), val('matplotlib'),      eval('python3 -c "import matplotlib; print(matplotlib.__version__)"'), topic: versions

    script:
    def chunksize_int = chunksize.toInteger()
    """
    apply_additional_filters.py \\
        --vcf $vcf \\
        --fai $genome_fai \\
        --min-depth $min_depth \\
        --max-depth-quantile $max_depth_quantile \\
        --chunksize $chunksize_int
    """

}
