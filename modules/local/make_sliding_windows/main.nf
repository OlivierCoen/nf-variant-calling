process MAKE_SLIDING_WINDOWS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e65b93cd90735298b93ea6864cb7072e839b88c9509225ed876b674a2c16666/data':
        'community.wave.seqera.io/library/polars_python:f73377f543756137' }"

    input:
    tuple val(meta), path(vcf), path(pvalue_file)
    val(window_size)

    output:
    tuple val(meta), path("*.parquet"),                                                                                       emit: grouped_variants
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                         topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),         topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.variant_type}.grouped_by_windows"
    """
    # limiting number of threads
    export POLARS_MAX_THREADS=${task.cpus}

    make_siding_windows.py \\
        --vcf $vcf \\
        --out ${prefix}.parquet \\
        --pvalues $pvalue_file \\
        --window-size $window_size
    """

}
