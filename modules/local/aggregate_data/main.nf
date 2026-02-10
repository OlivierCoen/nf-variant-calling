process AGGREGATE_DATA {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebf9af86c494e5b3e93706ae0d611d28e8788b2e99717f87a1f586ecf30b7067/data':
        'community.wave.seqera.io/library/polars_python_tqdm:ca595df92ae9b061' }"

    input:
    tuple val(meta), path(variant_file), path(pvalue_file), path(ref_counts), path(alt_counts)
    path(design_file)
    val(window_size)

    output:
    tuple val(meta), path("*.parquet"),                                                                         emit: grouped_variants
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions
    tuple val("${task.process}"), val('tqdm'),   eval('python3 -c "import tqdm; print(tqdm.__version__)"'),     topic: versions

    script:
    prefix = task.ext.prefix ?: "${meta.variant_type}.grouped_variants"
    """
    # limiting number of threads
    export POLARS_MAX_THREADS=${task.cpus}

    aggregate_data.py \\
        --variants $variant_file \\
        --RO $ref_counts \\
        --AO $alt_counts \\
        --design $design_file \\
        --out ${prefix}.parquet \\
        --pvalues $pvalue_file \\
        --window-size $window_size
    """

}
