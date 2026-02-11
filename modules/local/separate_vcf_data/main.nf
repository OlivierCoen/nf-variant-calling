process SEPARATE_VCF_DATA {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e65b93cd90735298b93ea6864cb7072e839b88c9509225ed876b674a2c16666/data':
        'community.wave.seqera.io/library/polars_python:f73377f543756137' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${prefix}.variants.parquet"),                                                        emit: variants
    tuple val(meta), path("${prefix}.RO_counts.parquet"),                                                       emit: ref_counts
    tuple val(meta), path("${prefix}.AO_counts.parquet"),                                                       emit: alt_counts
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"),                 topic: versions
    tuple val("${task.process}"), val('polars'), eval('python3 -c "import polars; print(polars.__version__)"'), topic: versions

    script:
    prefix = "${meta.variant_type}"
    """
    # limiting number of threads
    export POLARS_MAX_THREADS=${task.cpus}

    separate_vcf_data.py \\
        --vcf $vcf

    mv variants.parquet ${prefix}.variants.parquet
    mv RO_counts.parquet ${prefix}.RO_counts.parquet
    mv AO_counts.parquet ${prefix}.AO_counts.parquet
    """

}
