process EVALUATE_EFFECT_OF_FILTERS {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/01/01098af425b76fe5d38303e0462456beff5cc9c3eee56b4a7ebacb0d0bad7754/data':
        'community.wave.seqera.io/library/htslib_matplotlib_polars_python:b8381ca8e371940c' }"

    input:
    tuple val(meta), path(vcf), path(filtered_vcf)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("${prefix}.png"),                                                                                   emit: png
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                         topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),         topic: versions
    tuple val("${task.process}"), val('matplotlib'),   eval('python3 -c "import matplotlib; print(matplotlib.__version__)"'), topic: versions

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.variant_type}.effect_of_filters_on_variants"
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads
    export POLARS_MAX_THREADS=${task.cpus}

    export MPLCONFIGDIR=\$PWD

    evaluate_effect_of_filters.py \\
        --vcf $vcf \\
        --filtered-vcf $filtered_vcf \\
        --fai $fai \\
        --out ${prefix}.png
    """

}
