process ADDITIONAL_FILTERING {

    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/01/01098af425b76fe5d38303e0462456beff5cc9c3eee56b4a7ebacb0d0bad7754/data':
        'community.wave.seqera.io/library/htslib_matplotlib_polars_python:b8381ca8e371940c' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), optional: true,                                                                emit: vcf
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                         topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),         topic: versions
    tuple val("${task.process}"), val('matplotlib'),   eval('python3 -c "import matplotlib; print(matplotlib.__version__)"'), topic: versions

    script:
    def args = task.ext.args ?: ""
    prefix = task.ext.prefix ?: "${meta.id}.filtered"
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    export MPLCONFIGDIR=\$PWD

    apply_additional_filters.py \\
        --vcf $vcf \\
        --fai $fai \\
        --out ${prefix}.vcf \\
        ${args}

    if [[ -f ${prefix}.vcf ]]; then
        bgzip ${prefix}.vcf
    fi
    """

}
