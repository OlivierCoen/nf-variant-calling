process GET_ALLELE_COUNTS {

    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e65b93cd90735298b93ea6864cb7072e839b88c9509225ed876b674a2c16666/data':
        'community.wave.seqera.io/library/polars_python:f73377f543756137' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("RO_counts.csv"), path("AO_counts.csv"),                                                            emit: counts
    tuple val("${task.process}"), val('python'),       eval("python3 --version | sed 's/Python //'"),                         topic: versions
    tuple val("${task.process}"), val('polars'),       eval('python3 -c "import polars; print(polars.__version__)"'),         topic: versions

    script:
    def is_using_containers = workflow.containerEngine ? true : false
    """
    # limiting number of threads when using conda / micromamba
    if [ "${is_using_containers}" == "false" ]; then
        export POLARS_MAX_THREADS=${task.cpus}
    fi

    get_allele_counts.py \\
        --vcf $vcf
    """

}
