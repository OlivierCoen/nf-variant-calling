# nf-variant-calling

[![GitHub Actions CI Status](https://github.com/OlivierCoen/nf-variant-calling/actions/workflows/nf-test.yml/badge.svg)](https://github.com/OlivierCoen/nf-variant-calling/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/OlivierCoen/nf-variant-calling/actions/workflows/linting.yml/badge.svg)](https://github.com/OlivierCoen/nf-variant-calling/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/OlivierCoen/nf-variant-calling)

## Introduction

**nf-variant-calling** is a bioinformatics pipeline that performs variant calling on using short reads and a genome fasta file.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Read samplesheet

First, prepare a samplesheet listing your input reads and their associated metadata.

`samplesheet.csv`:

```csv
sample,population,phenotype,lane,fastq_1,fastq_2
control,pop1,phenotype1,1,control_R1.fastq.gz,control_R2.fastq.gz
control,pop1,phenotype1,2,control_R1.fastq.gz,control_R2.fastq.gz
sample1,pop2,phenotype2,1,sample1_R1.fastq.gz,sample1_R2.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

> [!NOTE]
> The `lane` and `fastq_2̀` fields are optional. If `fastq_2̀` is not provided, files will be interpreted as single-end reads.

### Pipeline parameters

Then, you must prepare a `params.yaml` file with **at least** the following content:

```yaml
input: <PATH TO SAMPLESHEET>
genome: <PATH TO GENOME FASTA FILE>
outdir: <OUTDIR>
```

> [!NOTE]
> The full list of parameters can be obtained with `nextflow run OlivierCoen/nf-variant-calling --help`.

### Advanced usage: custom config

Finally, you can prepare a `custom.config` file if you want to tweak the pipeline configuration.
See [the configuration documentation](docs/configuration.md) for more details.


### Run the pipeline

Now, you can run the pipeline using:

```bash
nextflow run OlivierCoen/nf-variant-calling \
   -latest \
   -profile <docker/apptainer/singularity/conda/...> \
   -params-file params.yaml \
   [-config custom.config \]
   -resume
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

OlivierCoen/nf-variant-calling was originally written by Olivier Coen.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use OlivierCoen/nf-variant-calling for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
