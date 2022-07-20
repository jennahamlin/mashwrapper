# ![nf-core/mashwrapper](docs/images/nf-core-mashwrapper_logo_light.png#gh-light-mode-only) ![nf-core/mashwrapper](docs/images/nf-core-mashwrapper_logo_dark.png#gh-dark-mode-only)

[![GitHub Actions CI Status](https://github.com/jennahamlin/nf-core-mashwrapper/workflows/nf-core%20CI/badge.svg)](https://github.com/jennahamlin/mashwrapper/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/jennahamlin/nf-core-mashwrapper/workflows/nf-core%20linting/badge.svg)](https://github.com/jennahamlin/nf-core-mashwrapper/actions?query=workflow%3A%22nf-core+linting%22)
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/mashwrapper/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23mashwrapper-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/mashwrapper)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/mashwrapper** is a wrapper around the program [Mash](https://mash.readthedocs.io/en/latest/) and the [NCBI Datasets command line tools](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/download-and-install/). Ultimately, the tool identifies the most likely species from a pair of gzipped fastq reads. The database that the reads are tested against can either be generated via the tool from a text file (--get_database) or an already built database supplied to the tool (--use_database). The output is a text file of the top five hits associated with the input reads along with standard Mash output. 

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/mashwrapper/results).

## Pipeline summary

1. Confirm input sample sheet
2. OPTIONAL (--get_database): Confirm input organism sheet
3. OPTIONAL (--get_database): Download genomes from NCBI using [NCBI datasets command line tool](https://www.ncbi.nlm.nih.gov/datasets/)
4. OPTIONAL (--get_database): Format donwloaded genomes to be Genus_Species_GenebankIdentifier.fna using [NCBI dataformat command line tool](https://www.ncbi.nlm.nih.gov/datasets/docs/v1/quickstarts/command-line-tools/#install-using-curl)
5. OPTIONAL (--get_database): Build individual [Mash sketches](https://mash.readthedocs.io/en/latest/) for all genomes downloaded
6. OPTIONAL (--get_database): Build [Mash database](https://mash.readthedocs.io/en/latest/) for all mash sketches
7. Test fastq.gz reads against either an optionally built Mash database or one provided by the user
8. Collate results from each isolate of interest tested against the Mash database

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Clone the pipeline and test it on a minimal dataset:

 >  A test dataset is available once you git clone this repo and includes [five files](https://github.com/jennahamlin/mashwrapper/tree/main/test-data):
 > - inputDB.txt - text file of species to download when using the test profile testGet (-profile testGet). Does not include a header
 > - inputReads.csv - csv file with single pair of reads listed and includes the header: sample,fastq_1,fastq_2
 > - myMashDatabase.msh - prebuilt Mash database using the same isolates listed in the inputDB.txt file. Will be used if test profile is set to testUse (-profile testUse)
 > - subERR125190_(1,2).fastq.gz - subset reads of *Legionella fallonii* to only 45000 reads

*For CDC users, you need to include the flag --custom_config_base and point to the conf subdirectory to supply the certificate information for singularity use, via this [custom config file](https://github.com/jennahamlin/mashwrapper/tree/main/conf) (e.g., --custom_config_base /scicomp/home-pure/ptx4/mashwrapper/conf)*
 

   ```console
    ## test out download of database
    nextflow run mashwrapper -profile testGet,YOURPROFILE
    
    ## test out using pre-built database
    nextflow run mashwrapper -profile testUse,YOURPROFILE 
   ```
   
4. Start running your own analysis!

  ```console
   ## Use your already built database
   nextflow run nf-core/mashwrapper -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --use_database myMashDatabase.msh --custom_config_base /scicomp/home-pure/ptx4/mashwrapper/conf

   ## Download and built your database for organism(s) of interest
   nextflow run nf-core/mashwrapper -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input samplesheet.csv --get_database organismsheet.txt --custom_config_base /scicomp/home-pure/ptx4/mashwrapper/conf
  ```

## Documentation

The nf-core/mashwrapper pipeline comes with documentation about the pipeline [usage](https://nf-co.re/mashwrapper/usage), [parameters](https://nf-co.re/mashwrapper/parameters) and [output](https://nf-co.re/mashwrapper/output).

## Credits

nf-core/mashwrapper was originally written by Jenna Hamlin.

We thank the following people for their extensive assistance in the development of this pipeline:

- Sateeshe Peri
- Micheal Cipriano

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#mashwrapper` channel](https://nfcore.slack.com/channels/mashwrapper) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/mashwrapper for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
