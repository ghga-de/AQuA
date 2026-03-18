# ghga-de/aqua: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

The ghga-de/aqua pipeline is designed for comprehensive quality control of genomic and transcriptomic data. It supports multiple input formats, including raw sequencing reads (FASTQ) and aligned sequences (BAM/CRAM), as well as variant calls (VCF). The pipeline adaptively selects QC tools based on the input type and sequencing method specified in the samplesheet.

## Metadata JSON input

If you are a GHGA user, you can also provide a metadata JSON file directly. The pipeline will automatically convert this into the required samplesheet format:

```bash
nextflow run main.nf --metadata ./metadata.json --outdir ./results -profile docker
```

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location.

```bash
--input '[path to samplesheet file]'
```

The samplesheet is a comma-separated (CSV) file. While the pipeline is flexible, every row must contain at least a sample name and exactly one type of data input (FASTQ pair, BAM, CRAM, or VCF).

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

## Flexible Data Inputs

The pipeline supports three primary "steps" or entry points based on the files provided:

1. Raw Reads (Step 1): Provide fastq_1 (and optionally fastq_2).
2. Alignments (Step 2): Provide bam and bai (or cram and crai).
3. Variants (Step 3): Provide vcf.

## Entry Points and Parameters

The pipeline adjusts its behavior based on the content of your samplesheet:

WGS/WXS: Runs Picard, Samtools, and Mosdepth.

RNA-Seq: Runs RSeQC and Picard RNA metrics.

Long-reads: If experiment_method is nanopore or pacbio, specialized tools like NanoPlot or FastpLong are triggered.

### Full samplesheet

A samplesheet containing a mix of raw data, mapped bams, and variant files would look like this:

```csv title="samplesheet.csv"
sample,lane,individual_id,sex,experiment_method,fastq_1,fastq_2,bam,bai,vcf
SAMPLE_FASTQ,L001,ind_1,MALE,wgs,s1_R1.fastq.gz,s1_R2.fastq.gz,,,
SAMPLE_BAM,L001,ind_2,FEMALE,wgs,,,s2.bam,s2.bam.bai,
SAMPLE_VCF,L001,ind_3,NA,wgs,,,,,s3.vcf.gz
```

| Column                | Description                                                                                                                                                                        |
| :-------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`              | **Required.** Custom sample name. This identifier is used to group multiple sequencing runs (lanes) from the same sample. Spaces are automatically converted to underscores (`_`). |
| `lane`                | **Required.** identifier for the sequencing lane or library (e.g., L001, L002). Must not contain spaces.                                                                           |
| `individual_id`       | Identifier for the individual (patient/subject).                                                                                                                                   |
| `sex`                 | Biological sex of the individual (e.g., MALE, FEMALE, NA).                                                                                                                         |
| `status`              | Disease status as an integer: `0` (Normal/Control) or `1` (Tumor/Case).                                                                                                            |
| `phenotype`           | Phenotypic terms or descriptions associated with the individual.                                                                                                                   |
| `sample_type`         | The type of sample (e.g., GENOMIC_DNA, TOTAL_RNA).                                                                                                                                 |
| `disease_status`      | Text description of the disease status (e.g., Healthy, Tumor).                                                                                                                     |
| `case_control_status` | Status in the study design (e.g., CASE, CONTROL).                                                                                                                                  |
| `tissue`              | The source tissue of the specimen (e.g., blood, tissue).                                                                                                                           |
| `experiment_method`   | The sequencing method used. Supported values: `wgs`, `wes`, `rna`, `atac`, `nanopore`, `pacbio`.                                                                                   |
| `analysis_method`     | The type of analysis performed (e.g., `varcall`).                                                                                                                                  |
| `fastq_1`             | Path to the Read 1 FastQ file. Must end in `.fastq.gz` or `.fq.gz`.                                                                                                                |
| `fastq_2`             | Path to the Read 2 FastQ file for paired-end data. Optional for single-end.                                                                                                        |
| `single_end`          | Boolean (`true`/`false`) indicating if the sequencing is single-end.                                                                                                               |
| `bam`                 | Path to the aligned BAM file.                                                                                                                                                      |
| `bai`                 | Path to the corresponding BAM index file.                                                                                                                                          |
| `cram`                | Path to the aligned CRAM file.                                                                                                                                                     |
| `crai`                | Path to the corresponding CRAM index file.                                                                                                                                         |
| `vcf`                 | Path to the Variant Call Format file. Must end in `.vcf` or `.vcf.gz`.                                                                                                             |
| `data_files`          | Semicolon-separated list of any other relevant data files not covered by specific columns.                                                                                         |

An [example samplesheet](../tests/samplesheets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run main.nf --input ./samplesheet.csv --outdir ./results --genome GRCh38 -profile singularity
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run ghga-de/aqua -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull ghga-de/aqua
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [ghga-de/aqua releases page](https://github.com/ghga-de/aqua/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
