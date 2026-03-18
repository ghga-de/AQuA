# GHGA AQuA Nextflow Pipeline

[![GitHub Actions CI Status](https://github.com/ghga-de/aqua/actions/workflows/nf-test.yml/badge.svg)](https://github.com/ghga-de/aqua/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/ghga-de/aqua/actions/workflows/linting.yml/badge.svg)](https://github.com/GHGA/qc/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<p align="center">
    <img title="GHGA_logo" src="docs/GHGA_short_Logo_orange.png" width=50%>
</p>

## Introduction

**GHGA AQuA (Automatic Quality Assessment) Pipeline** is a bioinformatics pipeline that performs basic quality control over input datasets without altering the raw data. It accepts three main input types:

1. Raw FastQ files
2. Aligned BAM/CRAM files
3. Variant called VCF/BCF files

The pipeline automatically selects the appropriate quality control tools based on your provided analysis method (e.g., WGS, RNA-Seq, Nanopore) and compiles the results into a single MultiQC report.

## Tool & Analysis Matrix

The following table details which tools are executed based on the analysis method and input data type provided in the samplesheet.

| Analysis Method     | Read QC (FastQ)      | Alignment QC (BAM/CRAM)                                   | Variant QC (VCF) |
| :------------------ | :------------------- | :-------------------------------------------------------- | :--------------- |
| **WGS / WES / TES** | FastQC, FastP, SeqFU | Mosdepth, Samtools Stats, Picard, VerifyBamID, NGS-Bits\* | BCFTools Stats   |
| **ATAC / ChIP-Seq** | FastQC, FastP, SeqFU | Mosdepth, Samtools Stats, Picard                          | BCFTools Stats   |
| **RNA-Seq / smRNA** | FastQC, FastP, SeqFU | RSeQC                                                     | BCFTools Stats   |
| **Nanopore**        | FastQC, NanoPlot     | -                                                         | BCFTools Stats   |
| **PacBio**          | FastPLong            | -                                                         | BCFTools Stats   |
| **MethylSeq**       | FastQC, FastP, SeqFU | -                                                         | BCFTools Stats   |

> \* _NGS-Bits SampleGender is run for WGS if `predict_sex` is enabled._

## Usage

The pipeline can be started in two ways: by providing a manual samplesheet or by providing GHGA-compliant metadata.

### Option A: Manual Samplesheet (CSV)

Create a samplesheet.csv with your data. The pipeline auto-detects the starting step based on which columns are populated.

You must create a `samplesheet.csv` containing your input data. The structure requires a `step` column to tell the pipeline which type of file you are providing:

- **step 1**: FastQ files (Read QC)
- **step 2**: BAM/CRAM files (Alignment QC)
- **step 3**: VCF files (Variant QC)

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

### Option B: GHGA Metadata (JSON)

If you already have a metadata.json following the GHGA metadata model, you can provide it directly. The pipeline will automatically convert the JSON into the required internal format, eliminating the need to create a manual samplesheet.

### 2. Run the Pipeline

Run the pipeline using the command below with input samplesheet.csv.

```bash
nextflow run main.nf \
   -profile docker/singularity/conda \
   --input samplesheet.csv \
   --outdir ./results
```

or using the command below with input metadata.json.

```bash
nextflow run main.nf \
   -profile docker/singularity/conda \
   --metadata metadata.json \
   --outdir ./results
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Supported Tools

### Read QC

- [**FastQC**](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/): Comprehensive quality control checks on raw sequence data.
- [**FastP**](https://github.com/OpenGene/fastp): All-in-one FASTQ preprocessor (used here for QC metrics).
- [**SeqFU**](https://telatin.github.io/seqfu2/tools/metadata.html): Sequence statstics
- [**NanoPlot**](https://github.com/wdecoster/NanoPlot): Plotting tool for long read sequencing data and alignments.
- [**FastPLong**](https://github.com/OpenGene/fastplong): Quality control for long read data (PacBio).

### Alignment QC

- [**Mosdepth**](https://github.com/brentp/mosdepth): Fast BAM/CRAM depth calculation.
- [**Samtools Stats**](http://www.htslib.org/doc/samtools.html): General statistics for alignment files.
- [**Picard CollectMultipleMetrics**](https://broadinstitute.github.io/picard/): Collects multiple classes of metrics from alignment files.
- [**RSeQC**](http://rseqc.sourceforge.net/): Quality control for RNA-seq experiments.
- [**NGS-Bits SampleGender**](https://github.com/imgag/ngs-bits): Sex determination based on coverage.
- [**VerifyBamID**](https://github.com/Griffan/VerifyBamID): A robust tool for DNA contamination estimation from sequence reads using ancestry-agnostic method.

### Variant QC

- [**BCFTools Stats**](http://samtools.github.io/bcftools/bcftools.html): Statistics for VCF/BCF files.

### Reporting

- [**MultiQC**](http://multiqc.info/): Aggregates results from all tools into a single HTML report.

## Credits

GHGA AQuA nextflow pipeline was originally written by Kubra Narci @kubranarci.

Current development team: - Manuel Kösters - Virag Sharma - Ruchi Tanavade

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
