# ghga-de/aqua: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to the top-level results directory. The pipeline organizes outputs primarily by sample name, with subdirectories for each quality control tool executed.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data through various quality control steps. The results are summarized in a final [MultiQC](#multiqc) report.

## Sample-specific Outputs

Each sample has a dedicated directory named after its sample identifier. Inside these folders, you will find tool-specific results depending on the input type (FastQ, BAM, or VCF) and analysis method.

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the report data and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content, adapter contamination and overrepresented sequences.

### Fastp

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/fastp/`
  - `*.fastp.html`: Visual report of sequencing stats and adapter trimming.
  - `*.fastp.json`: Machine-readable summary of the fastp results.
  - `*.fastp.log`: Execution logs for the tool.

</details>

[Fastp](https://github.com/OpenGene/fastp) is an ultra-fast tool that performs basic quality control and provides detailed metrics on read quality before and after any filtering steps.

### Picard CollectMultipleMetrics

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/picard/`
  - `*.alignment_summary_metrics`: Detailed alignment statistics.
  - `*.insert_size_metrics`: Calculated insert size statistics for paired-end data.
  - `*.quality_distribution.pdf`: Plot showing the distribution of base qualities.
  - `*.base_distribution_by_cycle.pdf`: Plot showing nucleotide percentages across read positions.

</details>

[Picard](https://broadinstitute.github.io/picard/) is used for alignment quality assessment. It generates comprehensive metrics regarding mapping rates, duplication, and read distribution.

### Mosdepth

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/mosdepth/`
  - `*.mosdepth.summary.txt`: Summary of average coverage across chromosomes.
  - `*.mosdepth.global.dist.txt`: Cumulative coverage distribution data.

</details>

[Mosdepth](https://github.com/brentp/mosdepth) provides extremely fast coverage calculations. It determines how well the target regions were covered by sequencing and identifies potential gaps in data.

### BCFtools Stats

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/bcftools/`
  - `*.bcftools_stats.txt`: General statistics for variant call files.

</details>

[BCFtools](http://samtools.github.io/bcftools/bcftools.html) generates statistics for VCF files, including the number of SNPs, indels, transitions, and transversions.

### NanoPlot (Long Reads)

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/nanoplot/`
  - `NanoPlot-report.html`: Visual report of long-read sequencing quality.
  - `NanoStats.txt`: General statistics for the sequencing run.
  - `*.png`: Various plots including length vs quality scatter plots.

</details>

[NanoPlot](https://github.com/wdecoster/NanoPlot) is used specifically for long-read data (Nanopore or PacBio) to visualize read lengths and quality scores.

### SeqFU

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/seqfu/`
  - `*.tsv`: Tab-separated file containing sequence statistics like N50 and read counts.

</details>

### RSeQC

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/rseqc/`
  - `*.bam_stat.txt`: Summary of mapping statistics for RNA-Seq data.

</details>

### Samtools

<details markdown="1">
<summary>Output files</summary>

- `<SAMPLE_ID>/samtools/`
  - `*.bam.bai`: Indices for input BAM files.
  - `*.fastq.gz`: FastQ files extracted from BAM inputs where applicable.
  - `genome.fa.fai`: Index for the reference genome used.

</details>

## Global Aggregate Outputs

These directories contain reports and files that combine data from across the entire run.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: The combined visual report for all samples.
  - `multiqc_data/`: Directory containing parsed statistics used to generate the report.
  - `multiqc_plots/`: Static images of report plots in PDF, PNG, and SVG formats.

</details>

[MultiQC](http://multiqc.info) aggregates results from every QC tool used in the pipeline into a single report. This is the primary file to check for an overall assessment.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `params_*.json`: Records of the parameters used for each run.
  - `qc_software_mqc_versions.yml`: List of all tool versions for reproducibility.

</details>
