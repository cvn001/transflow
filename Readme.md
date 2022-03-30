# 💎 TransFlow: Tuberculosis Transmission Analysis Workflow

## 🧉 Introduction

TransFlow is a modular, flexible and user-friendly tuberculosis transmission analysis workflow based on whole genome sequencing (WGS) of *Mycobacterium tuberculosis* complex (MTBC).

The workflow filters non-MTBC samples using Kraken, then preforms quality contol (QC) using both FastQC and MultiQC. After that, it uses the PANPASCO workflow to do pan-genome mapping and relative pairwise SNP distance calculation for transmission analysis. Next, it infers tranmission clusters and networks using transcluster and SeqTrack, separately. Finally, it detects risk factors that significantly associate with transmission clustering.

## 🐍 Workflow

![Workflow](https://github.com/cvn001/transflow/blob/master/flowchart/flowchart.jpg)

## ⚙️ Installation

### Environment

### Conda

Conda can function as a package manager and is available [here](https://docs.conda.io/en/latest/miniconda.html). If you have conda make sure the bioconda and conda-forge channels are added:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Clone the repository

```bash
git clone https://gitee.com/cwmda/transflow.git
```

### Create the environment

```bash
conda env create -f workflows/envs/env.yaml
```

### Activate the environment

```bash
conda activate transflow
```

### Install the **transcluster** in **R**

```R
devtools::install_github("JamesStimson/transcluster")
```

### Kraken database

A pre-built 8 GB database [MiniKraken DB_8GB](https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_8GB.tgz) is the suggested reference database for TransFlow. It is constructed from complete bacterial, archaeal, and viral genomes in RefSeq.

## ⚙️ Set up configuration

To run the complete workflow do the following:

* Create an metadata file for all your samples with one ID per line and corresponding collection date and other characteristics.
* Copy all `.fastq.gz` files of your samples into one directory (ID_1.fastq.gz, ID_2.fastq.gz)
* Customize the workflow based on your need in: `config/configfile.yaml`
  * `metadata_file` to `/path/to/metadata_file`
  * `genome_file` to `/path/to/pangenome/pangenome_consensus.fasta`
  * `genomegaps_file` to `/path/to/pangenome/pangenome_genomegaps.bed`
  * `kraken_db` to `/path/to/minikraken_20171019_8GB`

## 🔧 Quick start
### Example

We provide an example dataset including WGS and corresponding epidemiological data which is publicly hosted at Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6345888.svg)](https://doi.org/10.5281/zenodo.6345888)

Further information is provided in the database section below.

```bash
wget https://zenodo.org/record/6345888/files/example_data.zip
upzip example_data.zip
```

Run whole pipeline in just one command:

```bash
bash transflow.sh --configfile config.yaml -j 40
```

## 🔧 Step by step

### 1. Quality control

```bash
snakemake quality_control.snakefile --configfile config.yaml -j 40
```

### 2. MTBC filtering

```bash
snakemake mtbc_filtering.snakefile --configfile config.yaml -j 40
```

### 3. Reads mapping and variant calling

```bash
snakemake variant_calling.snakefile --configfile config.yaml -j 40
```

### 4. transmission analysis

```bash
snakemake transmission_analysis.snakefile --configfile config.yaml -j 40
```

### 5. generating summary report

```bash
snakemake generating_report.snakemake --configfile config.yaml -j 40
```

### Note

#### Crashed and burned (Unlocking)

After the workflow was killed (Snakemake didn’t shutdown), the workflow directory will be still locked. If you are sure, that snakemake is no longer running `(ps aux | grep snake)`.

Unlock the working directory:

```bash
snakemake *.snakemake --unlock
```

### Rerun incomplete

If Snakemake marked a file as incomplete after a crash, delete and produce it again.

```bash
snakemake *.snakemake --ri
```

## 🏗️ Parameters

These parameters are set in the configuration files.

| Parameter | Description | Default |
| :------ | :------ | :---------- |
| genome_file | Reference genome file | *pan-genome reference provided in this repository* |
| metadata_file | File with list of samples with one ID per line| - |
| fastqdir | Directory with .fastq.gz files named ID_1.fastq.gz, ID_2.fastq.gz| fastq |
| fastqpostfix | Specification of fastq.gz format; e.g. for the format sample_R1.fastq.gz put "fastqpostfix: R" | - |
| kraken_cutoff | Threshold of MTBC reads percentage | 90 |
| allele_frequency_threshold | Allele frequency threshold for definition of high-quality SNPs  | 0.75 |
| mapping_quality_threshold | Minimum Mapping Quality for reads to be used with GATK HaplotypeCaller | 10 |  
| depth_threshold | Minimum coverage depth for high-quality regions | 5 |
| flash_overlap | Number of overlapping basepairs of paired reads for FLASH to merge reads | 10 |
| trimmomatic_read_minimum_length | Minimum length of reads to be kept after trimming | 50 |
| trimmomatic_qual_slidingwindow | Minimum quality in sliding window in trimming | 15 |
| trimmomatic_adapter_file | File with adapters for trimming | *provided with installation* |
| output_prefix | Prefix for all distance files | all_samples |
| genomegaps_file| File in .bed format with gaps in WGA of pan-genome| *provided in this repository* |
| exclude_regions_file | File in .bed format with positions that should be excluded from distance analysis | - |
| method | Transmission clustering method [`SNP` or `trans`] | trans |
| snp_threshold | SNP distance threshold for transmission clustering | 12 |
| transmission_threshold | Transmission threshold for transmission clustering | 15 |
| clock_rate | Clock rate for MTBC samples (SNPs/genome/year) | 1.0 |
| transmission_rate | The rate at which the estimated number of intermediate transmissions must be | 2.0 |
| coordinate | Using coordinate to reconstruct transmission network [`True` or `False`] | True |
| characteristics | Epidemiological characteristics for risk factor inference | - |
| sample_threads | Number of threads for each sample run | 4 |

## 🍾 License

The code is available under the [GNU GPLv3 license](https://choosealicense.com/licenses/gpl-3.0/). The text and data are availabe under the [CC-BY license](https://choosealicense.com/licenses/cc-by-4.0/).

## 📲 Questions and Issues

For contacting the developer and issue reports please go to [Issues](https://gitee.com/cwmda/transflow/issues).
