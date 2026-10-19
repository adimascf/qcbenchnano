# QC Benchmark for Nanopore Reads

## Evaluating Nanopore read trimming and quality filtering in bacterial genomics

This repository contains the benchmarking pipeline and analysis code for our paper: **Evaluating Nanopore read trimming and quality filtering in bacterial genomics**.

The study systematically assesses the performance of Oxford Nanopore Technologies (ONT) adapter trimming and quality-processing tools, evaluating their impact on key downstream bacterial genomics applications such as *de novo* genome assembly and variant calling across 13 diverse bacterial isolates.

- [Citation](#citation)
- [Interactive Dashboard](#interactive-dashboard)
- [Data and Accessions](#data-and-accessions)
- [Data Preparation: Creating Single-Kit Datasets](#data-preparation-creating-single-kit-datasets)
- [Clair3 Model Setup](#clair3-model-setup)
- [Pipeline Prerequisites](#pipeline-prerequisites)
- [Running the Pipeline on HPC (Slurm)](#running-the-pipeline-on-hpc-slurm)
- [Repository Structure](#repository-structure)

---

## Citation

> Furqon, A. D. C., Wick, R. R., Roberts, L. W., & Hall, M. B. (2026). Evaluating Nanopore read trimming and quality filtering in bacterial genomics. *bioRxiv*, doi: [10.64898/2026.10.12.758701](https://doi.org/10.64898/2026.10.12.758701).

---

## Interactive Dashboard

An interactive dashboard built with Quarto and Shinylive allows users to customize scoring weights and explore tool rankings for *de novo* assembly:

**[https://adimascf.github.io/qcbenchnano/](https://adimascf.github.io/qcbenchnano/)**

---

## Data and Accessions

All accessions, Figshare DOIs, and reference assemblies are catalogued in [`config/accessions.csv`](config/accessions.csv).

The dataset comprises 13 bacterial isolates sequenced with Oxford Nanopore R10.4.1 (E8.2) flow cells. While the original public raw reads included runs with mixed sequencing kits, each sample was filtered to only include reads from a single designated barcoding kit:

| Sample | Species | Assigned Kit | BioSample | Figshare POD5 DOI | Reference Assembly |
|:---|:---|:---|:---|:---|:---|
| `ATCC_10708__202309` | *Salmonella enterica* | `SQK-RBK114-96` | SAMN38321309 | [10.26188/25521883](https://doi.org/10.26188/25521883) | CP149507–CP149508 |
| `ATCC_17802__202309` | *Vibrio parahaemolyticus* | `SQK-RBK114-96` | SAMN38321311 | [10.26188/25495063](https://doi.org/10.26188/25495063) | CP149505–CP149506 |
| `ATCC_19119__202309` | *Listeria ivanovii* | `SQK-RBK114-96` | SAMN38321312 | [10.26188/25495057](https://doi.org/10.26188/25495057) | CP149497 |
| `ATCC_25922__202309` | *Escherichia coli* | `SQK-RBK114-96` | SAMN38321313 | [10.26188/25521892](https://doi.org/10.26188/25521892) | CP149500–CP149504 |
| `ATCC_35897__202309` | *Listeria welshimeri* | `SQK-RBK114-96` | SAMN38321316 | [10.26188/25495081](https://doi.org/10.26188/25495081) | CP149496 |
| `ATCC_BAA-679__202309` | *Listeria monocytogenes* | `SQK-RBK114-96` | SAMN38321317 | [10.26188/25495069](https://doi.org/10.26188/25495069) | CP149495 |
| `ATCC_33560__202309` | *Campylobacter jejuni* | `SQK-NBD114-96` | SAMN38321314 | [10.26188/25495054](https://doi.org/10.26188/25495054) | CP149499 |
| `ATCC_35221__202309` | *Campylobacter lari* | `SQK-NBD114-96` | SAMN38321315 | [10.26188/25493905](https://doi.org/10.26188/25493905) | CP149498 |
| `AJ292__202310` | *Klebsiella variicola* | `SQK-NBD114-96` | SAMN40453079 | [10.26188/25495048](https://doi.org/10.26188/25495048) | CP149491 |
| `BPH2947__202310` | *Staphylococcus aureus* | `SQK-NBD114-96` | SAMN40453078 | [10.26188/25495075](https://doi.org/10.26188/25495075) | CP149492–CP149494 |
| `RDH275__202311` | *Streptococcus pyogenes* | `SQK-NBD114-96` | SAMN40453081 | [10.26188/25495072](https://doi.org/10.26188/25495072) | CP149486 |
| `MMC234__202311` | *Streptococcus dysgalactiae* | `SQK-NBD114-96` | SAMN40453082 | [10.26188/25495066](https://doi.org/10.26188/25495066) | CP149485 |
| `AMtb_1__202402` | *Mycobacterium tuberculosis* | `SQK-NBD114-96` | SAMN40453083 | [10.26188/25495045](https://doi.org/10.26188/25495045) | CP149484 |

*(Note: An additional sample from the original raw download, `KPC2__202310`, was excluded because its available read depth was <100×).*

---

## Data Preparation: Creating Single-Kit Datasets

In the raw Figshare archives, nine isolates contained a mixture of reads from both the native barcoding kit (`SQK-NBD114-96`) and the rapid barcoding kit (`SQK-RBK114-96`). To avoid confounding trimmer evaluations with mixed adapter signatures, each sample was filtered to only include reads from a single kit before basecalling (see the initial investigation in [Logbook Week 11](https://adimascf.github.io/biox7021/logbook/semester1/week11.html) and [Week 12](https://adimascf.github.io/biox7021/logbook/semester1/week12.html)).

The complete workflow to regenerate the single-kit FASTQ datasets from raw POD5 files is as follows:

### Step 1: Download raw POD5 archives
Download and unpack all raw POD5 tarballs with verified checksums:
```bash
bash scripts/download_pod5.sh /path/to/raw_pod5
```

### Step 2: Check kit composition per sample
Scan the directory of unpacked POD5 files to inspect the percentage of each sequencing kit based on read metadata:
```bash
python scripts/extract_kit_pd5.py /path/to/raw_pod5/<sample_dir>
```

### Step 3: Strictly filter POD5 files to pure single-kit files
Filter the POD5 directory to move only `.pod5` files where **100% of the reads** belong strictly to the target kit (ignoring mixed files and files from other kits):
```bash
python scripts/move_pod5.py \
    --source /path/to/raw_pod5/<sample_dir> \
    --kit SQK-RBK114-96 \
    --dest /path/to/clean_pod5/<sample_dir>
```
*(Use `SQK-RBK114-96` or `SQK-NBD114-96` matching the table above).*

### Step 4: Basecall with Dorado
Basecall the filtered single-kit POD5 files with Dorado (v5.2.0) for both HAC and SUP models. Run once with default adapter trimming and once with `--no-trim` (for raw/untrimmed controls):

```bash
# Dorado trimmed (default trimming enabled; kit name automatically passed)
bash scripts/basecall.sh /path/to/clean_pod5 sup /path/to/bam/sup_trimmed

# Dorado untrimmed (--no-trim)
bash scripts/basecall.sh /path/to/clean_pod5 sup /path/to/bam/sup_untrimmed --no-trim

# Repeat for HAC model:
bash scripts/basecall.sh /path/to/clean_pod5 hac /path/to/bam/hac_trimmed
bash scripts/basecall.sh /path/to/clean_pod5 hac /path/to/bam/hac_untrimmed --no-trim
```

### Step 5: Convert BAM to FASTQ
Convert the Dorado BAM files to FASTQ:
```bash
bash scripts/bam_to_fastq.sh /path/to/bam/<model_dir> /path/to/fastqs/<model>
```

---

## Clair3 Model Setup

Variant calling in the benchmark is performed using **Clair3** (v2.0.1) in an Apptainer/Singularity container (`docker://hkubal/clair3:v2.0.1`). 

Because our reads are basecalled with Dorado v5.2.0, Clair3 requires the converted PyTorch models matching Dorado v5.2.0 from HKU-BAL's Rerio catalog. You can download and set up both models directly using `wget`:

```bash
# Set your model destination directory
CLAIR3_MODELS_DIR="/path/to/data/clair3_models"

# 1. Download SUP model (r1041_e82_400bps_sup_v520)
mkdir -p "${CLAIR3_MODELS_DIR}/r1041_e82_400bps_sup_v520"
wget -P "${CLAIR3_MODELS_DIR}/r1041_e82_400bps_sup_v520" \
    https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_sup_v520/full_alignment.pt
wget -P "${CLAIR3_MODELS_DIR}/r1041_e82_400bps_sup_v520" \
    https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_sup_v520/pileup.pt

# 2. Download HAC model (r1041_e82_400bps_hac_v520)
mkdir -p "${CLAIR3_MODELS_DIR}/r1041_e82_400bps_hac_v520"
wget -P "${CLAIR3_MODELS_DIR}/r1041_e82_400bps_hac_v520" \
    https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_hac_v520/full_alignment.pt
wget -P "${CLAIR3_MODELS_DIR}/r1041_e82_400bps_hac_v520" \
    https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/r1041_e82_400bps_hac_v520/pileup.pt
```

Ensure the path to these models is set in `workflow/scripts/calling/clair3.sh` (or configured via parameters).

---

## Pipeline Prerequisites

The benchmarking workflow is implemented in [Snakemake](https://snakemake.readthedocs.io/) and uses software deployment via Conda and Apptainer/Singularity.

To run the pipeline on an HPC cluster, make sure you have:
- **Snakemake** (>= 9.0) with **peppy** installed:
  ```bash
  conda create -n snakemake -c conda-forge -c bioconda snakemake peppy
  conda activate snakemake
  ```
- **Apptainer** or **Singularity**
- **Conda** or **Mamba**
- A Slurm cluster environment (or Slurm executor plugin for Snakemake)

---

## Running the Pipeline on HPC (Slurm)

### 1. Configuration
Check and update the configuration files before launching:
- [`config/config.yaml`](config/config.yaml): Defines basecalling models (`sup`, `hac`), target read depths (`20`, `100`), trimming tools (`porechop_abi`, `dorado`, `barbell`), and quality filtering tools (`fastplong`, `filtlong`, `chopper`, `seqkit`).
- [`config/pep/project_config.yaml`](config/pep/project_config.yaml): Sets the source directories for `fastqs`, `reference`, `mutreference`, and `truth_vcf`.
- [`workflow/Snakefile`](workflow/Snakefile): Set `RESULTS` and `DATA` root paths.

### 2. Launch via Slurm script
The master Snakemake workflow is submitted to Slurm using [`scripts/slurm_script.sh`](scripts/slurm_script.sh). This script automatically mounts scratch paths, sets up Slurm resource defaults, enables Conda and Apptainer execution (`--sdm conda apptainer`), and monitors job status:

```bash
# Ensure log directory exists
mkdir -p logs

# Submit the master Snakemake process
bash scripts/slurm_script.sh
```

You can pass additional Snakemake arguments directly through the script, for example to perform a dry-run:
```bash
bash scripts/slurm_script.sh -n
```

---

## Repository Structure

```text
├── config/
│   ├── config.yaml               # Main pipeline settings (models, depths, tools)
│   ├── accessions.csv            # Sample metadata and data accessions
│   └── pep/                      # PEP schema and sample specifications
├── workflow/
│   ├── Snakefile                 # Main Snakemake workflow
│   ├── rules/                    # Modular rule definitions (trim, quality, call, assembly, etc.)
│   ├── envs/                     # Conda environment YAML files
│   └── scripts/                  # Downstream assessment and plotting scripts
├── scripts/                      # Pre-pipeline utilities:
│   ├── download_pod5.sh          # Download raw Figshare POD5 files
│   ├── extract_kit_pd5.py        # Check kit composition across reads
│   ├── move_pod5.py              # Filter strictly pure single-kit POD5s
│   ├── basecall.sh               # Dorado basecalling (HAC/SUP, trimmed & untrimmed)
│   ├── bam_to_fastq.sh           # BAM to FASTQ extraction
│   └── slurm_script.sh           # HPC Slurm Snakemake submission script
└── docs/                         # Additional project documentation
```
