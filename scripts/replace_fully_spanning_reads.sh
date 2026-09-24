#!/usr/bin/env bash
set -euo pipefail

IDS_FILE="highcov_paradox/BPH2947_100x_to_plasmid.seqkit-dorado.spanning_ids.txt"
# IDS_FILE="highcov_paradox/BPH2947_100x_to_plasmid.seqkit-dorado.all_mapped_ids.txt"

DORADO_FASTQ="/scratch/project/bug_seq_scratch/qc_bench/results/QC/downsampling/seqkit-dorado/100x/sup/BPH2947__202310.seqkit-dorado.rasusa.fastq"
RAW_FASTQ="data/fastqs/notrim/sup/BPH2947__202310.fastq"


DORADO_FILTERED="highcov_paradox/BPH2947_100x.dorado_nospanning.fastq"
DORADO_SPANNING="highcov_paradox/BPH2947_100x.dorado_spanning_only.fastq"
RAW_SPANNING="highcov_paradox/BPH2947_raw_spanning_only.fastq"
HYBRID_FASTQ="/scratch/project/bug_seq_scratch/qc_bench/results/QC/downsampling/seqkit-dorado/100x/sup/BPH2947_100x.hybrid_experiment.fastq"

# Removing fully spanning reads from the Dorado dataset..."
seqkit grep -v -f "$IDS_FILE" "$DORADO_FASTQ" -o "$DORADO_FILTERED"

# Extract the fully spanning reads as well
seqkit grep -f "$IDS_FILE" "$DORADO_FASTQ" -o "$DORADO_SPANNING"

# Extracting the exact fully spanning reads from the Raw dataset..."
seqkit grep -f "$IDS_FILE" "$RAW_FASTQ" -o "$RAW_SPANNING"

# Merging into hybrid FASTQ..."
# cat "$DORADO_FILTERED" "$RAW_SPANNING" > "$HYBRID_FASTQ"

