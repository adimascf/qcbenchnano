#!/usr/bin/env bash
set -euo pipefail

RAW_FASTQ="highcov_paradox/BPH2947_raw_spanning_only.fastq"
DORADO_FASTQ="highcov_paradox/BPH2947_100x.dorado_spanning_only.fastq"
OUT_TSV="highcov_paradox/trimming_comparison.tsv"

# Extracting lengths from raw and Dorado FASTQ files..."
seqkit fx2tab -l -n -i "$RAW_FASTQ" | sort -k1,1 > raw_len.tmp
seqkit fx2tab -l -n -i "$DORADO_FASTQ" | sort -k1,1 > dorado_len.tmp

# Calculating trimmed amounts..."

echo -e "ReadID\tRawLength\tTrimmedLength\tTrimmedAmount" > "$OUT_TSV"

# Join on the Read ID (Column 1), then subtract Dorado length from Raw length
join -1 1 -2 1 raw_len.tmp dorado_len.tmp | \
    awk -v OFS="\t" '{print $1, $2, $3, $2 - $3}' >> "$OUT_TSV"

rm raw_len.tmp dorado_len.tmp

