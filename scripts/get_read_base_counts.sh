#!/usr/bin/env bash
set -euo pipefail

# Use the directory provided as an argument, or default to the current directory
DIR="${1:-.}"

seqkit stats -T "$DIR"/*.f*q* | \
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 {print "File_Name", "Read_Count", "Base_Count"} NR>1 {print $1, $4, $5}' | \
column -t
