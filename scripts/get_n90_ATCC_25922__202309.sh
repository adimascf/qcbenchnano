#!/usr/bin/env bash
set -euo pipefail


output_csv=$1

echo "combo,depth,model,n90" > "$output_csv"

for file in /scratch/project/bug_seq_scratch/qc_bench/results/QC/downsampling/*/*/*/ATCC_25922__202309.*.rasusa.fastq; do
    
    # Extract metadata by parsing the parent directories
    model=$(basename $(dirname "$file"))
    depth=$(basename $(dirname $(dirname "$file")))
    combo=$(basename $(dirname $(dirname $(dirname "$file"))))
    
    # Calculate N90 and remove the comma (e.g., converting "1,868" to "1868")
    n90_val=$(seqkit stat -N 90 "$file" | awk '{print $NF}' | tail -n 1 | tr -d ',')
    
    echo "$combo,$depth,$model,$n90_val" >> "$output_csv"

done

