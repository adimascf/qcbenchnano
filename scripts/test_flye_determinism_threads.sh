#!/usr/bin/env bash
set -euo pipefail

READS=$1
REF=$2
PLASMID=$3
ITERATIONS=${4:-4} # Defaulting to 3 runs per thread count to save time

THREAD_COUNTS=(1 2 4 8 16 32)

# Initialise the summary table
depth=$(echo "$READS" | awk -v FS="/" '{print $10}')
SUMMARY_FILE="flye_determinism_summary_${PLASMID}.${depth}.tsv"

echo -e "Threads\tRun\tPlasmid\tCoverage(%)\tMeanDepth" > "$SUMMARY_FILE"

echo "Testing Flye determinism for $PLASMID..."
echo "Reads: $READS"
echo "Threads to test: ${THREAD_COUNTS[*]}"
echo "Iterations per thread count: $ITERATIONS"
echo "------------------------------------------------------------"

# Outer loop: Iterate through the different thread counts
for THREADS in "${THREAD_COUNTS[@]}"; do
    
    # Inner loop: Run the set number of iterations for the current thread count
    for i in $(seq 1 "$ITERATIONS"); do
        echo "Starting Flye (Threads: $THREADS, Run: $i/$ITERATIONS)..."
        
        OUT_DIR="flye_t${THREADS}_iter_${i}"

        flye --debug --min-overlap 1000 --nano-hq "$READS" --out-dir "$OUT_DIR" --threads "$THREADS" 

        echo "  Assessing assembly..."
        minimap2 -a -x asm5 "$REF" "$OUT_DIR/assembly.fasta" 2>/dev/null | \
        samtools sort 2>/dev/null | \
        samtools coverage - > "$OUT_DIR/coverage.tsv"

        # Extract the target metrics (Column 1 is name, Col 6 is coverage, Col 7 is meandepth)
        METRICS=$(grep -w "$PLASMID" "$OUT_DIR/coverage.tsv" | awk '{print $6 "\t" $7}')

        # Append to the summary table with the current thread count
        echo -e "${THREADS}\t${i}\t${PLASMID}\t${METRICS}" >> "$SUMMARY_FILE"
    done
done

echo "------------------------------------------------------------"
echo "All runs complete. Results saved to: $SUMMARY_FILE"
cat "$SUMMARY_FILE"
