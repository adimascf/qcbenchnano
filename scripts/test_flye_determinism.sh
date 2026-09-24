#!/usr/bin/env bash
set -euo pipefail

READS=$1
REF=$2
PLASMID=$3
ITERATIONS=${4:-10} # Defaults to 10 if a 4th parameter isn't provided
THREADS=16

# Initialise the summary table
depth=$(echo $READS | awk -v FS="/" '{print $10}')
SUMMARY_FILE="flye_determinism_summary_${PLASMID}.${depth}.tsv"
echo -e "Run\tPlasmid\tCoverage(%)\tMeanDepth" > "$SUMMARY_FILE"

echo "Testing Flye determinism for $PLASMID over $ITERATIONS runs..."
echo "Reads: $READS"
echo "------------------------------------------------------------"

# Start the loop
for i in $(seq 1 $ITERATIONS); do
    echo "Starting Flye run $i/$ITERATIONS..."
    OUT_DIR="flye_iter_$i"

    flye --debug --nano-hq "$READS" --out-dir "$OUT_DIR" --threads "$THREADS" 

    # Run the assessment pipeline
    echo "  Assessing assembly $i..."
    minimap2 -a -x asm5 "$REF" "$OUT_DIR/assembly.fasta" 2>/dev/null | \
    samtools sort 2>/dev/null | \
    samtools coverage - > "$OUT_DIR/coverage.tsv"

    # Extract the target metrics (Column 1 is name, Col 6 is coverage, Col 7 is meandepth)
    METRICS=$(grep -w "$PLASMID" "$OUT_DIR/coverage.tsv" | awk '{print $6 "\t" $7}')

    # Append to the summary table
    echo -e "$i\t$PLASMID\t$METRICS" >> "$SUMMARY_FILE"
done

echo "------------------------------------------------------------"
echo "All runs complete. Results saved to: $SUMMARY_FILE"
cat "$SUMMARY_FILE"
