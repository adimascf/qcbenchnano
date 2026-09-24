#!/usr/bin/env bash
set -euo pipefail

REF_PLASMID="data/reference/BPH2947__202310.plasmid_2.fna"

READS_20X="/scratch/project/bug_seq_scratch/qc_bench/results/QC/downsampling/seqkit-dorado/20x/sup/BPH2947__202310.seqkit-dorado.rasusa.fastq"
READS_50X="/scratch/project/bug_seq_scratch/qc_bench/results/QC/downsampling/seqkit-dorado/50x/sup/BPH2947__202310.seqkit-dorado.rasusa.fastq"
READS_100X="/scratch/project/bug_seq_scratch/qc_bench/results/QC/downsampling/seqkit-dorado/100x/sup/BPH2947__202310.seqkit-dorado.rasusa.fastq"

mkdir -p highcov_paradox

OUT_20X="highcov_paradox/BPH2947_20x_to_plasmid.seqkit-dorado.bam"
OUT_50X="highcov_paradox/BPH2947_50x_to_plasmid.seqkit-dorado.bam"
OUT_100X="highcov_paradox/BPH2947_100x_to_plasmid.seqkit-dorado.bam"

PAF_20X="highcov_paradox/BPH2947_20x_to_plasmid.seqkit-dorado.paf"
PAF_50X="highcov_paradox/BPH2947_50x_to_plasmid.seqkit-dorado.paf"
PAF_100X="highcov_paradox/BPH2947_100x_to_plasmid.seqkit-dorado.paf"

THREADS=4

# Mapping 20x dataset
minimap2 -ax map-ont -t "$THREADS" "$REF_PLASMID" "$READS_20X" | \
    samtools sort -@ "$THREADS" -o "$OUT_20X" -
samtools index "$OUT_20X"
minimap2 --secondary=no -x map-ont -t "$THREADS" "$REF_PLASMID" "$READS_20X" > "$PAF_20X"

# Mapping 50x dataset 
minimap2 -ax map-ont -t "$THREADS" "$REF_PLASMID" "$READS_50X" | \
    samtools sort -@ "$THREADS" -o "$OUT_50X" -
samtools index "$OUT_50X"
minimap2 --secondary=no -x map-ont -t "$THREADS" "$REF_PLASMID" "$READS_50X" > "$PAF_50X"

# Mapping 100x dataset
minimap2 -ax map-ont -t "$THREADS" "$REF_PLASMID" "$READS_100X" | \
    samtools sort -@ "$THREADS" -o "$OUT_100X" -
samtools index "$OUT_100X"
minimap2 --secondary=no -x map-ont -t "$THREADS" "$REF_PLASMID" "$READS_100X" > "$PAF_100X"


PLASMID_LEN=3011

echo ""
echo "=== Read Length Statistics (Based on Total Read Length) ==="

for PAF in "$PAF_20X" "$PAF_50X" "$PAF_100X"; do
    echo "Dataset: $PAF"
    
    # Define output filenames based on the PAF filename
    OUT_IDS="${PAF%.paf}.spanning_ids.txt"
    ALL_IDS="${PAF%.paf}.all_mapped_ids.txt"
    BAM="${PAF%.paf}.bam"
    SPAN_BAM="${PAF%.paf}.spanning.bam"

    # Grab the total read length ($2) for each unique read ID ($1)
    awk -v plen="$PLASMID_LEN" -v out_ids="$OUT_IDS" -v all_ids="$ALL_IDS" '
        {
            # A read might have multiple alignments, but its total length ($2) is the same.
            read_len[$1] = $2
        }
        END {
            spanning = 0
            inner = 0
            total = 0
            for (read in read_len) {
                total++
                
                # Save ALL mapped read IDs to the new output file
                print read > all_ids
                
                if (read_len[read] >= plen) {
                    spanning++
                    # Save the fully spanning read ID to the output file
                    print read > out_ids
                } else {
                    inner++
                }
            }
            print "  Total unique mapped reads: " total
            print "  Fully span reads (total read length >= " plen "bp): " spanning
            print "  Inner reads (total read length < " plen "bp):    " inner
            print "  -> All mapped read IDs saved to: " all_ids
            print "  -> Spanning read IDs saved to: " out_ids
        }' "$PAF"
        
    # Extract only the spanning reads into a new BAM file and index it for IGV
    samtools view -N "$OUT_IDS" -o "$SPAN_BAM" "$BAM"
    samtools index "$SPAN_BAM"
    
    echo "  -> Spanning BAM saved to: $SPAN_BAM"
    echo ""
done
