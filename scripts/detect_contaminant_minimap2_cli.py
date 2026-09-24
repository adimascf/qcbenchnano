import subprocess
import os
import argparse
import sys
import glob

def run_minimap2(assembly_fa, barcodes_fa):
    """
    Run minimap with paramter specified for this task,
    treat barcode_fa as the target, and assembly_fa as the query.
    """
    cmd = [
        "minimap2", "-c", "-x", "sr",
        "-k", "9", "-w", "1", "-n", "1", "-m", "14",
        barcodes_fa, assembly_fa
    ]

    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    return process

def calc_overlap(aln_start, aln_end, reg_start, reg_end):
    """
    Calculate coverage or overlap per region in the barcode fasta sequence.
    """
    overlap_start = max(aln_start, reg_start)
    overlap_end = min(aln_end, reg_end)
    overlap_len = max(0, overlap_end - overlap_start)

    region_len = reg_end - reg_start
    if region_len == 0:
        return 0
    return overlap_len / region_len

def parse_filter(process, barcodes_fa, min_region_cov=0.80, min_identity=0.90):
    """
    Filtering strategies for the raw minimap2 PAF output.
    Returns a list of valid hit dictionaries.
    """
    filename = os.path.basename(barcodes_fa).lower()
    if "rapid" in filename:
        kit_type = "rapid"
    elif "native_left" in filename:
        kit_type = "native_left"
    elif "native_right" in filename:
        kit_type = "native_right"
    else:
        kit_type = "unknown"

    valid_hits = []

    for line in process.stdout:
        fields = line.strip().split('\t')
        if len(fields) < 12:
            continue

        contig_name = fields[0]
        contig_length = int(fields[1])
        contig_start = int(fields[2])
        contig_end = int(fields[3])

        strand = fields[4]

        barcode_name = fields[5]
        barcode_length = int(fields[6])
        barcode_start = int(fields[7])
        barcode_end = int(fields[8])

        matches = int(fields[9])
        aln_block_length = int(fields[10])
        mapping_quality = int(fields[11])

        identity = matches / aln_block_length if aln_block_length > 0 else 0

        dp_score = 0
        for tag in fields[12:]:
            if tag.startswith("AS:i:"):
                dp_score = int(tag.split(":")[-1])
                break

        # Calculate coverage using barcode_start and barcode_end
        if kit_type == "rapid":
            left_cov = calc_overlap(barcode_start, barcode_end, 0, 16)
            bar_cov = calc_overlap(barcode_start, barcode_end, 16, 40)
            right_cov = calc_overlap(barcode_start, barcode_end, 40, 90)
            breakdown = f"[Left: {left_cov:.0%} | Bar: {bar_cov:.0%} | Right: {right_cov:.0%}]"
            is_hit = max(left_cov, bar_cov, right_cov) >= min_region_cov
            cov_score = left_cov + bar_cov + right_cov

        elif kit_type == "native_left":
            adapt_cov = calc_overlap(barcode_start, barcode_end, 0, 36)
            bar_flank_cov = calc_overlap(barcode_start, barcode_end, 36, 76)
            breakdown = f"[Adapter: {adapt_cov:.0%} | Bar+Flank: {bar_flank_cov:.0%}]"
            is_hit = max(adapt_cov, bar_flank_cov) >= min_region_cov
            cov_score = adapt_cov + bar_flank_cov

        elif kit_type == "native_right":
            bar_flank_cov = calc_overlap(barcode_start, barcode_end, 0, 40)
            adapt_cov = calc_overlap(barcode_start, barcode_end, 40, 67)
            breakdown = f"[Adapter: {adapt_cov:.0%} | Bar+Flank: {bar_flank_cov:.0%}]"
            is_hit = max(adapt_cov, bar_flank_cov) >= min_region_cov
            cov_score = adapt_cov + bar_flank_cov

        else:
            all_cov = calc_overlap(barcode_start, barcode_end, 0, barcode_length)
            breakdown = f"[Total: {all_cov:.0%}]"
            is_hit = all_cov >= min_region_cov
            cov_score = all_cov

        # Positional logic uses the contig coordinates
        if contig_start <= 100:
            position = "start"
        elif contig_end >= (contig_length - 100):
            position = "end"
        else:
            position = "middle"

        if is_hit and identity >= min_identity:
            valid_hits.append({
                'barcode_name': barcode_name,
                'kit_type': kit_type,
                'contig_name': contig_name,
                'strand': strand,
                'c_start': contig_start,
                'c_end': contig_end,
                'identity': identity,
                'dp_score': dp_score,
                'mapq': mapping_quality,
                'position': position,
                'breakdown': breakdown,
                'cov_score': cov_score,
                'aln_len': contig_end - contig_start
            })

    process.wait()
    return valid_hits

def resolve_global_overlaps(all_hits):
    """
    Takes a combined list of hits from all kits, sorts them, 
    and resolves overlaps to find the single best matching barcode.
    """
    if not all_hits:
        return []

    # Sort globally by contig, then by start position
    all_hits.sort(key=lambda x: (x['contig_name'], x['c_start']))

    final_hits = []
    for hit in all_hits:
        if not final_hits:
            final_hits.append(hit)
            continue

        prev = final_hits[-1]

        is_same_contig = hit['contig_name'] == prev['contig_name']
        is_overlapping = hit['c_start'] < prev['c_end']

        if is_same_contig and is_overlapping:
            # using the dynamic programming score first from minimap2 result
            if hit['dp_score'] > prev['dp_score']:
                final_hits[-1] = hit
            elif hit['dp_score'] == prev['dp_score']:
                if hit['cov_score'] > prev['cov_score']:
                    final_hits[-1] = hit
                elif hit['cov_score'] == prev['cov_score']:
                    if hit['identity'] > prev['identity']:
                        final_hits[-1] = hit
        else:
            final_hits.append(hit)

    return final_hits

def main():
    parser = argparse.ArgumentParser(
        description="Find barcode contaminants in assemblies across multiple kits using minimap2."
    )

    parser.add_argument("-a", "--assembly", required=True, help="Path to the assembly FASTA file.")
    parser.add_argument("-d", "--barcode-dir", required=True, help="Directory containing the barcodes FASTA files.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output TSV file.")
    parser.add_argument("--min-cov", type=float, default=0.90, help="Minimum region coverage threshold (default: 0.90).")
    parser.add_argument("--min-id", type=float, default=0.90, help="Minimum identity threshold (default: 0.90).")

    args = parser.parse_args()

    if not os.path.isfile(args.assembly):
        sys.exit(f"Error: Assembly file '{args.assembly}' not found.")

    if not os.path.isdir(args.barcode_dir):
        sys.exit(f"Error: Barcode directory '{args.barcode_dir}' not found.")

    # Find all fasta/fa files in the directory
    barcode_files = glob.glob(os.path.join(args.barcode_dir, "*.fasta")) + glob.glob(os.path.join(args.barcode_dir, "*.fa"))

    if not barcode_files:
        sys.exit(f"Error: No .fasta or .fa files found in '{args.barcode_dir}'.")

    all_valid_hits = []

    for barcodes_fa in barcode_files:
        output = run_minimap2(args.assembly, barcodes_fa)

        kit_hits = parse_filter(
            output,
            barcodes_fa=barcodes_fa,
            min_region_cov=args.min_cov,
            min_identity=args.min_id
        )
        all_valid_hits.extend(kit_hits)

    final_resolved_hits = resolve_global_overlaps(all_valid_hits)

    # Write final output to TSV file
    print(f"Detected {len(final_resolved_hits)} contaminants.")
    print(f"Writing results to {args.output}...")

    with open(args.output, "w") as out_file:
        out_file.write("query\tkit\ttarget_contig\tstrand\tc_start\tc_end\tidentity\tdp_score\tmapq\tposition\tregion_breakdown\n")
        for hit in final_resolved_hits:
            id_str = f"{hit['identity']:.2f}"
            out_file.write(f"{hit['barcode_name']}\t{hit['kit_type']}\t{hit['contig_name']}\t{hit['strand']}\t{hit['c_start']}\t{hit['c_end']}\t{id_str}\t{hit['dp_score']}\t{hit['mapq']}\t{hit['position']}\t{hit['breakdown']}\n")

    print("Done!")

if __name__ == '__main__':
    main()
