import subprocess
import os
import sys
import glob
from pathlib import Path

# Redirect both stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

def run_minimap2(assembly_fa, barcodes_fa):
    """
    Run minimap with parameter specified for this task,
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

def load_circularity(info_file):
    """
    Parse flye assembly_info.txt andreturn a dictionary of contig circularity.
    """
    circularity = {}
    if not os.path.exists(info_file):
        return circularity

    with open(info_file, 'r') as file:
        for line in file:
            if line.startswith("#"): # skip the header
                continue
            cols = line.strip().split("\t")
            if len(cols) >= 4:
                # cols[0] is seq_name, cols[3] is circ. (Y/N)
                circularity[cols[0]] = cols[3]
    return circularity

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

        dist_closest_end = min(contig_start, contig_length - contig_end)

        if is_hit and identity >= min_identity:
            valid_hits.append({
                'barcode_name': barcode_name,
                'kit_type': kit_type,
                'contig_name': contig_name,
                'strand': strand,
                'c_start': contig_start,
                'c_end': contig_end,
                'dist_end': dist_closest_end,
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

    I define overalp where two barcodes (usually come from different files)
    aliged overlap between each other relative to contig in the assembly
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


# --- SNAKEMAKE EXECUTION BLOCK ---
assembly_fasta = snakemake.input.assembly
assembly_info = snakemake.input.info
contaminants_dir = snakemake.input.contaminants
output_tsv = snakemake.output.contaminants

min_cov = getattr(snakemake.params, 'min_cov', 0.90)
min_id = getattr(snakemake.params, 'min_id', 0.90)
kit_name = getattr(snakemake.params, 'kit', "unknown")

os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True) 
assembly_label = Path(assembly_fasta).name

print(f"Starting contaminant detection for: {assembly_label}")
print(f"Library Kit specified: {kit_name}")

# Map the kit name to logical flags
is_native = kit_name == "SQK-NBD114-96"
is_rapid = kit_name == "SQK-RBK114-96"

circularity = load_circularity(assembly_info)

# Find all fasta/fa files in the directory
barcode_files = glob.glob(os.path.join(contaminants_dir, "*.fasta")) + glob.glob(os.path.join(contaminants_dir, "*.fa"))

if not barcode_files:
    print(f"Error: No .fasta or .fa files found in '{contaminants_dir}'.")
    sys.exit(1)

all_valid_hits = []

for barcodes_fa in barcode_files:
    filename = os.path.basename(barcodes_fa).lower()
    
    # SKIP LOGIC BASED ON KIT
    if is_native and "rapid" in filename:
        print(f"Skipping {os.path.basename(barcodes_fa)} (Kit is Native)...")
        continue
    elif is_rapid and "native" in filename:
        print(f"Skipping {os.path.basename(barcodes_fa)} (Kit is Rapid)...")
        continue

    print(f"Running minimap2 against {os.path.basename(barcodes_fa)}...")
    output = run_minimap2(assembly_fasta, barcodes_fa)
    
    kit_hits = parse_filter(
        output,
        barcodes_fa=barcodes_fa,
        min_region_cov=min_cov,
        min_identity=min_id
    )
    all_valid_hits.extend(kit_hits)

final_resolved_hits = resolve_global_overlaps(all_valid_hits)

# Write final output to TSV file
print(f"Detected {len(final_resolved_hits)} contaminants.")
print(f"Writing results to {output_tsv}...")

with open(output_tsv, "w") as out_file:
    out_file.write("query\tkit\ttarget_contig\tis_circular\tstrand\tdist_end\tidentity\tdp_score\tmapq\tposition\tregion_breakdown\n")
    for hit in final_resolved_hits:
        id_str = f"{hit['identity']:.2f}"
        is_circ = circularity.get(hit['contig_name'], 'Not found')
        out_file.write(f"{hit['barcode_name']}\t{hit['kit_type']}\t{hit['contig_name']}\t{is_circ}\t{hit['strand']}\t{hit['dist_end']}\t{id_str}\t{hit['dp_score']}\t{hit['mapq']}\t{hit['position']}\t{hit['breakdown']}\n")

print("Done!")
