import sys
import os
import pandas as pd
from pathlib import Path

# redirect both stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

data_frames = []
summary_records = []

# The exact column names written by the new minimap2 script
col_names = [
    'query', 'kit', 'target_contig', 'is_circular', 'strand', 
    'dist_end', 'identity', 'dp_score', 'mapq', 
    'position', 'region_breakdown'
]

print(f"Aggregating {len(snakemake.input.contam_report)} contamination files...")

for sample_report in snakemake.input.contam_report:
    p = Path(sample_report)
    
    # Path: .../contaminant/<combo>/<depth>/<model>/<sample>.file.tsv
    model = p.parts[-2]
    depth = p.parts[-3]
    combo = p.parts[-4]

    sample = p.name.split('.')[0]

    # Handle potentially completely empty files just in case, 
    # but primarily rely on pandas to read the header
    if os.path.getsize(p) == 0:
        df = pd.DataFrame(columns=col_names)
    else:
        # Read the TSV. It automatically picks up the header written by minimap2.
        df = pd.read_csv(p, sep='\t')

    # If the file only had a header and no data rows, length will be 0
    if len(df) == 0:
        summary_records.append({
            'combo': combo, 'depth': depth, 'sample': sample, 'model': model, 'contamination_count': 0
        })
    else:
        df['combo'] = combo
        df['depth'] = depth
        df['sample'] = sample
        df['model'] = model
        
        data_frames.append(df)
        
        # Record the actual number of contaminant hits
        summary_records.append({
            'combo': combo, 'depth': depth, 'sample': sample, 'model': model, 'contamination_count': len(df)
        })

if data_frames:
    combined_df = pd.concat(data_frames, ignore_index=True)
    # Reorder columns to put our pipeline metadata first
    final_cols = ['combo', 'depth', 'sample', 'model'] + col_names
    combined_df = combined_df[final_cols]
    
    combined_df.to_csv(snakemake.output.details, index=False)
    print(f"Saved detailed results to: {snakemake.output.details}")

else:
    # no contamination is found anywhere across all files
    pd.DataFrame(columns=['combo', 'depth', 'sample', 'model'] + col_names).to_csv(snakemake.output.details, index=False)
    print("No detailed contaminations found. Empty template created.")


if summary_records:
    summary_df = pd.DataFrame(summary_records)
    
    # apply custom sorting logic
    depth_order = {"100x": 1, "50x": 2, "20x": 3}
    summary_df['depth_sort'] = summary_df['depth'].map(depth_order)

    model_order = {"sup": 1, "hac": 2}
    summary_df['model_sort'] = summary_df['model'].map(model_order).fillna(3)
    
    trimmer_order = {t: i for i, t in enumerate(sorted(summary_df['combo'].unique())) if t != 'untrimmed'}
    trimmer_order['untrimmed'] = 999 # force to the end
    summary_df['trimmer_sort'] = summary_df['combo'].map(trimmer_order)
    
    # Sort the summary_df: sample -> model -> depth -> combo
    summary_df = summary_df.sort_values(
        by=['sample', 'model_sort', 'depth_sort', 'trimmer_sort']
    ).drop(columns=['depth_sort', 'model_sort', 'trimmer_sort'])
    
    summary_df.reset_index(drop=True, inplace=True)
    
    # Save as CSV
    summary_df.to_csv(snakemake.output.summary, index=False)
    print(f"Saved summary counts to: {snakemake.output.summary}")

print("Aggregation complete.")
