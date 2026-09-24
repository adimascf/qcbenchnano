import sys
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
from pathlib import Path

# Redirect all prints and errors to the Snakemake log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# --- Configuration ---
baseline_combo = "unprocessed-untrimmed" 
out_dir = Path(snakemake.output[0]).parent
out_dir.mkdir(parents=True, exist_ok=True)

# Define which metrics to test and whether a HIGHER score is better
# The script will only test metrics that actually exist in the loaded CSV
metric_properties = {
    "F1_SCORE": True,
    "PREC": True,
    "RECALL": True,
    "NGA50_norm": True,
    "Mismatches per 100kbp": False,
    "Indels per 100kbp": False
}

print(f"Loading data from {snakemake.input.csv}...")
df = pd.read_csv(snakemake.input.csv)

# Dynamically determine grouping columns based on the file type
if "VAR_TYPE" in df.columns:
    group_cols = ["VAR_TYPE", "model", "depth"]
    # Filter out 'ALL' and 'SV' for variant data
    df = df.query("VAR_TYPE not in ('ALL', 'SV')").copy()
    
else:
    # Assembly files do not have VAR_TYPE
    group_cols = ["model", "depth"]
    if "NGA50_norm" in df.columns:
        dataix = df.groupby(["combo", "depth", "sample", "model"])["NGA50_norm"].idxmax()
        df = df.loc[dataix].copy()

for metric, higher_is_better in metric_properties.items():
    if metric not in df.columns:
        continue # Skip metrics that aren't in this specific CSV
        
    print(f"Processing {metric}...")
    results = []
    groups = df.groupby(group_cols)
    
    for group_keys, group_data in groups:
        # Create a dictionary for the group keys to easily add them to the results
        group_dict = dict(zip(group_cols, group_keys if isinstance(group_keys, tuple) else (group_keys,)))
        
        baseline_data = group_data[group_data["combo"] == baseline_combo]
        if baseline_data.empty:
            continue
            
        for combo in group_data["combo"].unique():
            if combo == baseline_combo:
                continue 
                
            target_data = group_data[group_data["combo"] == combo]
            merged = pd.merge(
                baseline_data[["sample", metric]], 
                target_data[["sample", metric]], 
                on="sample", 
                suffixes=('_base', '_target')
            )
            
            if len(merged) < 2:
                p_value = np.nan
                median_diff = np.nan
                direction = "Unknown"
            else:
                diff = merged[f"{metric}_target"] - merged[f"{metric}_base"]
                median_diff = diff.median()
                
                # Dynamic direction logic based on the metric properties
                if median_diff == 0:
                    direction = "No Change"
                elif (median_diff > 0 and higher_is_better) or (median_diff < 0 and not higher_is_better):
                    direction = "Improved"
                else:
                    direction = "Worsened"
                
                if np.all(diff == 0):
                    p_value = 1.0 
                else:
                    try:
                        stat, p_value = wilcoxon(merged[f"{metric}_base"], merged[f"{metric}_target"])
                    except ValueError:
                        p_value = np.nan
                        
            result_row = {
                "combo": combo,
                **group_dict,
                "median_diff": median_diff,
                "direction": direction,
                "p_value": p_value
            }
            results.append(result_row)
            
    if not results:
        continue
        
    results_df = pd.DataFrame(results)
    
    # Multiple Testing Correction (Benjamini-Hochberg FDR)
    valid_mask = results_df['p_value'].notna()
    results_df['p_adj_BH'] = np.nan 
    
    if valid_mask.any():
        # Adjust p-values to control the False Discovery Rate
        _, p_adj, _, _ = multipletests(results_df.loc[valid_mask, 'p_value'], method='fdr_bh')
        results_df.loc[valid_mask, 'p_adj_BH'] = p_adj
        
    # Flag significant results based on the ADJUSTED p-value
    results_df['significant'] = results_df['p_adj_BH'] < 0.05
    
    # Sort for readability
    sort_cols = group_cols + ["p_adj_BH"]
    results_df = results_df.sort_values(by=sort_cols)
    
    # Format filename safely
    metric_safe_name = metric.replace(" ", "_").replace("/", "_").lower()
    output_file = out_dir / f"combo_{metric_safe_name}_significance.csv"
    results_df.to_csv(output_file, index=False)
    print(f"Saved statistical results to {output_file.name}")

print("Statistical testing complete.")
