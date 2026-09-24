import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Redirect stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Summarising assembly performance metrics...")

MASTER_CSV = snakemake.input.master_csv
OUTPUT_CSV = snakemake.output.csv

# 1. Load compiled master data
df = pd.read_csv(MASTER_CSV)

# 2. Filter for active model
target_model = getattr(snakemake.wildcards, 'model', 'sup')
print(f"Analysing assembly summary metrics for model: {target_model}")
df = df.query("model == @target_model").copy()

# 3. Handle and sort sampling depths dynamically (e.g. 100x and 20x)
if hasattr(snakemake, 'config') and 'depth' in snakemake.config:
    config_depths = snakemake.config['depth']
    depth_order = [f"{str(d).rstrip('x')}x" for d in sorted([int(str(d).rstrip('x')) for d in config_depths], reverse=True)]
else:
    unique_depths = df["depth"].unique()
    depth_ints = [int(str(d).replace('x', '')) for d in unique_depths]
    depth_order = [f"{d}x" for d in sorted(depth_ints, reverse=True)]

# Standardise depth format and filter active depths
df["depth"] = df["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
df = df[df["depth"].isin(depth_order)].copy()

# 4. Calculate per-sample composite error and contiguity metrics
df["total_errors_per_100kbp"] = df["Mismatches per 100kbp"] + df["Indels per 100kbp"]

aunga_col = "auNGA_ratio" if "auNGA_ratio" in df.columns else ("auNGA_norm" if "auNGA_norm" in df.columns else "auNGA")
df["aunga_score"] = np.maximum(0.0, 1.0 - np.abs(df[aunga_col] - 1.0))

nga50_norm_col = "NGA50_norm" if "NGA50_norm" in df.columns else "nga50_norm"
missed_full_col = "full_missed" if "full_missed" in df.columns else "full_missed_count"
missed_partial_col = "partial_missed" if "partial_missed" in df.columns else "partial_missed_count"
missed_total_col = "total_missed" if "total_missed" in df.columns else "total_missed_count"
contam_col = "contamination_count" if "contamination_count" in df.columns else "contaminant_count"

# 5. Summarise metrics across samples: means for errors & contiguity, totals for contig loss & contamination
summarised_df = df.groupby(["combo", "model", "depth"], as_index=False, observed=True).agg(
    sample_count=("sample", "count"),
    total_errors_per_100kbp_mean=("total_errors_per_100kbp", "mean"),
    mismatches_per_100kbp_mean=("Mismatches per 100kbp", "mean"),
    indels_per_100kbp_mean=("Indels per 100kbp", "mean"),
    aunga_score_mean=("aunga_score", "mean"),
    aunga_ratio_mean=(aunga_col, "mean"),
    nga50_normalised_mean=(nga50_norm_col, "mean"),
    full_missed_total=(missed_full_col, "sum"),
    partial_missed_total=(missed_partial_col, "sum"),
    total_missed_total=(missed_total_col, "sum"),
    contamination_count_total=(contam_col, "sum")
)

# 6. Apply clean ordering and type formatting
summarised_df["depth"] = pd.Categorical(summarised_df["depth"], categories=depth_order, ordered=True)

# Separate untrimmed / unprocessed to the end if present, otherwise sort alphabetically
unique_combos = sorted(summarised_df["combo"].unique().tolist())
combo_order = {c: i for i, c in enumerate(c for c in unique_combos if not c.startswith("unprocessed"))}
max_idx = len(combo_order)
for c in unique_combos:
    if c not in combo_order:
        combo_order[c] = max_idx
        max_idx += 1

summarised_df["combo_rank"] = summarised_df["combo"].map(combo_order)
summarised_df = summarised_df.sort_values(by=["combo_rank", "depth"]).drop(columns=["combo_rank"]).reset_index(drop=True)

# Ensure integer types for count metrics
integer_columns = ["sample_count", "full_missed_total", "partial_missed_total", "total_missed_total", "contamination_count_total"]
for col in integer_columns:
    if col in summarised_df.columns:
        summarised_df[col] = summarised_df[col].astype(int)

# Round floating point averages cleanly
float_columns = [
    "total_errors_per_100kbp_mean", "mismatches_per_100kbp_mean", "indels_per_100kbp_mean",
    "aunga_score_mean", "aunga_ratio_mean", "nga50_normalised_mean"
]
for col in float_columns:
    if col in summarised_df.columns:
        summarised_df[col] = summarised_df[col].round(6)

# 7. Write summarised results to CSV
Path(OUTPUT_CSV).parent.mkdir(parents=True, exist_ok=True)
summarised_df.to_csv(OUTPUT_CSV, index=False)
print(f"Successfully saved summarised assembly metrics table ({len(summarised_df)} rows) to: {OUTPUT_CSV}")
