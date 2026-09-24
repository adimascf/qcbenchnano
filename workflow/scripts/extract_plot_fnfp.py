import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List

# Redirect both stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Reading metrics from {snakemake.input.csv}...")

# colourblind-friendly palette from colour universal design (CUD)
named_colors = {
    "black": "#000000",
    "orange": "#e69f00",
    "skyblue": "#56b4e9",
    "vermilion": "#d55e00",
    "bluish green": "#009e73",
    "yellow": "#f0e442",
    "blue": "#0072b2",
    "reddish purple": "#cc79a7",
}
cud_palette = list(named_colors.values())

def cud(n: int = len(cud_palette), start: int = 0) -> List[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]

# Load data from the pre compiled csv
df = pd.read_csv(snakemake.input.csv)

# Group by relevant columns and calculate the absolute SUM of errors across all samples
df_grouped = df.groupby(['combo', 'model', 'depth', 'VAR_TYPE'])[['TRUTH_FN', 'QUERY_FP']].sum().reset_index()

# Save aggregated sums to the csv file
df_grouped.to_csv(snakemake.output.csv, index=False)
print(f"Saved aggregated FN/FP sums to {snakemake.output.csv}")


sns.set_theme(style="whitegrid")

metrics_to_plot = {
    "TRUTH_FN": snakemake.output.fn_plot,
    "QUERY_FP": snakemake.output.fp_plot
}

rows = ["SNP", "INDEL"]
models = ["sup", "hac"]

# Fetch the depths dynamically from the Snakemake config
config_depths = snakemake.config['depth']

# sort the depths numerically, and add 'x'
hue_order = [f"{d}x" for d in sorted([int(d) for d in config_depths], reverse=True)]

for metric_col, output_file in metrics_to_plot.items():
    
    # Calculate overall performance for sorting
    combo_totals = df_grouped.groupby('combo')[[metric_col]].sum().sum(axis=1)
    order = combo_totals.sort_values(ascending=True).index.tolist()
    
    fig, axes = plt.subplots(
        nrows=len(rows), 
        ncols=len(models), 
        figsize=(14, 10), 
        dpi=300, 
        sharex=True, 
        sharey=True
    )

    for r, vartype in enumerate(rows):
        for c, model in enumerate(models):
            ax = axes[r, c]
            df_sub = df_grouped.query("VAR_TYPE == @vartype and model == @model").copy()

            if not df_sub.empty:
                sns.barplot(
                    data=df_sub, x="combo", y=metric_col, hue="depth",
                    order=order, hue_order=hue_order,
                    palette=cud(len(hue_order), start=2),
                    ax=ax,
                    dodge=True,
                    edgecolor="black", # Adds a clean border to the bars
                    linewidth=0.5
                )

            # Clean formatting
            error_name = "False Negatives" if metric_col == "TRUTH_FN" else "False Positives"
            
            if c == 0:
                ax.set_ylabel(f"{vartype} Total {error_name}")
            else:
                ax.set_ylabel("")
                
            ax.set_title(f"{model} model" if r == 0 else "")
            ax.set_xlabel("")
            
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(order, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
            ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

            # Legend handling (only top right)
            if r == 0 and c == 1 and ax.get_legend() is not None:
                ax.legend(title="Depth", loc="upper left", bbox_to_anchor=(1.02, 1))
            elif ax.get_legend() is not None:
                ax.get_legend().remove()

    fig.tight_layout()
    fig.savefig(output_file, bbox_inches='tight')
    plt.close(fig)
    
    print(f"Saved {error_name} plot to {output_file}")
