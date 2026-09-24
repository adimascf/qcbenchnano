import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Reading contamination summary from {snakemake.input.count}...")

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

# Read input directly from Snakemake
contamination_count = pd.read_csv(snakemake.input.csv)

aggregated_counts = contamination_count.groupby(
    ["combo", "model", "depth"], as_index=False
)["contamination_count"].sum()

sns.set_theme(style="whitegrid")
models = ["sup", "hac"]

# Fetch the depths dynamically from the Snakemake config
config_depths = snakemake.config['depth']

# sort the depths numerically, and add 'x'
hue_order = [f"{d}x" for d in sorted([int(d) for d in config_depths], reverse=True)]

combo_perf = aggregated_counts.groupby("combo")["contamination_count"].sum()
order_abs = combo_perf.sort_values(ascending=True).index.tolist()

fig, axes = plt.subplots(nrows=1, ncols=len(models), figsize=(14, 5.5), dpi=300, sharex=True, sharey=True)

for c, model in enumerate(models):
    ax = axes[c]
    df_sub = aggregated_counts.query("model == @model").copy()

    if not df_sub.empty:
        sns.stripplot(
            data=df_sub, x="combo", y="contamination_count", hue="depth", 
            order=order_abs, hue_order=hue_order,
            palette=cud(len(hue_order), start=2), 
            ax=ax, alpha=0.9, dodge=True, linewidth=0.5, edgecolor="black"
        )

    # Force standard numbers instead of scientific notation
    ax.yaxis.set_major_formatter(plt.ScalarFormatter())

    ax.set_ylabel("Total Contamination Count" if c == 0 else "")
    ax.set_title(f"Total Contaminants - {model.upper()} Model")
    ax.set_xlabel("")

    ax.set_xticks(range(len(order_abs)))
    ax.set_xticklabels(order_abs, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
    ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

    if c == 1 and ax.get_legend() is not None:
        # Extract handles and labels to recreate it outside the plot
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
    elif ax.get_legend() is not None:
        # Remove legend from the first plot
        ax.get_legend().remove()

fig.tight_layout()

fig.savefig(snakemake.output.figure, bbox_inches='tight')
plt.close(fig)

print(f"Saved contamination plot to {snakemake.output.figure}")
