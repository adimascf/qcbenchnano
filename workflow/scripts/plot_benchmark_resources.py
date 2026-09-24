# Adapted from https://github.com/mbhall88/NanoVarBench/blob/main/workflow/scripts/benchmark_plot.py

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_rgba
from typing import List

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# colourblind-friendly palette from colour universal design (CUD)
named_colors = {
    "black": "#000000",
    "orange": "#e69f00",
    "skyblue": "#56b4e9",
    "bluish green": "#009e73",
    "yellow": "#f0e442",
    "blue": "#0072b2",
    "vermilion": "#d55e00",
    "reddish purple": "#cc79a7",
}
cud_palette = list(named_colors.values())
sns.set_theme(style="whitegrid")

def cud(n: int = len(cud_palette), start: int = 0) -> List[str]:
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]

# Gather all TSVs and tag them with their processing category
tsvs = []
if hasattr(snakemake.input, 'trimming_benchmark'):
    tsvs.extend([(Path(p), "Adapter Trimming") for p in snakemake.input.trimming_benchmark])
if hasattr(snakemake.input, 'quality_benchmark'):
    tsvs.extend([(Path(p), "Quality Processing") for p in snakemake.input.quality_benchmark])

# Parse Snakemake Benchmark Data
frames = []
for p, category in tsvs:
    df = pd.read_csv(p, sep="\t")
    
    sample = p.name.split('.')[0]
    model = p.parts[-2]
    combo_name = p.parts[-3] 
    
    base_tool = combo_name 
    if category == "Quality Processing" and "-" in combo_name:
        base_tool = combo_name.split("-")[0] 

    df["sample"] = sample
    df["model"] = model
    df["category"] = category
    df["combo_name"] = combo_name # if we find "filtlong-dorado"
    df["tool"] = base_tool # here we just exatract the quality processed tool, which is only "filtlong"
    
    time_col = "cpu_time" if "cpu_time" in df.columns else "s"
    df["time_s"] = df[time_col]
    
    frames.append(df)

df = pd.concat(frames)

# Save the raw aggregated benchmark data
df.to_csv(snakemake.output.csv, index=False)
print(f"Saved compiled metrics to {snakemake.output.csv}")

# Group the tools for the categories
trimming_tools = sorted(df.query("category == 'Adapter Trimming'")["tool"].unique())
quality_tools = sorted(df.query("category == 'Quality Processing'")["tool"].unique())

order = trimming_tools + quality_tools
palette = {tool: cud(len(order), start=1)[i] for i, tool in enumerate(order)}

# Set up Grid Parameters
categories = ["Adapter Trimming", "Quality Processing"]
metrics = ["max_rss", "time_s"] 
metric_titles = {"max_rss": "Maximum Memory Usage", "time_s": "Total CPU Time"}

# Create a 2x2 grid
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 9), dpi=300, sharey='row')
violin_alpha = 0.15
jitter = 0.2

for r, category in enumerate(categories):
    for c, metric in enumerate(metrics):
        ax = axes[r, c]
        
        # Subset data for this specific category row
        df_sub = df.query("category == @category")
        cat_tools = trimming_tools if category == "Adapter Trimming" else quality_tools
        
        kwargs = dict(
            data=df_sub, x=metric, y="tool", hue="tool",
            order=cat_tools, hue_order=cat_tools, palette=palette,
            orient="h", dodge=False, legend=False
        )

        sns.violinplot(**kwargs, cut=0, inner="quartile", ax=ax)
        for collection in ax.collections:
            if type(collection).__name__ == "PolyCollection":
                collection.set_facecolor(to_rgba(collection.get_facecolor(), alpha=violin_alpha))

        # Plot Stripplots
        sns.stripplot(
            **kwargs, 
            alpha=0.7, size=5, edgecolor="black", linewidth=0.5, 
            ax=ax, jitter=jitter
        )

        ax.set_xscale("log")
        
        # Grid Formatting
        ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)
        
        # Specific Formatting per Metric Column
        if metric == "max_rss": 
            ticks = [(100, "100MB"), (500, "500MB"), (1000, "1GB"), (2000, "2GB"), (4000, "4GB"), (8000, "8GB"), (16000, "16GB")]
        else: 
            ticks = [(1, "1s"), (10, "10s"), (60, "1m"), (600, "10m"), (3600, "1h"), (14400, "4h")]
        
        ax.set_xticks([t[0] for t in ticks])
        
        # X-Axis Formatting (Only show labels and titles on the bottom row to prevent clutter)
        if r == 1:
            ax.set_xticklabels([t[1] for t in ticks], fontsize=12)
            ax.set_xlabel("")
        else:
            ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
            ax.set_xlabel("")
            
        # Y-Axis Formatting (Only show tool names and category titles on the left column)
        if c == 0:
            ax.set_yticklabels(cat_tools, fontsize=12)
            ax.set_ylabel(category, fontsize=12, fontweight='bold', labelpad=15)
        else:
            ax.tick_params(axis='y', which='both', left=False, labelleft=False)
            ax.set_ylabel("")
        if r == 0:
            ax.set_title(metric_titles[metric], fontsize=13, fontweight='bold', pad=14)

# Force the y-axis labels in the leftmost column to align horizontally
fig.align_ylabels(axes[:, 0])

fig.tight_layout()
fig.savefig(snakemake.output.figure, bbox_inches='tight')
plt.close(fig)
print(f"Saved resource plot to {snakemake.output.figure}")
