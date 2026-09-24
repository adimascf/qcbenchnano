import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Reading compiled metrics from {snakemake.input.csv}...")

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

df = pd.read_csv(snakemake.input.csv)

# Plotting preparation
sns.set_theme(style="whitegrid")
models = ["sup", "hac"]

# Fetch the depths dynamically from the Snakemake config
config_depths = snakemake.config['depth']

# sort the depths numerically, and add 'x'
hue_order = [f"{d}x" for d in sorted([int(d) for d in config_depths], reverse=True)]

plot_types = ["barplot", "overlay", "pointplot"]

# Calculate overall performance by summing both error types, then finding the global mean
df["Total_Errors"] = df["Mismatches per 100kbp"] + df["Indels per 100kbp"]
combo_perf = df.groupby("combo")["Total_Errors"].mean()
order_abs = combo_perf.sort_values(ascending=True).index.tolist()

metrics_abs = ["Mismatches per 100kbp", "Indels per 100kbp"]
out_dir = Path(snakemake.output.figures[0]).parent

for p_type in plot_types:
    fig, axes = plt.subplots(nrows=len(metrics_abs), ncols=len(models), figsize=(14, 10), dpi=300, sharex=True, sharey='row')

    for r, metric in enumerate(metrics_abs):
        for c, model in enumerate(models):
            ax = axes[r, c]
            df_sub = df.query("model == @model").copy()

            if not df_sub.empty:
                if p_type == "barplot":
                    sns.barplot(
                        data=df_sub, x="combo", y=metric, hue="depth", 
                        order=order_abs, hue_order=hue_order, 
                        palette=cud(len(hue_order), start=2),
                        ax=ax, dodge=True, edgecolor="black", linewidth=0.5,
                        estimator='mean' # bar height represents the mean
                    )
                elif p_type == "pointplot":
                    sns.pointplot(
                        data=df_sub, x="combo", y=metric, hue="depth", 
                        order=order_abs, hue_order=hue_order, 
                        palette=cud(len(hue_order), start=2),
                        ax=ax, dodge=0.3, errorbar=('pi', 100), capsize=0.1, 
                        err_kws={'linewidth': 1}, linewidth=1, markersize=4, estimator='mean' # point height represents the mean
                    ) 
                elif p_type == "overlay":
                    sns.stripplot(
                        data=df_sub, x="combo", y=metric, hue="depth", 
                        order=order_abs, hue_order=hue_order, 
                        palette=cud(len(hue_order), start=2),
                        ax=ax, alpha=0.4, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=3
                    )
                    sns.pointplot(
                        data=df_sub, x="combo", y=metric, hue="depth", 
                        order=order_abs, hue_order=hue_order, 
                        palette=cud(len(hue_order), start=2),
                        ax=ax, dodge=0.3, errorbar=('ci', 95), capsize=0.1, 
                        err_kws={'linewidth': 1}, linewidth=1, markersize=4, estimator='mean', legend=False, zorder=2
                    )

            # Formatting, making the plot prettier
            ax.set_ylabel(metric if c == 0 else "")
            ax.set_title(f"{model} model" if r == 0 else "")
            ax.set_xlabel("")
            ax.set_xticks(range(len(order_abs)))
            ax.set_xticklabels(order_abs, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
            ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

            # Legend handling
            if r == 0 and c == 1:
                if ax.get_legend() is not None:
                    ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
            else:
                if ax.get_legend() is not None:
                    ax.get_legend().remove()
    fig.tight_layout()
    
    # Save dynamically based on the current plot type in the loop
    new_output_path = out_dir / f"combo_assembly_errors_per_100kbp_{p_type}.png"
    fig.savefig(new_output_path, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {p_type} to {new_output_path}")

print("Finished plotting assembly metrics.")
