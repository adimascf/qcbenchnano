import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Reading data for variant main figure...")

CSV_IN = snakemake.input.csv
OUT_FIG = snakemake.output.figure

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

def cud(n: int = len(cud_palette), start: int = 0):
    remainder = cud_palette[:start]
    palette = cud_palette[start:] + remainder
    return palette[:n]

def main():
    sns.set_theme(style="whitegrid")
    
    # 1. Load data
    model = getattr(snakemake.wildcards, 'model', 'sup')
    print(f"Generating main variant calling figure for model: {model}")

    df = pd.read_csv(CSV_IN)
    # Filter for target model and SNP/INDEL
    df = df.query("model == @model and VAR_TYPE in ('SNP', 'INDEL')").copy()

    # Sort depths dynamically
    if hasattr(snakemake, 'config') and 'depth' in snakemake.config:
        config_depths = snakemake.config['depth']
        hue_order = [f"{str(d).rstrip('x')}x" for d in sorted([int(str(d).rstrip('x')) for d in config_depths], reverse=True)]
    else:
        unique_depths = df["depth"].unique()
        depth_ints = [int(str(d).replace('x', '')) for d in unique_depths]
        hue_order = [f"{d}x" for d in sorted(depth_ints, reverse=True)]

    # Standardize depth format and filter for active depths
    df["depth"] = df["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df = df[df["depth"].isin(hue_order)].copy()
        
    palette = cud(len(hue_order), start=2)

    # 2. Setup Figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=300)
    axes = axes.flatten()

    def plot_metric(ax, data, metric, title, order):
        cap = 0.99999
        data_plot = data.copy()
        data_plot.loc[:, metric] = data_plot[metric].apply(lambda v: cap if v > cap else v)

        sns.stripplot(
            data=data_plot, x="combo", y=metric, hue="depth",
            order=order, hue_order=hue_order,
            palette=palette, ax=ax, alpha=0.4, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=4
        )
        sns.pointplot(
            data=data_plot, x="combo", y=metric, hue="depth",
            order=order, hue_order=hue_order,
            palette=palette, ax=ax, dodge=0.3, errorbar=('ci', 95), capsize=0.1,
            err_kws={'linewidth': 1}, linewidth=1, markersize=5, estimator=np.mean, legend=False, zorder=2
        )

        ax.set_title(title, fontsize=14, pad=10)
        ylabel = {"F1_SCORE": "F1 score", "PREC": "Precision", "RECALL": "Recall"}.get(metric, metric)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_xlabel("")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
        ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

        # Logit scale with 50% bottom cap
        yticks = [0.5, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
        yticklabels = [f"{yval:.2%}" if yval < cap else "100%" for yval in yticks]
        
        ax.set_yscale("logit", nonpositive="clip")
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(bottom=0.5)

        # Fix legend
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    # Get data subsets
    df_snp = df.query("VAR_TYPE == 'SNP'")
    df_indel = df.query("VAR_TYPE == 'INDEL'")

    # Calculate sort order based on mean F1_SCORE
    order_snp = df_snp.groupby("combo")["F1_SCORE"].mean().sort_values(ascending=False).index.tolist()
    order_indel = df_indel.groupby("combo")["F1_SCORE"].mean().sort_values(ascending=False).index.tolist()

    # --- Plot A: SNP F1_SCORE ---
    plot_metric(axes[0], df_snp, "F1_SCORE", "A. SNP F1 Score", order_snp)

    # --- Plot B: INDEL F1_SCORE ---
    plot_metric(axes[1], df_indel, "F1_SCORE", "B. INDEL F1 Score", order_indel)

    # Add a single unified legend for Depth
    from matplotlib.patches import Patch
    custom_handles = [Patch(facecolor=palette[i], edgecolor='black', label=hue_order[i]) for i in range(len(hue_order))]
    fig.legend(handles=custom_handles, title="Depth", loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=len(hue_order), fontsize=12, title_fontsize=14)

    plt.tight_layout()
    fig.savefig(OUT_FIG, bbox_inches='tight')
    plt.close(fig)
    print(f"Main figure saved to {OUT_FIG}")

if __name__ == "__main__":
    main()
