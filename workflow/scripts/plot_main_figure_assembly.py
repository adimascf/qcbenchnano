import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import seaborn as sns

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Reading data for assembly main figure...")

# File paths from snakemake
QUAST_CSV = snakemake.input.quast_csv
MISSED_CSV = snakemake.input.missed_csv
CONTAM_CSV = snakemake.input.contam_csv
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
    print(f"Generating main assembly figure for model: {model}")

    df_quast = pd.read_csv(QUAST_CSV).query("model == @model").copy()
    df_missed = pd.read_csv(MISSED_CSV).query("model == @model").copy()
    df_contam = pd.read_csv(CONTAM_CSV).query("model == @model").copy()

    # Sort depths dynamically using config or from data
    if hasattr(snakemake, 'config') and 'depth' in snakemake.config:
        config_depths = snakemake.config['depth']
        hue_order = [f"{str(d).rstrip('x')}x" for d in sorted([int(str(d).rstrip('x')) for d in config_depths], reverse=True)]
    else:
        unique_depths = df_quast["depth"].unique()
        depth_ints = [int(str(d).replace('x', '')) for d in unique_depths]
        hue_order = [f"{d}x" for d in sorted(depth_ints, reverse=True)]

    # Standardize depth format and filter for active depths
    df_quast["depth"] = df_quast["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df_missed["depth"] = df_missed["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df_contam["depth"] = df_contam["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")

    df_quast = df_quast[df_quast["depth"].isin(hue_order)].copy()
    df_missed = df_missed[df_missed["depth"].isin(hue_order)].copy()
    df_contam = df_contam[df_contam["depth"].isin(hue_order)].copy()
        
    palette = cud(len(hue_order), start=2)

    # Calculate Total Errors
    df_quast["Total_Errors"] = df_quast["Mismatches per 100kbp"] + df_quast["Indels per 100kbp"]
    df_quast["auNGA_score"] = np.maximum(0, 1 - np.abs(df_quast["auNGA_ratio"] - 1.0))
    
    # Aggregated Contamination
    df_contam_agg = df_contam.groupby(["combo", "depth"], as_index=False)["contamination_count"].sum()

    # 2. Setup Figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300)
    axes = axes.flatten()

    def plot_metric(ax, data, x, y, ylabel, title, order, estimator=np.mean, is_bar=False, is_strip_only=False, is_logit=False):
        if is_bar:
            sns.barplot(
                data=data, x=x, y=y, hue="depth",
                order=order, hue_order=hue_order,
                palette=palette, ax=ax, edgecolor="black", linewidth=0.5
            )
        else:
            sns.stripplot(
                data=data, x=x, y=y, hue="depth",
                order=order, hue_order=hue_order,
                palette=palette, ax=ax, alpha=0.4 if not is_strip_only else 0.8, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=4
            )
            if not is_strip_only:
                sns.pointplot(
                    data=data, x=x, y=y, hue="depth",
                    order=order, hue_order=hue_order,
                    palette=palette, ax=ax, dodge=0.3, errorbar=('ci', 95), capsize=0.1,
                    err_kws={'linewidth': 1}, linewidth=1, markersize=5, estimator=estimator, legend=False, zorder=2
                )

        ax.set_title(title, fontsize=14, pad=10)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.set_xlabel("")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
        ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)
        
        if is_logit:
            cap = 0.99999
            yticks = [0.5, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
            yticklabels = [f"{yval * 100:g}%" if yval < cap else "100%" for yval in yticks]
            ax.set_yscale("logit", nonpositive="clip")
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels)

        # Fix legend
        if ax.get_legend() is not None:
            ax.get_legend().remove()

    # --- Plot A: Assembly Errors ---
    order_errors = df_quast.groupby("combo")["Total_Errors"].mean().sort_values(ascending=True).index.tolist()
    plot_metric(
        axes[0], df_quast, "combo", "Total_Errors", 
        "Total Errors per 100kbp\n(Mismatches + Indels)", "A. Assembly Errors", 
        order_errors, np.mean
    )

    # --- Plot B: auNGA ---
    order_aunga = df_quast.groupby("combo")["auNGA_score"].mean().sort_values(ascending=False).index.tolist()
    cap = 0.99999
    df_quast['plot_metric_aunga'] = df_quast['auNGA_score'].apply(lambda v: cap if v >= cap else (0.00001 if v <= 0 else v))
    plot_metric(
        axes[1], df_quast, "combo", "plot_metric_aunga", 
        "auNGA Score (1 - |auNGA ratio - 1|)", "B. Assembly Contiguity", 
        order_aunga, np.mean, is_logit=True
    )

    # --- Plot C: Missed Contigs ---
    if "total_missed" not in df_missed.columns and "full_missed" in df_missed.columns:
        df_missed["total_missed"] = df_missed["full_missed"] + df_missed.get("partial_missed", 0)

    df_missed_agg = df_missed.groupby(["combo", "depth"], as_index=False)["total_missed"].sum()
    order_missed = df_missed_agg.groupby("combo")["total_missed"].sum().sort_values(ascending=True).index.tolist()
    plot_metric(
        axes[2], df_missed_agg, "combo", "total_missed", 
        "Total Missed Contigs", "C. Missed Contigs", 
        order_missed, np.mean, is_bar=True
    )
    axes[2].yaxis.set_major_locator(MaxNLocator(integer=True))

    # --- Plot D: Contamination ---
    order_contam = df_contam_agg.groupby("combo")["contamination_count"].sum().sort_values(ascending=True).index.tolist()
    plot_metric(
        axes[3], df_contam_agg, "combo", "contamination_count", 
        "Total Contaminants", "D. Contamination Count", 
        order_contam, np.mean, is_bar=True
    )
    axes[3].yaxis.set_major_locator(MaxNLocator(integer=True))

    # Add a single unified legend for Depth
    from matplotlib.patches import Patch
    custom_handles = [Patch(facecolor=palette[i], edgecolor='black', label=hue_order[i]) for i in range(len(hue_order))]
    fig.legend(handles=custom_handles, title="Depth", loc="upper center", bbox_to_anchor=(0.5, 1.03), ncol=len(hue_order), fontsize=12, title_fontsize=14)

    plt.tight_layout()
    fig.savefig(OUT_FIG, bbox_inches='tight')
    plt.close(fig)
    print(f"Main figure saved to {OUT_FIG}")

if __name__ == "__main__":
    main()
