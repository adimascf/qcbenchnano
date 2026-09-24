import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Reading data from {snakemake.input.csv}...")

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

# Normalise auNGA_ratio for accurate ranking/sorting and also for logit
df["auNGA_score"] = np.maximum(0, 1 - np.abs(df["auNGA_ratio"] - 1.0))

sns.set_theme(style="whitegrid")
models = ["sup", "hac"]

# Fetch the depths dynamically from the Snakemake config
config_depths = snakemake.config['depth']

# sort the depths numerically, and add 'x'
hue_order = [f"{d}x" for d in sorted([int(d) for d in config_depths], reverse=True)]

plot_types = ["pointplot", "overlay"]
estimators = ["mean", "median"] 
scales = ["linear", "logit"] 

# Extract output directory dynamically from Snakemake
out_dir = Path(snakemake.output.figures[0]).parent

for est in estimators:

    if est == "mean":
        combo_perf = df.groupby("combo")["auNGA_score"].mean()
        est_func = np.mean
    else:
        combo_perf = df.groupby("combo")["auNGA_score"].median()
        est_func = np.median

    order_abs = combo_perf.sort_values(ascending=False).index.tolist()

    for p_type in plot_types:
        for scale in scales: 
            fig, axes = plt.subplots(nrows=1, ncols=len(models), figsize=(14, 5.5), dpi=300, sharex=True, sharey=False)

            for c, model in enumerate(models):
                ax = axes[c]
                df_sub = df.query("model == @model").copy()

                if not df_sub.empty:
                    
                    # Decide what goes on the Y-axis based on the scale
                    if scale == "logit":
                        cap = 0.99999
                        # Cap at 0.99999 strictly to prevent logit from crashing
                        df_sub['plot_metric'] = df_sub['auNGA_ratio'].apply(lambda v: cap if v >= cap else (0.00001 if v <= 0 else v))
                    else:
                        # For linear, use the exact raw values so we can see values > 1
                        df_sub['plot_metric'] = df_sub['auNGA_ratio']

                    if p_type == "pointplot":
                        sns.pointplot(
                            data=df_sub, x="combo", y="plot_metric", hue="depth",
                            order=order_abs, hue_order=hue_order,
                            palette=cud(len(hue_order), start=2),
                            ax=ax, dodge=0.3, errorbar=('pi', 100), capsize=0.1,
                            err_kws={'linewidth': 1}, linewidth=1, markersize=4, estimator=est_func
                        )
                    elif p_type == "overlay":
                        sns.stripplot(
                            data=df_sub, x="combo", y="plot_metric", hue="depth",
                            order=order_abs, hue_order=hue_order,
                            palette=cud(len(hue_order), start=2),
                            ax=ax, alpha=0.4, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=3
                        )
                        sns.pointplot(
                            data=df_sub, x="combo", y="plot_metric", hue="depth",
                            order=order_abs, hue_order=hue_order,
                            palette=cud(len(hue_order), start=2),
                            ax=ax, dodge=0.3, errorbar=('ci', 95), capsize=0.1,
                            err_kws={'linewidth': 1}, linewidth=1, markersize=4, estimator=est_func, legend=False, zorder=2
                        )

                # Dynamic Axis Scaling and Zoom
                if scale == "logit":
                    ax.set_yscale("logit", nonpositive="clip")

                    # Custom ticks strictly for the zoomed logit region
                    yticks = [0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
                    yticklabels = [f"{yval:.2%}" if yval < cap else "100%" for yval in yticks]
                    ax.set_yticks(yticks)
                    ax.set_yticklabels(yticklabels)

                    
                    ax.set_ylabel(f"[auNGA ratio (Distance from 1) | {est.title()}]" if c == 0 else "")
                    ax.set_title(f"auNGA ratio - {model} model")
                    ax.set_xlabel("")
                elif scale == "linear":
                    ax.set_yscale("linear")
                    # ax.set_ylim(bottom=-0.05) 
                    ax.set_ylim(0.95, 1.002)
                    ax.set_ylabel(f"[auNGA ratio (Raw) | {est.title()}]" if c == 0 else "")
                    ax.set_title(f"auNGA ratio - {model} model")
                    ax.set_xlabel("")

                ax.set_xticks(range(len(order_abs)))
                ax.set_xticklabels(order_abs, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
                ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

                # Legend handling
                if c == 1 and ax.get_legend() is not None:
                    ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
                elif ax.get_legend() is not None:
                    ax.get_legend().remove()

            fig.tight_layout()

            # mutate the output files
            new_output_path = out_dir / f"combo_assembly_aunga_normalised_{scale}_{p_type}_{est}.png"
            fig.savefig(new_output_path, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved {new_output_path.name}")
