# originally from Michael Hall https://github.com/mbhall88/NanoVarBench/blob/main/workflow/scripts/best_f1.py

import sys
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
from typing import List


# redirect both stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

# colourblind-friendly palette from colour universal design (CUD)
# https://jfly.uni-koeln.de/color/
# https://github.com/mbhall88/cud
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

sns.set_theme(style="whitegrid")

out_dir = Path(snakemake.output.figures[0]).parent

# Load all TSV files 
frames = []
for p in map(Path, snakemake.input.pr):
    df = pd.read_csv(p, sep="\t")
    # Path parsing: .../assess/call/<combo>/<depth>x/<model>/<sample>.precision-recall.tsv
    df["model"] = p.parts[-2]
    df["depth"] = p.parts[-3]
    df["combo"] = p.parts[-4]
    
    df["sample"] = p.name.split('.')[0]
    frames.append(df)

pr_df = pd.concat(frames)
pr_df.reset_index(inplace=True, drop=True)

metrics = ["F1_SCORE", "PREC", "RECALL"]
vartypes = ["SNP", "INDEL"]
models = ["sup", "hac"]
estimators = ["mean", "median"]
plot_types = ["point", "overlay"] # I want to generate pointplot and overlay (point + strip)

dataix = pr_df.groupby(["combo", "depth", "VAR_TYPE", "sample", "model"])["F1_SCORE"].idxmax()
data = pr_df.iloc[dataix].copy()

# Fetch the depths dynamically from the Snakemake config
config_depths = snakemake.config['depth']

# sort the depths numerically, and add 'x'
hue_order = [f"{d}x" for d in sorted([int(d) for d in config_depths], reverse=True)]

# Filter out 'ALL' and 'SV'
data = data.query("VAR_TYPE not in ('ALL', 'SV')")

# Generate the plots
for y in metrics:
    for est in estimators:

        # Calculate OVERALL performance to sort from best (highest) to worst (lowest)
        if est == "mean":
            combo_perf = data.groupby("combo")[y].mean()
            est_func = 'mean'
        else:
            combo_perf = data.groupby("combo")[y].median()
            est_func = 'median'

        order = combo_perf.sort_values(ascending=False).index.tolist()

        # Loop through both plot types for the current metric and estimator
        for p_type in plot_types:

            # Set up a 2x2 grid: rows are VAR_TYPE, cols are MODEL
            fig, axes = plt.subplots(
                nrows=len(vartypes),
                ncols=len(models),
                figsize=(15, 10), # Slightly wider to accommodate two columns comfortably
                dpi=300,
                sharex=True,
                sharey=False,
            )

            # Loop through rows (SNP, INDEL) and columns (sup, hac)
            for r, vartype in enumerate(vartypes):
                for c, model in enumerate(models):
                    ax = axes[r, c]
                    
                    # Filter data for this specific row/col combination
                    df = data.query("VAR_TYPE == @vartype and model == @model").copy()

                    cap = 0.99999
                    df.loc[:, y] = df[y].apply(lambda v: cap if v > cap else v)
                    yticks = [0.01, 0.1, 0.5, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
                    yticklabels = [f"{yval:.2%}" if yval < cap else "100%" for yval in yticks]

                    # Only plot if data exists for this combination
                    if not df.empty:
                        # Conditional rendering based on the current plot type in the loop
                        if p_type == "point":
                            sns.pointplot(
                                data=df, x="combo", y=y, hue="depth",
                                order=order, hue_order=hue_order,
                                palette=cud(len(hue_order), start=2),
                                ax=ax,
                                dodge=0.3,
                                errorbar=('pi', 100),
                                capsize=0.1,
                                err_kws={"linewidth": 1},
                                linewidth=1,
                                markersize=4,
                                estimator=est_func
                            )
                        elif p_type == "overlay":
                            sns.stripplot(
                                data=df, x="combo", y=y, hue="depth",
                                order=order, hue_order=hue_order,
                                palette=cud(len(hue_order), start=2),
                                ax=ax,
                                dodge=True,
                                alpha=0.4,
                                linewidth=0.5,
                                edgecolor="black",
                                zorder=1,
                                size=3
                            )
                            sns.pointplot(
                                data=df, x="combo", y=y, hue="depth",
                                order=order, hue_order=hue_order,
                                palette=cud(len(hue_order), start=2),
                                ax=ax,
                                dodge=0.3,
                                errorbar=('ci', 95),
                                capsize=0.1,
                                err_kws={"linewidth": 1},
                                linewidth=1,
                                markersize=4,
                                estimator=est_func,
                                zorder=2,
                                legend=False
                            )

                    # Format Y-axis (Shared across columns)
                    ax.set_yscale("logit", nonpositive="clip")
                    ax.set_yticks(yticks)
                    ax.set_yticklabels(yticklabels)

                    ylabel = {"F1_SCORE": "F1 score", "PREC": "Precision", "RECALL": "Recall"}[y]
                    
                    # Only put Y-axis labels on the leftmost column
                    if c == 0:
                        ax.set_ylabel(f"{vartype} {ylabel}")
                    else:
                        ax.set_ylabel("")

                    if r == 0:
                        ax.set_title(f"{model} model", fontsize=14)
                    else:
                        ax.set_title("")

                    ax.set_xlabel("")
                    ax.set_xticks(range(len(order)))
                    ax.set_xticklabels(order, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
                    ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

                    # Legend Handling: Show it only once, on the top-right plot
                    if r == 0 and c == 1:
                        if ax.get_legend() is not None:
                            # Move legend outside the plot to avoid overlapping data points
                            ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
                    else:
                        if ax.get_legend() is not None:
                            ax.get_legend().remove()

            fig.tight_layout()
            
            # Adjust the output file to include the plot type suffix
            metric_file = "f1" if y == "F1_SCORE" else "precision" if y == "PREC" else "recall"
            new_output_path = out_dir / f"combo_variant_{metric_file}_{est}_{p_type}plot.png"

            fig.savefig(new_output_path, bbox_inches='tight')
            plt.close(fig)

# Save the subsetted data to CSV
data.sort_values(by=["combo", "sample", "model", "depth"], inplace=True)
data.to_csv(snakemake.output.csv, index=False)
