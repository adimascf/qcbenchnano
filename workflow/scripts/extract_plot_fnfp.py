import sys
from pathlib import Path
from typing import List
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns

# Detect execution context (Snakemake vs Standalone)
if "snakemake" in locals():
    # Redirect stdout and stderr to Snakemake log file if specified
    if hasattr(snakemake, "log") and snakemake.log:
        sys.stderr = sys.stdout = open(snakemake.log[0], "w")

    print(f"Reading metrics from {snakemake.input.csv}...")
    csv_in = Path(snakemake.input.csv)
    csv_out = Path(snakemake.output.csv)

    # Models & Depths dynamically from Snakemake config
    models = snakemake.config.get("model", ["sup", "hac"])
    config_depths = snakemake.config.get("depth", [100, 20])
    hue_order = [f"{str(d).rstrip('x')}x" for d in sorted([int(str(d).rstrip('x')) for d in config_depths], reverse=True)]

    # Map model to output figure paths
    model_plot_map = {}
    if hasattr(snakemake.output, "plots"):
        for p in snakemake.output.plots:
            p_path = Path(p)
            for m in models:
                if f"combo_call_fnfp_{m}.png" in p_path.name or f"_{m}." in p_path.name:
                    model_plot_map[m] = p_path
    for m in models:
        if m not in model_plot_map:
            if hasattr(snakemake.output, f"plot_{m}"):
                model_plot_map[m] = Path(getattr(snakemake.output, f"plot_{m}"))
            else:
                model_plot_map[m] = Path(f"results/figures/assess/call/metrics/combo_call_fnfp_{m}.png")
else:
    # Standalone execution paths
    candidate_paths = [
        Path(sys.argv[1]) if len(sys.argv) > 1 else None,
        Path("../../results/tables/assess/call-latest-with-auNGA/metrics/combo_variant_fnfp.csv"),
        Path("/scratch/project/bug_seq_scratch/qc_bench/results/tables/assess/call/metrics/combo_variant_fnfp.csv"),
        Path("/scratch/project/bug_seq_scratch/qc_bench/results/tables/assess/call/metrics/combo_variant_summary.csv"),
        Path("results/tables/assess/call/metrics/combo_variant_fnfp.csv"),
        Path("results/tables/assess/call/metrics/combo_variant_summary.csv"),
    ]
    csv_in = None
    for p in candidate_paths:
        if p is not None and p.exists():
            csv_in = p
            break

    if csv_in is None:
        raise FileNotFoundError("Could not find input variant calling CSV file. Please provide path via CLI argument: python extract_plot_fnfp.py <path_to_csv>")

    out_dir = Path("figures")
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_out = out_dir / "aggregated_fnfp_sums.csv"
    models = ["sup", "hac"]
    hue_order = None
    model_plot_map = {m: out_dir / f"combo_call_fnfp_{m}.png" for m in models}
    print(f"Running standalone: reading metrics from {csv_in}...")

# Colourblind-friendly palette from colour universal design (CUD)
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

# Load data
df = pd.read_csv(csv_in)

# Group by relevant columns and calculate absolute SUM of errors across all samples
df_grouped = df.groupby(["combo", "model", "depth", "VAR_TYPE"])[["TRUTH_FN", "QUERY_FP"]].sum().reset_index()

# Standardize depth format (e.g., '100x', '20x')
df_grouped["depth"] = df_grouped["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")

# Save aggregated sums to the CSV file
csv_out.parent.mkdir(parents=True, exist_ok=True)
df_grouped.to_csv(csv_out, index=False)
print(f"Saved aggregated FN/FP sums to {csv_out}")

# Resolve hue_order dynamically if running standalone
if hue_order is None:
    unique_depths = df_grouped["depth"].dropna().unique().tolist()
    depth_ints = [int(str(d).replace("x", "")) for d in unique_depths]
    hue_order = [f"{d}x" for d in sorted(depth_ints, reverse=True)]

df_grouped = df_grouped[df_grouped["depth"].isin(hue_order)].copy()
palette = cud(len(hue_order), start=2)

# Unified Depth legend handles matching plot_supp_figure_call.py
custom_handles = [
    Patch(facecolor=palette[i], edgecolor="black", linewidth=0.5, label=hue_order[i])
    for i in range(len(hue_order))
]

sns.set_theme(style="whitegrid")

# 2x2 grid definitions:
# Row 0: False Negatives (Col 0: SNP, Col 1: INDEL)
# Row 1: False Positives (Col 0: SNP, Col 1: INDEL)
metric_configs = [
    ("TRUTH_FN", "Total False Negatives", ["A. SNP False Negatives", "B. INDEL False Negatives"]),
    ("QUERY_FP", "Total False Positives", ["C. SNP False Positives", "D. INDEL False Positives"]),
]
var_types = ["SNP", "INDEL"]

for model in models:
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(16, 12), dpi=300, sharex=False, sharey=False)

    for r, (metric_col, y_label, panel_titles) in enumerate(metric_configs):
        for c, vartype in enumerate(var_types):
            ax = axes[r, c]
            df_sub = df_grouped.query("VAR_TYPE == @vartype and model == @model").copy()

            if not df_sub.empty:
                # Calculate sorting order specific to this facet (fewest errors first on the left)
                combo_totals = df_sub.groupby("combo")[metric_col].sum()
                facet_order = combo_totals.sort_values(ascending=True).index.tolist()

                sns.barplot(
                    data=df_sub,
                    x="combo",
                    y=metric_col,
                    hue="depth",
                    order=facet_order,
                    hue_order=hue_order,
                    palette=palette,
                    ax=ax,
                    dodge=True,
                    edgecolor="black",
                    linewidth=0.5,
                )

                ax.set_xticks(range(len(facet_order)))
                ax.set_xticklabels(facet_order, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)

            # Labels and styling
            ax.set_title(panel_titles[c], fontsize=14, loc="center", pad=10)
            ax.set_ylabel(y_label, fontsize=12)
            ax.set_xlabel("")
            ax.xaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)

            # Remove individual facet legend (figure-level legend will be added)
            if ax.get_legend() is not None:
                ax.get_legend().remove()

    # Unified Depth legend at the top center
    fig.legend(
        handles=custom_handles,
        title="Depth",
        loc="upper center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=len(hue_order),
        fontsize=12,
        title_fontsize=14,
    )

    plt.tight_layout()
    output_file = model_plot_map[model]
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {model.upper()} False Negative/Positive plot to {output_file}")

print("\nFinished plotting variant calling False Negatives and False Positives.")
