import sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import pandas as pd
import seaborn as sns
import numpy as np

# Redirect all prints and errors to log file
if hasattr(snakemake, "log") and snakemake.log:
    sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Generating quality filtering parameter exploration figure...")

CSV_CALL = snakemake.input.variant_csv
CSV_QUAST = snakemake.input.quast_csv
CSV_MISSED = snakemake.input.missed_csv
CSV_CONTAM = snakemake.input.contam_csv
OUT_FIG = snakemake.output.figure

# Colourblind-friendly palette from Colour Universal Design (CUD)
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

# 100x -> skyblue (#56b4e9), 20x -> vermilion (#d55e00) matching main paper figures
hue_order = ["100x", "20x"]
palette = [named_colors["skyblue"], named_colors["vermilion"]]

def clean_tool_name(combo_series):
    return combo_series.astype(str).str.replace(r"-dorado$", "", regex=True)

def main():
    sns.set_theme(style="whitegrid")
    
    model = getattr(snakemake.wildcards, "model", "sup")
    print(f"Plotting model: {model}")

    # --- 1. Load and prepare variant calling data ---
    df_call = pd.read_csv(CSV_CALL)
    df_call = df_call.query("model == @model and VAR_TYPE in ('SNP', 'INDEL')").copy()
    df_call["depth"] = df_call["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df_call = df_call[df_call["depth"].isin(hue_order)].copy()
    df_call["tool"] = clean_tool_name(df_call["combo"])

    # --- 2. Load and prepare assembly QUAST data ---
    df_quast = pd.read_csv(CSV_QUAST).query("model == @model").copy()
    df_quast["depth"] = df_quast["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df_quast = df_quast[df_quast["depth"].isin(hue_order)].copy()
    df_quast["tool"] = clean_tool_name(df_quast["combo"])
    df_quast["Total_Errors"] = df_quast["Mismatches per 100kbp"] + df_quast["Indels per 100kbp"]
    df_quast["auNGA_score"] = np.maximum(0, 1 - np.abs(df_quast["auNGA_ratio"] - 1.0))

    # --- 3. Load and prepare missed contigs data ---
    df_missed = pd.read_csv(CSV_MISSED).query("model == @model").copy()
    df_missed["depth"] = df_missed["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df_missed = df_missed[df_missed["depth"].isin(hue_order)].copy()
    df_missed["tool"] = clean_tool_name(df_missed["combo"])
    if "total_missed" not in df_missed.columns and "full_missed" in df_missed.columns:
        df_missed["total_missed"] = df_missed["full_missed"] + df_missed.get("partial_missed", 0)
    df_missed_agg = df_missed.groupby(["tool", "depth"], as_index=False)["total_missed"].sum()

    # --- 4. Load and prepare contaminant data ---
    df_contam = pd.read_csv(CSV_CONTAM).query("model == @model").copy()
    df_contam["depth"] = df_contam["depth"].apply(lambda d: f"{str(d).rstrip('x')}x")
    df_contam = df_contam[df_contam["depth"].isin(hue_order)].copy()
    df_contam["tool"] = clean_tool_name(df_contam["combo"])
    df_contam_agg = df_contam.groupby(["tool", "depth"], as_index=False)["contamination_count"].sum()

    # Establish consistent tool display order (excluding 1000 bp settings)
    tools_in_data = sorted(list(set(
        df_call["tool"].unique().tolist() +
        df_quast["tool"].unique().tolist() +
        df_missed_agg["tool"].unique().tolist()
    )))
    canonical_order = [
        "chopper_extract100", "chopper_trim100", "fastplong_100",
        "filtlong_default", "filtlong_len", "filtlong_meanq",
        "seqkit_100", "unprocessed"
    ]
    order = [t for t in canonical_order if t in tools_in_data] or tools_in_data

    # Setup 2x3 grid
    fig, axes = plt.subplots(2, 3, figsize=(18, 11), dpi=300)
    axes = axes.flatten()

    def format_x_axis(ax):
        ax.set_xlabel("")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
        ax.xaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)

    # --- Panel A: SNP F1 Score ---
    ax_a = axes[0]
    df_snp = df_call.query("VAR_TYPE == 'SNP'").copy()
    cap = 0.99999
    df_snp["plot_val"] = df_snp["F1_SCORE"].apply(lambda v: cap if v > cap else v)
    sns.stripplot(
        data=df_snp, x="tool", y="plot_val", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_a,
        alpha=0.4, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=4
    )
    sns.pointplot(
        data=df_snp, x="tool", y="plot_val", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_a,
        dodge=0.3, errorbar=("ci", 95), capsize=0.1,
        err_kws={"linewidth": 1}, linewidth=1, markersize=5, estimator=np.mean, legend=False, zorder=2
    )
    ax_a.set_title("A. SNP F1 Score", fontsize=13, pad=8)
    ax_a.set_ylabel("SNP F1 score", fontsize=11)
    yticks = [0.5, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
    ax_a.set_yscale("logit", nonpositive="clip")
    ax_a.set_yticks(yticks)
    ax_a.set_yticklabels([f"{y*100:g}%" if y < cap else "100%" for y in yticks])
    ax_a.set_ylim(bottom=0.5)
    format_x_axis(ax_a)
    if ax_a.get_legend(): ax_a.get_legend().remove()

    # --- Panel B: INDEL F1 Score ---
    ax_b = axes[1]
    df_indel = df_call.query("VAR_TYPE == 'INDEL'").copy()
    df_indel["plot_val"] = df_indel["F1_SCORE"].apply(lambda v: cap if v > cap else v)
    sns.stripplot(
        data=df_indel, x="tool", y="plot_val", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_b,
        alpha=0.4, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=4
    )
    sns.pointplot(
        data=df_indel, x="tool", y="plot_val", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_b,
        dodge=0.3, errorbar=("ci", 95), capsize=0.1,
        err_kws={"linewidth": 1}, linewidth=1, markersize=5, estimator=np.mean, legend=False, zorder=2
    )
    ax_b.set_title("B. INDEL F1 Score", fontsize=13, pad=8)
    ax_b.set_ylabel("INDEL F1 score", fontsize=11)
    ax_b.set_yscale("logit", nonpositive="clip")
    ax_b.set_yticks(yticks)
    ax_b.set_yticklabels([f"{y*100:g}%" if y < cap else "100%" for y in yticks])
    ax_b.set_ylim(bottom=0.5)
    format_x_axis(ax_b)
    if ax_b.get_legend(): ax_b.get_legend().remove()

    # --- Panel C: Total Errors per 100 kbp ---
    ax_c = axes[2]
    sns.stripplot(
        data=df_quast, x="tool", y="Total_Errors", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_c,
        alpha=0.4, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=4
    )
    sns.pointplot(
        data=df_quast, x="tool", y="Total_Errors", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_c,
        dodge=0.3, errorbar=("ci", 95), capsize=0.1,
        err_kws={"linewidth": 1}, linewidth=1, markersize=5, estimator=np.mean, legend=False, zorder=2
    )
    ax_c.set_title("C. Assembly Errors per 100kbp\n(Mismatches + Indels)", fontsize=13, pad=8)
    ax_c.set_ylabel("Total errors per 100kbp", fontsize=11)
    format_x_axis(ax_c)
    if ax_c.get_legend(): ax_c.get_legend().remove()

    # --- Panel D: auNGA Contiguity Score ---
    ax_d = axes[3]
    df_quast["plot_aunga"] = df_quast["auNGA_score"].apply(lambda v: cap if v >= cap else (0.00001 if v <= 0 else v))
    sns.stripplot(
        data=df_quast, x="tool", y="plot_aunga", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_d,
        alpha=0.4, dodge=True, linewidth=0.5, edgecolor="black", zorder=1, size=4
    )
    sns.pointplot(
        data=df_quast, x="tool", y="plot_aunga", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_d,
        dodge=0.3, errorbar=("ci", 95), capsize=0.1,
        err_kws={"linewidth": 1}, linewidth=1, markersize=5, estimator=np.mean, legend=False, zorder=2
    )
    ax_d.set_title("D. Assembly Contiguity\n(auNGA Score)", fontsize=13, pad=8)
    ax_d.set_ylabel("auNGA score (1 - |ratio - 1|)", fontsize=11)
    yticks_aunga = [0.5, 0.7, 0.8, 0.9, 0.99, 0.999, 0.9999, cap]
    ax_d.set_yscale("logit", nonpositive="clip")
    ax_d.set_yticks(yticks_aunga)
    ax_d.set_yticklabels([f"{y*100:g}%" if y < cap else "100%" for y in yticks_aunga])
    ax_d.set_ylim(bottom=0.5)
    format_x_axis(ax_d)
    if ax_d.get_legend(): ax_d.get_legend().remove()

    # --- Panel E: Total Missed Contigs ---
    ax_e = axes[4]
    sns.barplot(
        data=df_missed_agg, x="tool", y="total_missed", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_e,
        edgecolor="black", linewidth=0.5
    )
    ax_e.set_title("E. Missed Contigs", fontsize=13, pad=8)
    ax_e.set_ylabel("Total missed contigs", fontsize=11)
    ax_e.yaxis.set_major_locator(MaxNLocator(integer=True))
    format_x_axis(ax_e)
    if ax_e.get_legend(): ax_e.get_legend().remove()

    # --- Panel F: Total Contamination Count ---
    ax_f = axes[5]
    sns.barplot(
        data=df_contam_agg, x="tool", y="contamination_count", hue="depth",
        order=order, hue_order=hue_order, palette=palette, ax=ax_f,
        edgecolor="black", linewidth=0.5
    )
    ax_f.set_title("F. Contamination Count", fontsize=13, pad=8)
    ax_f.set_ylabel("Total contaminants", fontsize=11)
    ax_f.yaxis.set_major_locator(MaxNLocator(integer=True))
    format_x_axis(ax_f)
    if ax_f.get_legend(): ax_f.get_legend().remove()

    plt.tight_layout()

    # Top unified legend
    custom_handles = [Patch(facecolor=palette[i], edgecolor="black", label=hue_order[i]) for i in range(len(hue_order))]
    fig.legend(
        handles=custom_handles, title="Depth", loc="lower center",
        bbox_to_anchor=(0.5, 1.02), ncol=len(hue_order), fontsize=12, title_fontsize=13
    )

    Path(OUT_FIG).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved exploration figure to: {OUT_FIG}")

if __name__ == "__main__":
    main()
