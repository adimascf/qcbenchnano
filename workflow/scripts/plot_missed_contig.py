import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
import seaborn as sns
from pathlib import Path
from typing import List

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Compiling full and partial missed contig data...")

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

# Compile the Data on the fly
records = []
input_files = snakemake.input.tsv if isinstance(snakemake.input.tsv, list) else [snakemake.input.tsv]

for f in input_files:
    path = Path(f)
    model = path.parent.name                 
    depth = path.parent.parent.name          
    combo = path.parent.parent.parent.name   
    sample = path.name.split('.')[0]         

    try:
        df = pd.read_csv(f, sep='\t')
        
        if not df.empty and 'coverage' in df.columns:
            # get the sequence name and length columns
            rname_col = '#rname' if '#rname' in df.columns else df.columns[0]
            len_col = 'endpos' if 'endpos' in df.columns else None

            full_mask = df['coverage'] == 0
            partial_mask = (df['coverage'] > 0) & (df['coverage'] < 50)

            full_missed = full_mask.sum()
            partial_missed = partial_mask.sum()

            # Extract names, lengths, and coverage for ALL contigs
            missed_info = []
            all_info = []

            for _, row in df.iterrows():
                name = row[rname_col]
                cov = row['coverage']
                length_str = f"{int(row[len_col])}bp" if len_col and pd.notna(row[len_col]) else "unknown length"

                # format string: e.g., "plasmid_3 (3173bp, 0% cov)"
                info_str = f"{name} ({length_str}, {cov}% cov)"
                
                # Append to the total record list
                all_info.append(info_str)
                
                # If it meets the threshold, also append it to the missed list
                if cov < 50:
                    missed_info.append(info_str)

            missed_names_str = "; ".join(missed_info)
            all_names_str = "; ".join(all_info)
        else:
            full_missed, partial_missed = 0, 0
            missed_names_str = ""
            all_names_str = ""
            
    except pd.errors.EmptyDataError:
        full_missed, partial_missed = 0, 0
        missed_names_str = ""
        all_names_str = ""

    records.append({
        "sample": sample,
        "combo": combo,
        "depth": depth,
        "model": model,
        "full_missed": full_missed,
        "partial_missed": partial_missed,
        "total_missed": full_missed + partial_missed,
        "missed_contig_details": missed_names_str,
        "all_contigs_coverage": all_names_str  # Added the new column here
    })

# Convert to pandas DataFrame and save the compiled table
missed_df = pd.DataFrame(records)
missed_df.to_csv(snakemake.output.table, index=False)
print(f"Saved compiled data table to {snakemake.output.table}")

# Plot 1: 2x2 Grid (Individual Sample Points)
melted_df = missed_df.melt(
    id_vars=["sample", "combo", "depth", "model"],
    value_vars=["full_missed", "partial_missed"],
    var_name="miss_type",
    value_name="count"
)

sns.set_theme(style="whitegrid")
models = ["sup", "hac"]

# Fetch the depths dynamically from the Snakemake config
config_depths = snakemake.config['depth']

# sort the depths numerically, and add 'x'
hue_order = [f"{d}x" for d in sorted([int(d) for d in config_depths], reverse=True)]

miss_types = ["full_missed", "partial_missed"]
miss_labels = {"full_missed": "Full Miss (0% Cov)", "partial_missed": "Partial Miss (<50% Cov)"}

# Sort combos by lowest mean total missed contigs (Best to Worst overall)
combo_perf = missed_df.groupby("combo")["total_missed"].mean()
order_abs = combo_perf.sort_values(ascending=True).index.tolist()

# Define the target samples and their unique markers
target_markers = {
    "ATCC_25922__202309": "^",  
    "BPH2947__202310": "s"      
}

fig, axes = plt.subplots(nrows=2, ncols=len(models), figsize=(16, 10), dpi=300, sharex=True, sharey=True)

for r, m_type in enumerate(miss_types):
    for c, model in enumerate(models):
        ax = axes[r, c]
        df_sub = melted_df.query("model == @model and miss_type == @m_type").copy()

        if not df_sub.empty:
            draw_legend = (r == 0 and c == len(models) - 1) 
            
            # 1. Base Layer: Plot all OTHER samples (lowered alpha to 0.3)
            sns.stripplot(
                data=df_sub[~df_sub["sample"].isin(target_markers.keys())], 
                x="combo", y="count", hue="depth", 
                order=order_abs, hue_order=hue_order,
                palette=cud(len(hue_order), start=2), 
                ax=ax, alpha=0.3, dodge=True, linewidth=0.5, edgecolor="black", 
                legend=draw_legend, marker="o"
            )

            # 2. Overlay Layer: Plot target samples (lowered alpha to 0.6)
            for sample_name, marker_shape in target_markers.items():
                target_df = df_sub[df_sub["sample"] == sample_name]
                if not target_df.empty:
                    sns.stripplot(
                        data=target_df,
                        x="combo", y="count", hue="depth",
                        order=order_abs, hue_order=hue_order,
                        palette=cud(len(hue_order), start=2),
                        ax=ax, alpha=0.6, dodge=True, linewidth=0.8, edgecolor="black", 
                        legend=False, marker=marker_shape, s=7 
                    )

        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        if c == 0:
            ax.set_ylabel(miss_labels[m_type])
        else:
            ax.set_ylabel("")

        if r == 0:
            ax.set_title(f"{model.upper()} Model")
        
        ax.set_xlabel("")
        ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

        if r == 1:
            ax.set_xticks(range(len(order_abs)))
            ax.set_xticklabels(order_abs, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
        else:
            ax.tick_params(labelbottom=False)

        # Build Custom Legend for Shapes & Colors with exact sample names
        if r == 0 and c == len(models) - 1 and ax.get_legend() is not None:
            handles, labels = ax.get_legend_handles_labels()
            
            base_handle = Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='black', markersize=8, label='Other Samples')
            sample1_handle = Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', markeredgecolor='black', markersize=8, label='ATCC_25922__202309')
            sample2_handle = Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', markeredgecolor='black', markersize=8, label='BPH2947__202310')
            
            ax.legend(handles=handles + [base_handle, sample1_handle, sample2_handle], title="Depth & Samples", bbox_to_anchor=(1.05, 1), loc='upper left')

fig.tight_layout()
fig.savefig(snakemake.output.figures, bbox_inches='tight')
plt.close(fig)
print(f"Saved original 2x2 plot to {snakemake.output.figures}")

# Plot 2: 1x2 Grid (Pooled across the samples and missed contig types)
# Calculate the sum of total missed contigs across all samples for each tool combo/depth
pooled_df = missed_df.groupby(["combo", "model", "depth"])["total_missed"].sum().reset_index()

fig2, axes2 = plt.subplots(nrows=1, ncols=len(models), figsize=(14, 6), dpi=300, sharex=True, sharey=True)

# Safety check for single axes
if not isinstance(axes2, (np.ndarray, list)):
    axes2 = [axes2]

for c, model in enumerate(models):
    ax = axes2[c]
    df_sub = pooled_df.query("model == @model").copy()

    if not df_sub.empty:
        sns.stripplot(
            data=df_sub, x="combo", y="total_missed", hue="depth", 
            order=order_abs, hue_order=hue_order,
            palette=cud(len(hue_order), start=2), 
            ax=ax, alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black"
        )

    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    
    if c == 0:
        ax.set_ylabel("Total Missed Contigs")
    else:
        ax.set_ylabel("")

    ax.set_title(f"Total Missed Contigs - {model.upper()} Model")
    ax.set_xlabel("")
    ax.xaxis.grid(True, linestyle='--', color='lightgrey', zorder=0)

    ax.set_xticks(range(len(order_abs)))
    ax.set_xticklabels(order_abs, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)

    # Manage legends
    if c == len(models) - 1 and ax.get_legend() is not None:
        ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
    elif ax.get_legend() is not None:
        ax.get_legend().remove()

fig2.tight_layout()
fig2.savefig(snakemake.output.figures_total, bbox_inches='tight')
plt.close(fig2)
print(f"Saved pooled totals plot to {snakemake.output.figures_total}")
