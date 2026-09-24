import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Parsing {len(snakemake.input.reports)} QUAST reports...")

assembly_metrics = []
for file_path in snakemake.input.reports:
    try:
        p = Path(file_path)
        
        # Path: ../quast/<combo>/<depth>x/<model>/<sample>.<combo>.report.tsv
        model = p.parts[-2]
        depth = p.parts[-3]
        combo = p.parts[-4]
        
        # Extract the sample name from the file name (everything before the first '.')
        sample = p.name.split('.')[0]

        # index_col=0 sets the metric names as the row index
        df_quast = pd.read_csv(p, sep='\t', index_col=0)
        col_name = df_quast.columns[0]

        mismatches = float(df_quast.loc['# mismatches per 100 kbp', col_name])
        indels = float(df_quast.loc['# indels per 100 kbp', col_name])
        nga = int(df_quast.loc['NGA50', col_name])
        misassemblies = int(df_quast.loc['# misassemblies', col_name])

        assembly_metrics.append({
            'combo': combo, 'depth': depth, 'sample': sample, 'model': model,
            'Mismatches per 100kbp': mismatches,
            'Indels per 100kbp': indels,
            'NGA50': nga,
            'misassemblies': misassemblies
        })

    except Exception as e:
        print(f"Skipping {file_path} due to error: {e}")

df = pd.DataFrame(assembly_metrics)
df.sort_values(by=["combo", "sample", "model", "depth"], inplace=True)
df.to_csv(snakemake.output.csv, index=False)


# Calcualte baseline using delta
baseline = df[df['combo'] == 'unprocessed-untrimmed'].set_index(['sample', 'depth', 'model'])[['Mismatches per 100kbp', 'Indels per 100kbp', 'NGA50']]

df_merged = df.set_index(['sample', 'depth', 'model']).join(baseline, rsuffix='_baseline').reset_index()

df_merged['Delta_Mismatches'] = df_merged['Mismatches per 100kbp_baseline'] - df_merged['Mismatches per 100kbp']
df_merged['Delta_Indels'] = df_merged['Indels per 100kbp_baseline'] - df_merged['Indels per 100kbp']

# Calculate Log2FC (Handling zeroes gracefully if NGA50 is ever exactly 0)
df_merged['Log2FC_NGA50'] = np.log2((df_merged['NGA50'] + 1) / (df_merged['NGA50_baseline'] + 1))

# Remove unprocessed-untrimmed from the delta plots
df_delta = df_merged[df_merged['combo'] != 'unprocessed-untrimmed'].copy()


## PLOTTING (ABSOLUTE SNP AND MISMATCHES AND NGA50)
sns.set_theme(style="whitegrid")
models = ["sup", "hac"]

# Fetch the depths dynamically from the Snakemake config
config_depths = snakemake.config['depth']

# sort the depths numerically, and add 'x'
hue_order = [f"{d}x" for d in sorted([int(d) for d in config_depths], reverse=True)]

order_abs = sorted(df["combo"].unique())
order_delta = sorted(df_delta['combo'].unique())

# Absolute Mismatches and Indels
metrics_abs = ["Mismatches per 100kbp", "Indels per 100kbp"]
fig1, axes1 = plt.subplots(nrows=len(metrics_abs), ncols=len(models), figsize=(14, 10), dpi=300, sharex=True, sharey='row')

for r, metric in enumerate(metrics_abs):
    for c, model in enumerate(models):
        ax = axes1[r, c]
        df_sub = df.query("model == @model")

        if not df_sub.empty:
            sns.boxplot(data=df_sub, x="combo", y=metric, hue="depth", order=order_abs, hue_order=hue_order, ax=ax, fliersize=0, gap=0.1)
            sns.stripplot(data=df_sub, x="combo", y=metric, hue="depth", order=order_abs, hue_order=hue_order, ax=ax, alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=False)

        ax.set_ylabel(metric if c == 0 else "")
        ax.set_title(f"{model} model" if r == 0 else "")
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelrotation=45, labelsize=11) if r == len(metrics_abs) - 1 else ax.tick_params(axis="x", bottom=False)

        if r == 0 and c == 1 and ax.get_legend() is not None:
            ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
        elif ax.get_legend() is not None:
            ax.get_legend().remove()

fig1.tight_layout()
fig1.savefig(snakemake.output.error_plot, bbox_inches='tight')

# Absolute NGA50
fig2, axes2 = plt.subplots(nrows=1, ncols=len(models), figsize=(14, 5.5), dpi=300, sharex=True, sharey=True)

for c, model in enumerate(models):
    ax = axes2[c]
    df_sub = df.query("model == @model")

    if not df_sub.empty:
        sns.boxplot(data=df_sub, x="combo", y="NGA50", hue="depth", order=order_abs, hue_order=hue_order, ax=ax, fliersize=0, gap=0.1)
        sns.stripplot(data=df_sub, x="combo", y="NGA50", hue="depth", order=order_abs, hue_order=hue_order, ax=ax, alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=False)

    ax.set_ylabel("NGA50" if c == 0 else "")
    ax.set_title(f"NGA50 - {model} model")
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelrotation=45, labelsize=11)

    if c == 1 and ax.get_legend() is not None:
        ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
    elif ax.get_legend() is not None:
        ax.get_legend().remove()

fig2.tight_layout()
fig2.savefig(snakemake.output.nga50_plot, bbox_inches='tight')


## IMPROVEMENT PLOT
# Using delta for mismatches and indels
metrics_delta = ["Delta_Mismatches", "Delta_Indels"]
fig3, axes3 = plt.subplots(nrows=len(metrics_delta), ncols=len(models), figsize=(14, 10), dpi=300, sharex=True, sharey='row')

for r, metric in enumerate(metrics_delta):
    for c, model in enumerate(models):
        ax = axes3[r, c]
        df_sub = df_delta.query("model == @model")

        ax.axhline(0, color='red', linestyle='--', linewidth=1.5, zorder=0, label="Untrimmed Baseline")

        if not df_sub.empty:
            sns.boxplot(data=df_sub, x='combo', y=metric, hue='depth', order=order_delta, hue_order=hue_order, ax=ax, fliersize=0, gap=0.1)
            sns.stripplot(data=df_sub, x='combo', y=metric, hue='depth', order=order_delta, hue_order=hue_order, ax=ax, alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=False)

        clean_name = metric.replace("Delta_", "Errors Prevented: ")
        ax.set_ylabel(clean_name if c == 0 else "")
        ax.set_title(f"{model} model" if r == 0 else "")
        ax.set_xlabel("")
        ax.tick_params(axis="x", labelrotation=45, labelsize=11) if r == len(metrics_delta) - 1 else ax.tick_params(axis="x", bottom=False)

        if r == 0 and c == 1 and ax.get_legend() is not None:
            ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
        elif ax.get_legend() is not None:
            ax.get_legend().remove()

fig3.tight_layout()
fig3.savefig(snakemake.output.delta_error_plot, bbox_inches='tight')

# Using Log2FC for NGA50
fig4, axes4 = plt.subplots(nrows=1, ncols=len(models), figsize=(14, 6.0), dpi=300, sharex=True, sharey=True)

# Find absolute max Log2FC globally to make the Y-axis perfectly symmetric for both sup and hac
max_abs_val = df_delta['Log2FC_NGA50'].abs().max()
y_limit = max_abs_val * 1.1

for c, model in enumerate(models):
    ax = axes4[c]
    df_sub = df_delta.query("model == @model")

    ax.axhline(0, color='red', linestyle='--', linewidth=1.5, zorder=0, label="Untrimmed Baseline (Log2FC = 0)")
    ax.set_ylim(-y_limit, y_limit)

    if not df_sub.empty:
        sns.boxplot(data=df_sub, x='combo', y='Log2FC_NGA50', hue='depth', order=order_delta, hue_order=hue_order, ax=ax, fliersize=0, gap=0.1)
        sns.stripplot(data=df_sub, x='combo', y='Log2FC_NGA50', hue='depth', order=order_delta, hue_order=hue_order, ax=ax, alpha=0.6, dodge=True, linewidth=0.5, edgecolor="black", legend=False)

    ax.set_ylabel("Log2 Fold Change" if c == 0 else "")
    ax.set_title(f"Symmetric NGA50 Improv/Damage - {model}")
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelrotation=45, labelsize=11)

    if c == 1 and ax.get_legend() is not None:
        ax.legend(title="Depth", bbox_to_anchor=(1.05, 1), loc='upper left')
    elif ax.get_legend() is not None:
        ax.get_legend().remove()

fig4.tight_layout()
fig4.savefig(snakemake.output.delta_nga50_plot, bbox_inches='tight')

print("Finished plotting all assembly metrics.")
