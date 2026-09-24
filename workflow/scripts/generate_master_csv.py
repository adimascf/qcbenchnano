import sys
from pathlib import Path
import pandas as pd

# Redirect all prints and errors to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print("Compiling all assembly metrics...")

quast_metrics = pd.read_csv(snakemake.input.quast)
contam_metrics = pd.read_csv(snakemake.input.contaminants)
miscontig_metrics = pd.read_csv(snakemake.input.missed_contigs)

print(f"QUAST entries: {len(quast_metrics)}, Contaminants entries: {len(contam_metrics)}, Missed contigs entries: {len(miscontig_metrics)}")

assembly_metrics = quast_metrics.merge(
    contam_metrics, on=['combo', 'depth', 'sample', 'model'], how='inner'
).merge(
    miscontig_metrics, on=['combo', 'depth', 'sample', 'model'], how='inner'
)

# Sort rows consistently
sort_cols = [col for col in ['combo', 'sample', 'model', 'depth'] if col in assembly_metrics.columns]
if sort_cols:
    assembly_metrics.sort_values(by=sort_cols, inplace=True)

# Ensure parent directory exists
Path(snakemake.output.master).parent.mkdir(parents=True, exist_ok=True)
assembly_metrics.to_csv(snakemake.output.master, index=False)
print(f"Saved master assembly metrics table ({len(assembly_metrics)} rows) to: {snakemake.output.master}")

