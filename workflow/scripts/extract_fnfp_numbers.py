import sys
import argparse
from pathlib import Path
import pandas as pd

if "snakemake" in locals() or "snakemake" in globals():
    input_file = Path(snakemake.input.csv)
    output_file = Path(snakemake.output.csv)
    if hasattr(snakemake, "log") and snakemake.log:
        sys.stderr = sys.stdout = open(snakemake.log[0], "w")
else:
    parser = argparse.ArgumentParser(description="Extract FN and FP numbers per sample.")
    parser.add_argument(
        "--input",
        "-i",
        default="/scratch/project/bug_seq_scratch/qc_bench/results/tables/assess/call/metrics/combo_variant_summary.csv",
        help="Input CSV path",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="/scratch/project/bug_seq_scratch/qc_bench/results/tables/assess/call/metrics/combo_variant_sample_fnfp.csv",
        help="Output CSV path",
    )
    args = parser.parse_args()
    input_file = Path(args.input)
    output_file = Path(args.output)

print(f"Reading metrics from {input_file}...")
df = pd.read_csv(input_file)

# Filter for SNP and INDEL
df = df[df["VAR_TYPE"].isin(["SNP", "INDEL"])].copy()

# Compute Total Error
df["TOTAL_ERROR"] = df["TRUTH_FN"] + df["QUERY_FP"]

# Define sort order keys
# Model order: sup (1), hac (2)
model_order = {"sup": 1, "hac": 2}
df["model_sort"] = df["model"].map(model_order).fillna(99)

# Depth order: descending numerical (e.g. 100x -> 1, 50x -> 2, 20x -> 3)
def get_depth_val(d):
    try:
        return int(str(d).rstrip("x"))
    except ValueError:
        return 0

unique_depths = sorted(df["depth"].unique(), key=get_depth_val, reverse=True)
depth_order = {d: i for i, d in enumerate(unique_depths)}
df["depth_sort"] = df["depth"].map(depth_order).fillna(99)

# VAR_TYPE order: SNP (1), INDEL (2)
var_type_order = {"SNP": 1, "INDEL": 2}
df["vartype_sort"] = df["VAR_TYPE"].map(var_type_order).fillna(99)

# Combo order: alphabetical, with unprocessed-untrimmed at the end
combos = sorted([c for c in df["combo"].unique() if c != "unprocessed-untrimmed"])
combo_order = {c: i for i, c in enumerate(combos)}
combo_order["unprocessed-untrimmed"] = 999
df["combo_sort"] = df["combo"].map(combo_order).fillna(999)

# Sort by: combo -> model -> depth -> VAR_TYPE -> TOTAL_ERROR (asc) -> TRUTH_FN (asc) -> QUERY_FP (asc) -> sample
df_sorted = df.sort_values(
    by=[
        "combo_sort",
        "model_sort",
        "depth_sort",
        "vartype_sort",
        "TOTAL_ERROR",
        "TRUTH_FN",
        "QUERY_FP",
        "sample",
    ]
)

# Select target columns
cols = [
    "combo",
    "model",
    "depth",
    "VAR_TYPE",
    "sample",
    "TRUTH_FN",
    "QUERY_FP",
    "TOTAL_ERROR",
]
df_final = df_sorted[cols].reset_index(drop=True)

output_file.parent.mkdir(parents=True, exist_ok=True)
df_final.to_csv(output_file, index=False)
print(f"Extracted per-sample FN/FP numbers and saved to {output_file}")
