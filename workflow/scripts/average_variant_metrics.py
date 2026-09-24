import sys
import pandas as pd

# redirect both stdout and stderr to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

input_file=snakemake.input.csv
output_file=snakemake.output.csv

df = pd.read_csv(input_file)

summary = df.groupby(["combo", "model", "depth", "VAR_TYPE"])[["PREC", "RECALL", "F1_SCORE"]].mean().reset_index()
summary = summary.round(7)

# sort 100x, 50x, and 20x
summary["depth"] = pd.Categorical(summary["depth"], categories=["100x", "50x", "20x"], ordered= True)
# sort sup then hac
summary["model"] = pd.Categorical(summary["model"], categories=["sup", "hac"], ordered=True)

# Aplly sort logic
summary = summary.sort_values(by=["combo", "model", "depth", "VAR_TYPE"])

summary.to_csv(output_file, index=False)
print(f"Successfully calculated averages and saved to: {output_file}")
