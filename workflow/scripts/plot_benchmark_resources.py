import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

# Redirect stdout and stderr to the snakemake log file
sys.stderr = sys.stdout = open(snakemake.log[0], "w")

def main():
    model = getattr(snakemake.wildcards, "model", "sup")
    print(f"Aggregating benchmark resources for model: {model}")

    # 1. Calculate input throughput (Mbp) per sample from raw FASTQ files
    sample_mbp = {}
    print(f"Computing input sequencing yields from {len(snakemake.input.fastqs)} FASTQ files...")
    for fq_file in snakemake.input.fastqs:
        fq_path = Path(fq_file)
        sample = fq_path.stem
        total_bases = 0
        with open(fq_path, "r") as f:
            for i, line in enumerate(f):
                if i % 4 == 1:
                    total_bases += len(line.strip())
        sample_mbp[sample] = total_bases / 1e6
        print(f"  {sample}: {sample_mbp[sample]:.2f} Mbp")

    # 2. Parse Trimming Benchmarks
    trim_bench = {}
    if hasattr(snakemake.input, "trimming_benchmark"):
        for p_str in snakemake.input.trimming_benchmark:
            p = Path(p_str)
            trimmer = p.parts[-3]
            sample = p.name.split(".")[0]
            df_t = pd.read_csv(p, sep="\t")
            time_col = "cpu_time" if "cpu_time" in df_t.columns else "s"
            trim_bench[(trimmer, sample)] = (
                float(df_t[time_col].iloc[0]),
                float(df_t["max_rss"].iloc[0]),
            )

    # 3. Parse Quality Processing Benchmarks
    qual_bench = {}
    if hasattr(snakemake.input, "quality_benchmark"):
        for p_str in snakemake.input.quality_benchmark:
            p = Path(p_str)
            combo = p.parts[-3]
            sample = p.name.split(".")[0]
            df_q = pd.read_csv(p, sep="\t")
            time_col = "cpu_time" if "cpu_time" in df_q.columns else "s"
            qual_bench[(combo, sample)] = (
                float(df_q[time_col].iloc[0]),
                float(df_q["max_rss"].iloc[0]),
            )

    # 4. Combine resources for all combinations per sample
    rows = []
    tools = ["chopper", "filtlong", "seqkit"]
    trimmers = ["porechop_abi", "dorado", "barbell", "untrimmed"]

    for sample in sorted(sample_mbp.keys()):
        mbp = sample_mbp[sample]

        # Fastplong-all (joint trimming and quality filtering)
        if ("fastplong-all", sample) in qual_bench:
            t_s, m_rss = qual_bench[("fastplong-all", sample)]
            rows.append({
                "sample": sample,
                "model": model,
                "combo": "fastplong-all",
                "quality_tool": "fastplong",
                "trimmer": "all",
                "input_mbp": mbp,
                "cpu_time_s": t_s,
                "cpu_time_per_mbp": t_s / mbp if mbp > 0 else 0.0,
                "max_rss_mb": m_rss,
                "max_rss_gb": m_rss / 1024.0,
            })

        # Quality tools x Trimmers
        for tool in tools:
            for trimmer in trimmers:
                cname = f"{tool}-{trimmer}"
                if trimmer == "untrimmed":
                    t_s, m_rss = qual_bench.get((cname, sample), (0.0, 0.0))
                    rows.append({
                        "sample": sample,
                        "model": model,
                        "combo": cname,
                        "quality_tool": tool,
                        "trimmer": trimmer,
                        "input_mbp": mbp,
                        "cpu_time_s": t_s,
                        "cpu_time_per_mbp": t_s / mbp if mbp > 0 else 0.0,
                        "max_rss_mb": m_rss,
                        "max_rss_gb": m_rss / 1024.0,
                    })
                else:
                    t_time, t_rss = trim_bench.get((trimmer, sample), (0.0, 0.0))
                    q_time, q_rss = qual_bench.get((cname, sample), (0.0, 0.0))
                    tot_time = t_time + q_time
                    peak_rss = max(t_rss, q_rss)
                    rows.append({
                        "sample": sample,
                        "model": model,
                        "combo": cname,
                        "quality_tool": tool,
                        "trimmer": trimmer,
                        "input_mbp": mbp,
                        "cpu_time_s": tot_time,
                        "cpu_time_per_mbp": tot_time / mbp if mbp > 0 else 0.0,
                        "max_rss_mb": peak_rss,
                        "max_rss_gb": peak_rss / 1024.0,
                    })

        # Unprocessed trimmers (trimming only)
        for trimmer in ["porechop_abi", "dorado", "barbell"]:
            cname = f"unprocessed-{trimmer}"
            t_time, t_rss = trim_bench.get((trimmer, sample), (0.0, 0.0))
            rows.append({
                "sample": sample,
                "model": model,
                "combo": cname,
                "quality_tool": "unprocessed",
                "trimmer": trimmer,
                "input_mbp": mbp,
                "cpu_time_s": t_time,
                "cpu_time_per_mbp": t_time / mbp if mbp > 0 else 0.0,
                "max_rss_mb": t_rss,
                "max_rss_gb": t_rss / 1024.0,
            })

        # Unprocessed-untrimmed (raw reference baseline)
        rows.append({
            "sample": sample,
            "model": model,
            "combo": "unprocessed-untrimmed",
            "quality_tool": "unprocessed",
            "trimmer": "untrimmed",
            "input_mbp": mbp,
            "cpu_time_s": 0.0,
            "cpu_time_per_mbp": 0.0,
            "max_rss_mb": 0.0,
            "max_rss_gb": 0.0,
        })

    df = pd.DataFrame(rows)
    df.sort_values(by=["combo", "sample"], inplace=True)
    df.to_csv(snakemake.output.csv, index=False)
    print(f"Saved compiled benchmark metrics to {snakemake.output.csv}")

    # 5. Plotting Setup with CUD Color Palette
    sns.set_theme(style="whitegrid")

    trimmer_palette = {
        "untrimmed": "#e69f00",       # Orange
        "dorado": "#0072b2",          # Blue
        "all": "#cc79a7",             # Reddish Purple (fastplong)
        "barbell": "#009e73",         # Bluish Green
        "porechop_abi": "#d55e00",    # Vermilion
    }
    trimmer_order = ["untrimmed", "dorado", "all", "barbell", "porechop_abi"]
    trimmer_labels = {
        "untrimmed": "Untrimmed",
        "dorado": "Dorado",
        "all": "Fastplong",
        "barbell": "Barbell",
        "porechop_abi": "Porechop ABI",
    }

    # Two stacked panels: Top = Memory, Bottom = CPU time per Mbp
    fig, (ax_mem, ax_time) = plt.subplots(nrows=2, ncols=1, figsize=(15, 11), dpi=300)

    # Independent sorting: least to highest
    order_mem = df.groupby("combo")["max_rss_gb"].mean().sort_values().index.tolist()
    order_time = df.groupby("combo")["cpu_time_per_mbp"].mean().sort_values().index.tolist()

    # Panel A: Maximum memory usage (symlog)
    sns.boxenplot(
        data=df, x="combo", y="max_rss_gb", hue="trimmer",
        order=order_mem, hue_order=trimmer_order, palette=trimmer_palette,
        dodge=False, ax=ax_mem, showfliers=False, k_depth=3,
        box_kws=dict(alpha=0.6, edgecolor="black", linewidth=0.8),
        line_kws=dict(color="black", linewidth=1.2),
    )
    sns.stripplot(
        data=df, x="combo", y="max_rss_gb", hue="trimmer",
        order=order_mem, hue_order=trimmer_order, palette=trimmer_palette,
        dodge=False, ax=ax_mem, alpha=0.6, jitter=0.2, size=5,
        edgecolor="black", linewidth=0.5,
    )
    ax_mem.set_title("A. Maximum memory usage", fontsize=14, fontweight="normal", pad=10, loc="left")
    ax_mem.set_ylabel("Peak memory (GB)", fontsize=12, fontweight="normal")
    ax_mem.set_xlabel("")
    ax_mem.set_xticks(range(len(order_mem)))
    ax_mem.set_xticklabels(order_mem, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
    ax_mem.set_yscale("symlog", linthresh=0.05, linscale=0.5)
    ax_mem.set_ylim(bottom=0, top=20)
    ax_mem.set_yticks([0, 0.1, 0.5, 1, 2, 4, 8, 16])
    ax_mem.set_yticklabels(["0", "0.1", "0.5", "1", "2", "4", "8", "16"])
    ax_mem.xaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)
    ax_mem.yaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)
    if ax_mem.get_legend():
        ax_mem.get_legend().remove()

    # Panel B: CPU time per Mbp (symlog)
    sns.boxenplot(
        data=df, x="combo", y="cpu_time_per_mbp", hue="trimmer",
        order=order_time, hue_order=trimmer_order, palette=trimmer_palette,
        dodge=False, ax=ax_time, showfliers=False, k_depth=3,
        box_kws=dict(alpha=0.6, edgecolor="black", linewidth=0.8),
        line_kws=dict(color="black", linewidth=1.2),
    )
    sns.stripplot(
        data=df, x="combo", y="cpu_time_per_mbp", hue="trimmer",
        order=order_time, hue_order=trimmer_order, palette=trimmer_palette,
        dodge=False, ax=ax_time, alpha=0.6, jitter=0.2, size=5,
        edgecolor="black", linewidth=0.5,
    )
    ax_time.set_title("B. CPU time per Mbp", fontsize=14, fontweight="normal", pad=10, loc="left")
    ax_time.set_ylabel("CPU time (s / Mbp)", fontsize=12, fontweight="normal")
    ax_time.set_xlabel("")
    ax_time.set_xticks(range(len(order_time)))
    ax_time.set_xticklabels(order_time, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
    ax_time.set_yscale("symlog", linthresh=0.001, linscale=0.5)
    ax_time.set_ylim(bottom=0, top=45)
    ax_time.set_yticks([0, 0.001, 0.01, 0.1, 1, 10])
    ax_time.set_yticklabels(["0", "0.001", "0.01", "0.1", "1", "10"])
    ax_time.xaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)
    ax_time.yaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)
    if ax_time.get_legend():
        ax_time.get_legend().remove()

    # Unified Legend at the Top (shifted very slightly higher)
    handles = [
        Patch(facecolor=trimmer_palette[tr], edgecolor="black", label=trimmer_labels[tr], alpha=0.7)
        for tr in trimmer_order
    ]
    fig.legend(
        handles=handles, title="Adapter Trimmer", loc="upper center",
        bbox_to_anchor=(0.5, 1.008), ncol=len(trimmer_order),
        fontsize=11, title_fontsize=12, frameon=True,
    )

    plt.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(snakemake.output.figure, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved benchmark figure to {snakemake.output.figure}")

if __name__ == "__main__" or "snakemake" in globals():
    main()

