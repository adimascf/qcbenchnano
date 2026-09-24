#!/usr/bin/env python3
"""
generate_master_csv.py

Compiles the complete master assembly metrics matrix by merging:
  1. Compiled QUAST metrics (with auNGA_ratio)
  2. Contaminant summary counts
  3. Missed contigs summary

Default output is saved to logbook/assembly_metrics.csv.
"""

import sys
import argparse
from pathlib import Path
import pandas as pd

# Default table locations
DEFAULT_RESULTS_DIR = Path("/scratch/project_mnt/S0256/qc_bench/results/tables/assess/assembly/metrics")
DEFAULT_QUAST = DEFAULT_RESULTS_DIR / "combo_quast_compiled_metrics.csv"
DEFAULT_CONTAM = DEFAULT_RESULTS_DIR / "combo_contaminant_summary_count.csv"
DEFAULT_MISSED = DEFAULT_RESULTS_DIR / "combo_assembly_missed_contigs.csv"
DEFAULT_OUTPUT = Path(__file__).resolve().parent.parent / "logbook" / "assembly_metrics.csv"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compile and recreate the master assembly metrics matrix CSV."
    )
    parser.add_argument(
        "--quast",
        type=Path,
        default=DEFAULT_QUAST,
        help=f"Path to compiled QUAST metrics CSV (default: {DEFAULT_QUAST})"
    )
    parser.add_argument(
        "--contaminants",
        type=Path,
        default=DEFAULT_CONTAM,
        help=f"Path to contaminant count CSV (default: {DEFAULT_CONTAM})"
    )
    parser.add_argument(
        "--missed-contigs",
        type=Path,
        default=DEFAULT_MISSED,
        help=f"Path to missed contigs CSV (default: {DEFAULT_MISSED})"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Path to save compiled master CSV (default: {DEFAULT_OUTPUT})"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Validate input files
    for name, path in [
        ("QUAST metrics", args.quast),
        ("Contaminants summary", args.contaminants),
        ("Missed contigs summary", args.missed_contigs),
    ]:
        if not path.exists():
            print(f"Error: {name} file not found at '{path}'", file=sys.stderr)
            sys.exit(1)

    print(f"Loading QUAST metrics from: {args.quast}")
    quast_metrics = pd.read_csv(args.quast)

    print(f"Loading contaminants from: {args.contaminants}")
    contam_metrics = pd.read_csv(args.contaminants)

    print(f"Loading missed contigs from: {args.missed_contigs}")
    miscontig_metrics = pd.read_csv(args.missed_contigs)

    print("Merging metrics on ['combo', 'depth', 'sample', 'model']...")
    assembly_metrics = quast_metrics.merge(
        contam_metrics, on=['combo', 'depth', 'sample', 'model'], how='inner'
    ).merge(
        miscontig_metrics, on=['combo', 'depth', 'sample', 'model'], how='inner'
    )

    # Sort rows consistently
    sort_cols = [col for col in ['combo', 'sample', 'model', 'depth'] if col in assembly_metrics.columns]
    if sort_cols:
        assembly_metrics.sort_values(by=sort_cols, inplace=True)

    # Ensure output parent directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    assembly_metrics.to_csv(args.output, index=False)
    print(f"Successfully compiled {len(assembly_metrics)} rows into master matrix: {args.output}")


if __name__ == "__main__":
    main()
