import argparse
import sys
from pathlib import Path
from typing import Optional
import pandas as pd

# Ensure 'src' is importable when running standalone or within Snakemake
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "src"))

from qc_scoring.paper import (
    generate_paper_composite_results,
    generate_paper_summary_dataframe,
    plot_paper_composite_figure,
)


def run_pipeline(
    master_csv: str,
    survey_csv: Optional[str] = None,
    fig_perdepth_sup: Optional[str] = None,
    fig_perdepth_hac: Optional[str] = None,
    scores_csv_perdepth: Optional[str] = None,
    summary_csv_perdepth: Optional[str] = None,
    fig_global_sup: Optional[str] = None,
    fig_global_hac: Optional[str] = None,
    scores_csv_global: Optional[str] = None,
    summary_csv_global: Optional[str] = None,
) -> None:
    # Explicitly refuse cross-depth or aggregate all-depth outputs
    if any([fig_global_sup, fig_global_hac, scores_csv_global, summary_csv_global]):
        sys.stderr.write(
            "Error: Cross-depth or aggregate all-depth rankings are not supported under scoring specification v1.0. "
            "The paper composite is evaluated within-depth and within-model only.\n"
        )
        sys.exit(1)

    print(f"Loading master assembly metrics from: {master_csv}")
    df = pd.read_csv(master_csv)

    print("Calculating canonical Community-balanced composite results...")
    paper_results = generate_paper_composite_results(df)
    summary_df = generate_paper_summary_dataframe(df)

    for label, csv_path in [
        ("summary ranking", summary_csv_perdepth),
        ("scores", scores_csv_perdepth),
    ]:
        if csv_path:
            print(f"Saving per-depth {label} table to: {csv_path}")
            Path(csv_path).parent.mkdir(parents=True, exist_ok=True)
            summary_df.to_csv(csv_path, index=False)

    if fig_perdepth_sup:
        print(f"Generating per-depth composite figure for SUP: {fig_perdepth_sup}")
        plot_paper_composite_figure(paper_results, model="sup", out_path=fig_perdepth_sup)

    if fig_perdepth_hac:
        print(f"Generating per-depth composite figure for HAC: {fig_perdepth_hac}")
        plot_paper_composite_figure(paper_results, model="hac", out_path=fig_perdepth_hac)

    print("Aggregate assembly assessment processing complete.")


def main() -> None:
    # Determine if invoked via Snakemake
    if "snakemake" in globals():
        snakemake_obj = globals()["snakemake"]
        if hasattr(snakemake_obj, "log") and snakemake_obj.log:
            sys.stderr = sys.stdout = open(snakemake_obj.log[0], "w")

        master_csv = snakemake_obj.input.master_csv
        survey_csv = getattr(snakemake_obj.input, "survey_csv", None)

        fig_perdepth_sup = getattr(snakemake_obj.output, "fig_perdepth_sup", None)
        fig_perdepth_hac = getattr(snakemake_obj.output, "fig_perdepth_hac", None)
        scores_csv_perdepth = getattr(snakemake_obj.output, "scores_csv_perdepth", None)
        summary_csv_perdepth = getattr(snakemake_obj.output, "summary_csv_perdepth", None)

        fig_global_sup = getattr(snakemake_obj.output, "fig_global_sup", None)
        fig_global_hac = getattr(snakemake_obj.output, "fig_global_hac", None)
        scores_csv_global = getattr(snakemake_obj.output, "scores_csv_global", None)
        summary_csv_global = getattr(snakemake_obj.output, "summary_csv_global", None)

        run_pipeline(
            master_csv=master_csv,
            survey_csv=survey_csv,
            fig_perdepth_sup=fig_perdepth_sup,
            fig_perdepth_hac=fig_perdepth_hac,
            scores_csv_perdepth=scores_csv_perdepth,
            summary_csv_perdepth=summary_csv_perdepth,
            fig_global_sup=fig_global_sup,
            fig_global_hac=fig_global_hac,
            scores_csv_global=scores_csv_global,
            summary_csv_global=summary_csv_global,
        )
    else:
        parser = argparse.ArgumentParser(
            description="Generate survey-weighted assembly composite rankings and figures (scoring v1.0)."
        )
        parser.add_argument("--master-csv", required=True, help="Path to assembly metrics master CSV.")
        parser.add_argument("--survey-csv", required=False, help="Path to microbial QC survey CSV (optional).")
        parser.add_argument("--fig-perdepth-sup", required=False, help="Output path for SUP faceted figure.")
        parser.add_argument("--fig-perdepth-hac", required=False, help="Output path for HAC faceted figure.")
        parser.add_argument("--scores-csv-perdepth", required=False, help="Output path for per-depth scores CSV.")
        parser.add_argument("--summary-csv-perdepth", required=False, help="Output path for per-depth summary CSV.")

        # Deprecated global targets to explicitly refuse if requested
        parser.add_argument("--fig-global-sup", required=False, help=argparse.SUPPRESS)
        parser.add_argument("--fig-global-hac", required=False, help=argparse.SUPPRESS)
        parser.add_argument("--scores-csv-global", required=False, help=argparse.SUPPRESS)
        parser.add_argument("--summary-csv-global", required=False, help=argparse.SUPPRESS)

        args = parser.parse_args()
        run_pipeline(
            master_csv=args.master_csv,
            survey_csv=args.survey_csv,
            fig_perdepth_sup=args.fig_perdepth_sup,
            fig_perdepth_hac=args.fig_perdepth_hac,
            scores_csv_perdepth=args.scores_csv_perdepth,
            summary_csv_perdepth=args.summary_csv_perdepth,
            fig_global_sup=args.fig_global_sup,
            fig_global_hac=args.fig_global_hac,
            scores_csv_global=args.scores_csv_global,
            summary_csv_global=args.summary_csv_global,
        )


if __name__ == "__main__":
    main()
elif "snakemake" in globals():
    main()
