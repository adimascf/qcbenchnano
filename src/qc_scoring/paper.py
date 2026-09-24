"""Paper scoring interface and visualization for scoring specification v1.0.

Provides canonical Community-balanced composite results and figures for paper analyses,
explicitly enforcing separate scenario reporting and survey-weighted preference alignment semantics.
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from qc_scoring.models import Scenario
from qc_scoring.preferences import GateConfig, WeightsConfig
from qc_scoring.scorer import (
    SCORING_VERSION,
    SOURCE_DATA_BASELINE_COMMIT,
    ScoringResult,
    score_benchmark,
)

ANALYSIS_DESCRIPTION = "survey-weighted decision analysis (preference alignment)"

COMMUNITY_BALANCED_WEIGHTS = WeightsConfig(
    accuracy=28.0,
    contiguity=20.0,
    residual=17.0,
    replicon=35.0,
)

COMMUNITY_BALANCED_GATES = GateConfig(
    complete_recovery=False,
    zero_residual_hits=False,
)

SUPPORTED_PAPER_SCENARIOS: List[Scenario] = [
    Scenario("hac", "20x"),
    Scenario("hac", "100x"),
    Scenario("sup", "20x"),
    Scenario("sup", "100x"),
]


def generate_paper_composite_results(
    df: pd.DataFrame,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> Dict[Scenario, ScoringResult]:
    """Generate canonical Community-balanced scoring results for each supported paper scenario.

    Explicitly rejects unsupported scenarios and aggregate all-depth or all-model requests.
    """
    if "model" not in df.columns or "depth" not in df.columns:
        raise ValueError("Benchmark dataframe must contain 'model' and 'depth' columns.")

    # Normalise depth column for checking (e.g. 20 -> 20x)
    df_depths = set(df["depth"].astype(str).str.rstrip("x") + "x")
    df_models = set(df["model"].astype(str).str.lower())

    valid_depths = {"20x", "100x"}
    valid_models = {"hac", "sup"}
    unsupported_depths = df_depths - valid_depths
    unsupported_models = df_models - valid_models
    if unsupported_depths or unsupported_models:
        raise ValueError(
            f"Unsupported scenario values detected: models={sorted(unsupported_models)}, depths={sorted(unsupported_depths)}. "
            f"Aggregate or unsupported scenarios are not permitted."
        )

    for sc in SUPPORTED_PAPER_SCENARIOS:
        if sc.model not in df_models or sc.depth not in df_depths:
            raise ValueError(
                f"Missing required scenario in benchmark data: model='{sc.model}', depth='{sc.depth}'"
            )

    results: Dict[Scenario, ScoringResult] = {}
    for sc in SUPPORTED_PAPER_SCENARIOS:
        res = score_benchmark(
            df=df,
            scenario=sc,
            weights=COMMUNITY_BALANCED_WEIGHTS,
            gates=COMMUNITY_BALANCED_GATES,
            source_data_commit=source_data_commit,
            preset="community_balanced",
        )
        results[sc] = res

    return results


def generate_paper_summary_dataframe(
    df: pd.DataFrame,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> pd.DataFrame:
    """Generate 68-record paper summary dataframe across the 4 supported scenarios.

    Includes canonical criterion scores, overall score, rank, near-tie status,
    raw evidence summaries, warnings, and complete provenance.
    """
    results = generate_paper_composite_results(df, source_data_commit=source_data_commit)
    rows: List[dict] = []

    for sc in SUPPORTED_PAPER_SCENARIOS:
        res = results[sc]
        for r in res.recommendations:
            rows.append({
                "scenario_model": sc.model,
                "scenario_depth": sc.depth,
                "combo": r.combo,
                "rank": r.rank if r.rank is not None else "",
                "overall_score": r.overall_score,
                "display_score": r.display_score if r.display_score is not None else "",
                "score_contiguity": r.score_contiguity,
                "score_accuracy": r.score_accuracy,
                "score_residual": r.score_residual,
                "score_replicon": r.score_replicon,
                "mean_auNGA_ratio": r.mean_auNGA_ratio,
                "mean_error_rate": r.mean_error_rate,
                "mean_mismatches": r.mean_mismatches,
                "mean_indels": r.mean_indels,
                "residual_total_hits": r.residual_total_hits,
                "residual_clean_isolates": r.residual_clean_isolates,
                "residual_affected_isolates": r.residual_affected_isolates,
                "replicon_total_missed": r.replicon_total_missed,
                "replicon_full_missed": r.replicon_full_missed,
                "replicon_partial_missed": r.replicon_partial_missed,
                "replicon_affected_isolates": r.replicon_affected_isolates,
                "mean_duplication_ratio": r.mean_duplication_ratio,
                "total_misassemblies": r.total_misassemblies,
                "is_near_tie": r.is_near_tie,
                "near_tie_label": r.near_tie_label or "",
                "warnings": "; ".join(r.warnings),
                "is_eligible": r.is_eligible,
                "ineligible_reason": r.ineligible_reason or "",
                "analysis_description": ANALYSIS_DESCRIPTION,
                "preset": res.preset or "community_balanced",
                "weight_accuracy": res.weights.accuracy,
                "weight_contiguity": res.weights.contiguity,
                "weight_residual": res.weights.residual,
                "weight_replicon": res.weights.replicon,
                "gate_complete_recovery": res.gates.complete_recovery,
                "gate_zero_residual_hits": res.gates.zero_residual_hits,
                "scoring_version": res.scoring_version,
                "source_data_commit": res.provenance.source_data_commit or "",
                "source_data_hash": res.provenance.source_data_hash,
                "generated_at": res.provenance.generated_at,
            })

    return pd.DataFrame(rows)


def export_paper_summary_csv(
    df: pd.DataFrame,
    out_path: Optional[Union[str, Path]] = None,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> str:
    """Export machine-readable paper summary to CSV string and optionally write to out_path."""
    summary_df = generate_paper_summary_dataframe(df, source_data_commit=source_data_commit)
    csv_str = summary_df.to_csv(index=False)
    if out_path is not None:
        target = Path(out_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(csv_str)
    return csv_str


def plot_paper_composite_figure(
    paper_results: Union[Dict[Scenario, ScoringResult], pd.DataFrame],
    model: str,
    out_path: Union[str, Path],
) -> plt.Figure:
    """Generate a per-depth faceted figure (100x and 20x panels) for the given basecalling model.

    Combinations are ordered independently within each depth panel by rank / overall score descending.
    """
    model_norm = model.lower()
    if model_norm not in ("hac", "sup"):
        raise ValueError(f"Model must be 'hac' or 'sup', got: {model}")

    if isinstance(paper_results, pd.DataFrame):
        results_dict = generate_paper_composite_results(paper_results)
    else:
        results_dict = paper_results

    sns.set_theme(style="whitegrid")

    depth_order = ["100x", "20x"]
    panel_letters = ["A", "B"]
    palette = ["#56b4e9", "#d55e00"]  # skyblue for 100x, vermilion for 20x

    fig, axes = plt.subplots(nrows=1, ncols=len(depth_order), figsize=(20, 7), dpi=300, sharey=True)
    if len(depth_order) == 1:
        axes = [axes]

    for d_idx, depth in enumerate(depth_order):
        ax = axes[d_idx]
        sc = Scenario(model_norm, depth)
        if sc not in results_dict:
            continue
        sc_res = results_dict[sc]
        recs = [r for r in sc_res.recommendations if r.is_eligible]
        # Sort combos independently for this depth panel (best to worst)
        recs.sort(key=lambda r: (-r.overall_score, r.combo))

        combos = [r.combo for r in recs]
        scores = [r.overall_score for r in recs]
        depth_color = palette[d_idx]

        # Plot barplot of preference alignment scores
        bars = ax.bar(
            range(len(combos)),
            scores,
            color=depth_color,
            edgecolor="black",
            linewidth=0.5,
            zorder=2,
        )

        ax.set_title(f"{panel_letters[d_idx]}. Depth: {depth}", fontsize=14, pad=10)
        if d_idx == 0:
            ax.set_ylabel("Preference Alignment Score (0–100)", fontsize=12)
        else:
            ax.set_ylabel("")

        ax.set_xlabel("")
        ax.set_xticks(range(len(combos)))
        ax.set_xticklabels(combos, rotation=45, ha="right", rotation_mode="anchor", fontsize=11)
        ax.set_xlim(-0.5, len(combos) - 0.5)
        ax.set_ylim(0, 105)
        ax.xaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)
        ax.yaxis.grid(True, linestyle="--", color="lightgrey", zorder=0)

    target_path = Path(out_path)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    fig.savefig(target_path, bbox_inches="tight")
    plt.close(fig)
    return fig
