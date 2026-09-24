import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple, Union
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import (
    GateConfig,
    PresetName,
    WeightsConfig,
    get_preset,
)
from qc_scoring.scorer import (
    SCORING_VERSION,
    SOURCE_DATA_BASELINE_COMMIT,
    _calculate_content_hash,
    score_benchmark,
)
from qc_scoring.validation import EXPECTED_SAMPLES

DEFAULT_SCENARIOS: List[Scenario] = [
    Scenario(model="hac", depth="20x"),
    Scenario(model="hac", depth="100x"),
    Scenario(model="sup", depth="20x"),
    Scenario(model="sup", depth="100x"),
]


@dataclass(frozen=True)
class DenominatorRecord:
    scenario_model: str
    scenario_depth: str
    denominator: float
    combo: str
    rank: Optional[int]
    overall_score: float
    score_replicon: float
    is_eligible: bool
    is_near_tie: bool
    is_winner: bool
    scoring_version: str
    source_data_hash: str
    source_data_commit: Optional[str]


@dataclass(frozen=True)
class ScenarioDenominatorSummary:
    scenario: Scenario
    winner_by_denominator: Dict[float, str]
    winner_is_stable: bool
    max_rank_shift: int
    rank_by_combo_and_denominator: Dict[str, Dict[float, Optional[int]]]
    scoring_version: str
    source_data_hash: str
    source_data_commit: Optional[str]


@dataclass
class DenominatorSensitivityResult:
    records: List[DenominatorRecord]
    scenario_summaries: List[ScenarioDenominatorSummary]

    def to_dataframe(self) -> pd.DataFrame:
        rows = [
            {
                "scenario_model": r.scenario_model,
                "scenario_depth": r.scenario_depth,
                "denominator": r.denominator,
                "combo": r.combo,
                "rank": r.rank,
                "overall_score": r.overall_score,
                "score_replicon": r.score_replicon,
                "is_eligible": r.is_eligible,
                "is_near_tie": r.is_near_tie,
                "is_winner": r.is_winner,
                "scoring_version": r.scoring_version,
                "source_data_hash": r.source_data_hash,
                "source_data_commit": r.source_data_commit or "",
            }
            for r in self.records
        ]
        return pd.DataFrame(rows)


def run_denominator_sensitivity(
    df: pd.DataFrame,
    scenarios: Optional[Sequence[Scenario]] = None,
    denominators: Sequence[float] = (2.0, 3.0, 4.0),
    weights: Optional[WeightsConfig] = None,
    gates: Optional[GateConfig] = None,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> DenominatorSensitivityResult:
    if scenarios is None:
        scenarios = DEFAULT_SCENARIOS
    if weights is None:
        weights = get_preset(PresetName.COMMUNITY_BALANCED).weights
    if gates is None:
        gates = get_preset(PresetName.COMMUNITY_BALANCED).gates

    content_hash = _calculate_content_hash(df)
    records: List[DenominatorRecord] = []
    summaries: List[ScenarioDenominatorSummary] = []

    for scenario in scenarios:
        winner_by_denom: Dict[float, str] = {}
        ranks_by_combo: Dict[str, Dict[float, Optional[int]]] = {}

        for denom in denominators:
            scoring_res = score_benchmark(
                df=df,
                scenario=scenario,
                weights=weights,
                gates=gates,
                source_data_commit=source_data_commit,
                replicon_denominator=float(denom),
            )

            # Find winner among eligible recommendations
            winner_combo = ""
            for rec in scoring_res.recommendations:
                if rec.is_eligible and rec.rank == 1:
                    winner_combo = rec.combo
                    break
            winner_by_denom[float(denom)] = winner_combo

            for rec in scoring_res.recommendations:
                if rec.combo not in ranks_by_combo:
                    ranks_by_combo[rec.combo] = {}
                ranks_by_combo[rec.combo][float(denom)] = rec.rank

                is_winner = rec.is_eligible and (rec.rank == 1)
                records.append(
                    DenominatorRecord(
                        scenario_model=scenario.model,
                        scenario_depth=scenario.depth,
                        denominator=float(denom),
                        combo=rec.combo,
                        rank=rec.rank,
                        overall_score=rec.overall_score,
                        score_replicon=rec.score_replicon,
                        is_eligible=rec.is_eligible,
                        is_near_tie=rec.is_near_tie,
                        is_winner=is_winner,
                        scoring_version=SCORING_VERSION,
                        source_data_hash=content_hash,
                        source_data_commit=source_data_commit,
                    )
                )

        # Calculate max rank shift for this scenario
        unique_winners = set(winner_by_denom.values())
        winner_is_stable = (len(unique_winners) == 1) and ("" not in unique_winners)

        max_shift = 0
        for combo, denom_ranks in ranks_by_combo.items():
            valid_ranks = [r for r in denom_ranks.values() if r is not None]
            if len(valid_ranks) > 1:
                shift = max(valid_ranks) - min(valid_ranks)
                if shift > max_shift:
                    max_shift = shift

        summaries.append(
            ScenarioDenominatorSummary(
                scenario=scenario,
                winner_by_denominator=winner_by_denom,
                winner_is_stable=winner_is_stable,
                max_rank_shift=max_shift,
                rank_by_combo_and_denominator=ranks_by_combo,
                scoring_version=SCORING_VERSION,
                source_data_hash=content_hash,
                source_data_commit=source_data_commit,
            )
        )

    return DenominatorSensitivityResult(
        records=records,
        scenario_summaries=summaries,
    )


@dataclass(frozen=True)
class LeaveOneOutRecord:
    scenario_model: str
    scenario_depth: str
    omitted_isolate: str
    combo: str
    rank: Optional[int]
    overall_score: float
    is_eligible: bool
    is_near_tie: bool
    is_winner: bool
    scoring_version: str
    source_data_hash: str
    source_data_commit: Optional[str]


@dataclass(frozen=True)
class ScenarioLeaveOneOutSummary:
    scenario: Scenario
    baseline_winner: str
    winner_by_omitted_isolate: Dict[str, str]
    winner_is_stable: bool
    winner_stability_frequency: float
    max_rank_shift: int
    combo_rank_ranges: Dict[str, Tuple[int, int]]
    scoring_version: str
    source_data_hash: str
    source_data_commit: Optional[str]


@dataclass
class LeaveOneOutResult:
    records: List[LeaveOneOutRecord]
    scenario_summaries: List[ScenarioLeaveOneOutSummary]

    def to_dataframe(self) -> pd.DataFrame:
        rows = [
            {
                "scenario_model": r.scenario_model,
                "scenario_depth": r.scenario_depth,
                "omitted_isolate": r.omitted_isolate,
                "combo": r.combo,
                "rank": r.rank,
                "overall_score": r.overall_score,
                "is_eligible": r.is_eligible,
                "is_near_tie": r.is_near_tie,
                "is_winner": r.is_winner,
                "scoring_version": r.scoring_version,
                "source_data_hash": r.source_data_hash,
                "source_data_commit": r.source_data_commit or "",
            }
            for r in self.records
        ]
        return pd.DataFrame(rows)


def run_leave_one_isolate_out(
    df: pd.DataFrame,
    scenarios: Optional[Sequence[Scenario]] = None,
    isolates: Optional[Sequence[str]] = None,
    weights: Optional[WeightsConfig] = None,
    gates: Optional[GateConfig] = None,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> LeaveOneOutResult:
    if scenarios is None:
        scenarios = DEFAULT_SCENARIOS
    if isolates is None:
        isolates = sorted(EXPECTED_SAMPLES)
    if weights is None:
        weights = get_preset(PresetName.COMMUNITY_BALANCED).weights
    if gates is None:
        gates = get_preset(PresetName.COMMUNITY_BALANCED).gates

    content_hash = _calculate_content_hash(df)
    records: List[LeaveOneOutRecord] = []
    summaries: List[ScenarioLeaveOneOutSummary] = []

    for scenario in scenarios:
        # First compute baseline with all 13 isolates
        base_res = score_benchmark(
            df=df,
            scenario=scenario,
            weights=weights,
            gates=gates,
            source_data_commit=source_data_commit,
        )
        base_winner = ""
        base_ranks: Dict[str, Optional[int]] = {}
        for rec in base_res.recommendations:
            base_ranks[rec.combo] = rec.rank
            if rec.is_eligible and rec.rank == 1 and not base_winner:
                base_winner = rec.combo

        winner_by_isolate: Dict[str, str] = {}
        combo_ranks: Dict[str, List[int]] = {combo: [] for combo in base_ranks}

        for omitted in isolates:
            subset_samples = EXPECTED_SAMPLES - {omitted}
            filtered_df = df[df["sample"] != omitted].copy()

            res = score_benchmark(
                df=filtered_df,
                scenario=scenario,
                weights=weights,
                gates=gates,
                source_data_commit=source_data_commit,
                expected_samples=subset_samples,
            )

            current_winner = ""
            for rec in res.recommendations:
                if rec.is_eligible and rec.rank == 1 and not current_winner:
                    current_winner = rec.combo
                if rec.rank is not None:
                    combo_ranks[rec.combo].append(rec.rank)

                records.append(
                    LeaveOneOutRecord(
                        scenario_model=scenario.model,
                        scenario_depth=scenario.depth,
                        omitted_isolate=omitted,
                        combo=rec.combo,
                        rank=rec.rank,
                        overall_score=rec.overall_score,
                        is_eligible=rec.is_eligible,
                        is_near_tie=rec.is_near_tie,
                        is_winner=(rec.is_eligible and rec.rank == 1),
                        scoring_version=SCORING_VERSION,
                        source_data_hash=content_hash,
                        source_data_commit=source_data_commit,
                    )
                )

            winner_by_isolate[omitted] = current_winner

        # Stability calculations
        matches_base = sum(1 for w in winner_by_isolate.values() if w == base_winner)
        freq = matches_base / len(isolates) if isolates else 0.0
        is_stable = (matches_base == len(isolates))

        max_shift = 0
        rank_ranges: Dict[str, Tuple[int, int]] = {}
        for combo, ranks in combo_ranks.items():
            if ranks:
                min_r = min(ranks)
                max_r = max(ranks)
                rank_ranges[combo] = (min_r, max_r)
                shift = max_r - min_r
                if shift > max_shift:
                    max_shift = shift
            else:
                rank_ranges[combo] = (0, 0)

        summaries.append(
            ScenarioLeaveOneOutSummary(
                scenario=scenario,
                baseline_winner=base_winner,
                winner_by_omitted_isolate=winner_by_isolate,
                winner_is_stable=is_stable,
                winner_stability_frequency=freq,
                max_rank_shift=max_shift,
                combo_rank_ranges=rank_ranges,
                scoring_version=SCORING_VERSION,
                source_data_hash=content_hash,
                source_data_commit=source_data_commit,
            )
        )

    return LeaveOneOutResult(
        records=records,
        scenario_summaries=summaries,
    )


@dataclass(frozen=True)
class RespondentProfileRecord:
    scenario_model: str
    scenario_depth: str
    respondent_id: int
    timestamp: str
    weight_accuracy: float
    weight_contiguity: float
    weight_residual: float
    weight_replicon: float
    combo: str
    rank: Optional[int]
    overall_score: float
    is_eligible: bool
    is_near_tie: bool
    is_winner: bool
    scoring_version: str
    source_data_hash: str
    survey_data_hash: str
    source_data_commit: Optional[str]


@dataclass(frozen=True)
class ScenarioRespondentSummary:
    scenario: Scenario
    total_respondents: int
    winner_counts: Dict[str, int]
    winner_frequencies: Dict[str, float]
    combo_rank_variation: Dict[str, Dict[str, float]]
    scoring_version: str
    source_data_hash: str
    survey_data_hash: str
    source_data_commit: Optional[str]


@dataclass
class RespondentProfilesResult:
    records: List[RespondentProfileRecord]
    scenario_summaries: List[ScenarioRespondentSummary]
    survey_provenance_hash: str

    def to_dataframe(self) -> pd.DataFrame:
        rows = [
            {
                "scenario_model": r.scenario_model,
                "scenario_depth": r.scenario_depth,
                "respondent_id": r.respondent_id,
                "timestamp": r.timestamp,
                "weight_accuracy": r.weight_accuracy,
                "weight_contiguity": r.weight_contiguity,
                "weight_residual": r.weight_residual,
                "weight_replicon": r.weight_replicon,
                "combo": r.combo,
                "rank": r.rank,
                "overall_score": r.overall_score,
                "is_eligible": r.is_eligible,
                "is_near_tie": r.is_near_tie,
                "is_winner": r.is_winner,
                "scoring_version": r.scoring_version,
                "source_data_hash": r.source_data_hash,
                "survey_data_hash": r.survey_data_hash,
                "source_data_commit": r.source_data_commit or "",
            }
            for r in self.records
        ]
        return pd.DataFrame(rows)


def run_respondent_profiles(
    df: pd.DataFrame,
    survey_df: pd.DataFrame,
    scenarios: Optional[Sequence[Scenario]] = None,
    gates: Optional[GateConfig] = None,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> RespondentProfilesResult:
    if scenarios is None:
        scenarios = DEFAULT_SCENARIOS
    if gates is None:
        gates = get_preset(PresetName.COMMUNITY_BALANCED).gates

    required_survey_cols = ["w_accuracy", "w_contiguity", "w_decontam", "w_replicon"]
    for col in required_survey_cols:
        if col not in survey_df.columns:
            raise ValueError(f"Missing required survey column: '{col}'")

    content_hash = _calculate_content_hash(df)
    survey_hash = _calculate_content_hash(survey_df)

    records: List[RespondentProfileRecord] = []
    summaries: List[ScenarioRespondentSummary] = []

    for scenario in scenarios:
        winner_counts: Dict[str, int] = {}
        combo_ranks: Dict[str, List[int]] = {}

        for respondent_id, row in enumerate(survey_df.to_dict(orient="records"), start=1):
            ts = str(row.get("timestamp", ""))
            w_acc = float(row["w_accuracy"])
            w_cont = float(row["w_contiguity"])
            w_res = float(row["w_decontam"])
            w_rep = float(row["w_replicon"])

            weights = WeightsConfig(
                accuracy=w_acc,
                contiguity=w_cont,
                residual=w_res,
                replicon=w_rep,
            )

            res = score_benchmark(
                df=df,
                scenario=scenario,
                weights=weights,
                gates=gates,
                source_data_commit=source_data_commit,
            )

            resp_winner = ""
            for rec in res.recommendations:
                if rec.combo not in combo_ranks:
                    combo_ranks[rec.combo] = []
                if rec.rank is not None:
                    combo_ranks[rec.combo].append(rec.rank)

                if rec.is_eligible and rec.rank == 1 and not resp_winner:
                    resp_winner = rec.combo

                is_winner = rec.is_eligible and (rec.rank == 1)
                records.append(
                    RespondentProfileRecord(
                        scenario_model=scenario.model,
                        scenario_depth=scenario.depth,
                        respondent_id=respondent_id,
                        timestamp=ts,
                        weight_accuracy=w_acc,
                        weight_contiguity=w_cont,
                        weight_residual=w_res,
                        weight_replicon=w_rep,
                        combo=rec.combo,
                        rank=rec.rank,
                        overall_score=rec.overall_score,
                        is_eligible=rec.is_eligible,
                        is_near_tie=rec.is_near_tie,
                        is_winner=is_winner,
                        scoring_version=SCORING_VERSION,
                        source_data_hash=content_hash,
                        survey_data_hash=survey_hash,
                        source_data_commit=source_data_commit,
                    )
                )

            if resp_winner:
                winner_counts[resp_winner] = winner_counts.get(resp_winner, 0) + 1

        n_resp = len(survey_df)
        winner_frequencies = {
            combo: (count / n_resp) for combo, count in winner_counts.items()
        }

        combo_variation: Dict[str, Dict[str, float]] = {}
        for combo, ranks in combo_ranks.items():
            if ranks:
                combo_variation[combo] = {
                    "min_rank": float(min(ranks)),
                    "max_rank": float(max(ranks)),
                    "mean_rank": float(sum(ranks) / len(ranks)),
                }
            else:
                combo_variation[combo] = {
                    "min_rank": 0.0,
                    "max_rank": 0.0,
                    "mean_rank": 0.0,
                }

        summaries.append(
            ScenarioRespondentSummary(
                scenario=scenario,
                total_respondents=n_resp,
                winner_counts=winner_counts,
                winner_frequencies=winner_frequencies,
                combo_rank_variation=combo_variation,
                scoring_version=SCORING_VERSION,
                source_data_hash=content_hash,
                survey_data_hash=survey_hash,
                source_data_commit=source_data_commit,
            )
        )

    return RespondentProfilesResult(
        records=records,
        scenario_summaries=summaries,
        survey_provenance_hash=survey_hash,
    )


@dataclass
class RobustnessEvidenceBundle:
    denominator_result: DenominatorSensitivityResult
    leave_one_out_result: LeaveOneOutResult
    respondent_profiles_result: RespondentProfilesResult
    output_dir: Path
    summary_markdown: str
    summary_json: Dict[str, Any]


def _format_scenario_name(s: Scenario) -> str:
    return f"{s.model.upper()} {s.depth}"


def _generate_markdown_summary(
    denom_res: DenominatorSensitivityResult,
    loo_res: LeaveOneOutResult,
    resp_res: RespondentProfilesResult,
    benchmark_hash: str,
    survey_hash: str,
    commit: Optional[str],
) -> str:
    lines = [
        "# Scoring v1 Robustness Evidence Summary",
        "",
        "## Provenance & Metadata",
        f"- **Scoring Version**: {SCORING_VERSION}",
        f"- **Benchmark Source Commit**: `{commit or 'N/A'}`",
        f"- **Benchmark Source Data Hash (SHA-256)**: `{benchmark_hash}`",
        f"- **Survey Data Hash (SHA-256)**: `{survey_hash}`",
        f"- **Supported Scenarios**: HAC 20×, HAC 100×, SUP 20×, SUP 100×",
        "",
        "## Methodological Scope and Boundary Disclosures",
        "- **Descriptive Non-Inferential Scope**: Respondent-specific results are explicitly descriptive of the 22 observed responses in the microbial QC survey. They make no inference that these 22 respondents represent a wider population of bioinformaticians or microbiologists.",
        "- **No Re-survey**: No respondents were re-contacted and no new survey was required.",
        "- **Superseded Formulations Excluded**: This analysis does not compare against the superseded dynamic min–max or geometric mean scoring systems.",
        "- **Fixed Warning Threshold**: The accepted 25-point isolate-variability warning threshold is preserved unchanged.",
        "- **Outcome Interpretation**: A changed winner under alternative denominators, sub-cohorts, or preferences is reported as an empirical robustness result, not treated as a failed test.",
        "",
        "## 1. Replicon Calibration Denominator Sensitivity ($B \\in \\{2, 3, 4\\}$)",
        "Rankings were recomputed with replicon-loss denominators $B = 2.0, 3.0, 4.0$ while holding every other scoring-v1 decision fixed under Community-balanced weights (accuracy: 28%, contiguity: 20%, residual: 17%, replicon: 35%).",
        "",
        "| Scenario | Winner ($B=2$) | Winner ($B=3$, Canonical) | Winner ($B=4$) | Winner Stable? | Max Rank Shift |",
        "| :--- | :--- | :--- | :--- | :---: | :---: |",
    ]

    for s_denom in denom_res.scenario_summaries:
        sc_label = _format_scenario_name(s_denom.scenario)
        w2 = s_denom.winner_by_denominator.get(2.0, "N/A")
        w3 = s_denom.winner_by_denominator.get(3.0, "N/A")
        w4 = s_denom.winner_by_denominator.get(4.0, "N/A")
        stable_str = "Yes" if s_denom.winner_is_stable else "No"
        lines.append(f"| {sc_label} | `{w2}` | `{w3}` | `{w4}` | {stable_str} | {s_denom.max_rank_shift} |")

    lines.extend([
        "",
        "## 2. Leave-One-Isolate-Out Stability (13 Isolates)",
        "Rankings were recomputed leaving out one of the 13 reference isolates at a time for all 4 supported scenarios (52 sub-cohort evaluations total).",
        "",
        "| Scenario | Canonical Winner | LOO Winner Match Frequency | Winner Stable Across All 13? | Max Rank Shift Across Combinations |",
        "| :--- | :--- | :---: | :---: | :---: |",
    ])

    for s_loo in loo_res.scenario_summaries:
        sc_label = _format_scenario_name(s_loo.scenario)
        freq_pct = s_loo.winner_stability_frequency * 100.0
        n_match = len([w for w in s_loo.winner_by_omitted_isolate.values() if w == s_loo.baseline_winner])
        stable_str = "Yes" if s_loo.winner_is_stable else "No"
        lines.append(
            f"| {sc_label} | `{s_loo.baseline_winner}` | {n_match}/13 ({freq_pct:.1f}%) | {stable_str} | {s_loo.max_rank_shift} |"
        )

    lines.extend([
        "",
        "## 3. Observed Respondent Weight Profiles (22 Respondents)",
        "The 22 observed respondent preference profiles from `microbial-qc-survey.csv` were applied directly to the canonical cohort criterion scores.",
        "",
        "| Scenario | Leading Combination(s) | Winner Frequency among Observed Respondents |",
        "| :--- | :--- | :--- |",
    ])

    for s_resp in resp_res.scenario_summaries:
        sc_label = _format_scenario_name(s_resp.scenario)
        parts = []
        for combo, count in sorted(s_resp.winner_counts.items(), key=lambda item: -item[1]):
            freq = s_resp.winner_frequencies[combo] * 100.0
            parts.append(f"`{combo}`: {count}/22 ({freq:.1f}%)")
        lines.append(f"| {sc_label} | {', '.join(parts)} |")

    lines.extend([
        "",
        "### Top Combinations Rank Variation Across Observed Profiles",
        "",
        "| Scenario | Combination | Min Rank | Max Rank | Mean Rank |",
        "| :--- | :--- | :---: | :---: | :---: |",
    ])

    for s_resp in resp_res.scenario_summaries:
        sc_label = _format_scenario_name(s_resp.scenario)
        sorted_combos = sorted(s_resp.combo_rank_variation.items(), key=lambda item: item[1]["mean_rank"])
        # Display top 5 combos per scenario
        for combo, var in sorted_combos[:5]:
            lines.append(
                f"| {sc_label} | `{combo}` | {int(var['min_rank'])} | {int(var['max_rank'])} | {var['mean_rank']:.2f} |"
            )

    lines.extend([
        "",
        "## Artifact Manifest",
        "- `denominator_sensitivity.csv`: 204 records covering all 4 scenarios × 3 denominators × 17 combinations.",
        "- `leave_one_out.csv`: 884 records covering all 4 scenarios × 13 isolate omissions × 17 combinations.",
        "- `respondent_profiles.csv`: 1,496 records covering all 4 scenarios × 22 observed survey respondents × 17 combinations.",
        "- `evidence_summary.json`: Machine-readable structured provenance and stability summaries.",
        "- `summary.md`: Human-readable summary report (this file).",
    ])

    return "\n".join(lines) + "\n"


def _generate_json_summary(
    denom_res: DenominatorSensitivityResult,
    loo_res: LeaveOneOutResult,
    resp_res: RespondentProfilesResult,
    benchmark_hash: str,
    survey_hash: str,
    commit: Optional[str],
) -> Dict[str, Any]:
    return {
        "scoring_version": SCORING_VERSION,
        "source_data_commit": commit,
        "source_data_hash": benchmark_hash,
        "survey_data_hash": survey_hash,
        "denominator_sensitivity": [
            {
                "scenario_model": s.scenario.model,
                "scenario_depth": s.scenario.depth,
                "winner_by_denominator": {str(k): v for k, v in s.winner_by_denominator.items()},
                "winner_is_stable": s.winner_is_stable,
                "max_rank_shift": s.max_rank_shift,
                "ranks_by_combo": s.rank_by_combo_and_denominator,
            }
            for s in denom_res.scenario_summaries
        ],
        "leave_one_out": [
            {
                "scenario_model": s.scenario.model,
                "scenario_depth": s.scenario.depth,
                "baseline_winner": s.baseline_winner,
                "winner_by_omitted_isolate": s.winner_by_omitted_isolate,
                "winner_is_stable": s.winner_is_stable,
                "winner_stability_frequency": s.winner_stability_frequency,
                "max_rank_shift": s.max_rank_shift,
                "combo_rank_ranges": {k: list(v) for k, v in s.combo_rank_ranges.items()},
            }
            for s in loo_res.scenario_summaries
        ],
        "respondent_profiles": [
            {
                "scenario_model": s.scenario.model,
                "scenario_depth": s.scenario.depth,
                "total_respondents": s.total_respondents,
                "winner_counts": s.winner_counts,
                "winner_frequencies": s.winner_frequencies,
                "combo_rank_variation": s.combo_rank_variation,
            }
            for s in resp_res.scenario_summaries
        ],
    }


def generate_robustness_evidence(
    benchmark_path: Union[str, Path] = "assets/data/assembly_metrics.csv",
    survey_path: Union[str, Path] = "microbial-qc-survey.csv",
    output_dir: Union[str, Path] = "results/robustness",
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
) -> RobustnessEvidenceBundle:
    df = pd.read_csv(benchmark_path)
    survey_df = pd.read_csv(survey_path)

    benchmark_hash = _calculate_content_hash(df)
    survey_hash = _calculate_content_hash(survey_df)

    denom_res = run_denominator_sensitivity(
        df=df,
        source_data_commit=source_data_commit,
    )
    loo_res = run_leave_one_isolate_out(
        df=df,
        source_data_commit=source_data_commit,
    )
    resp_res = run_respondent_profiles(
        df=df,
        survey_df=survey_df,
        source_data_commit=source_data_commit,
    )

    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    denom_res.to_dataframe().to_csv(out_path / "denominator_sensitivity.csv", index=False)
    loo_res.to_dataframe().to_csv(out_path / "leave_one_out.csv", index=False)
    resp_res.to_dataframe().to_csv(out_path / "respondent_profiles.csv", index=False)

    summary_md = _generate_markdown_summary(
        denom_res=denom_res,
        loo_res=loo_res,
        resp_res=resp_res,
        benchmark_hash=benchmark_hash,
        survey_hash=survey_hash,
        commit=source_data_commit,
    )
    (out_path / "summary.md").write_text(summary_md, encoding="utf-8")

    summary_json = _generate_json_summary(
        denom_res=denom_res,
        loo_res=loo_res,
        resp_res=resp_res,
        benchmark_hash=benchmark_hash,
        survey_hash=survey_hash,
        commit=source_data_commit,
    )
    (out_path / "evidence_summary.json").write_text(
        json.dumps(summary_json, indent=2), encoding="utf-8"
    )

    return RobustnessEvidenceBundle(
        denominator_result=denom_res,
        leave_one_out_result=loo_res,
        respondent_profiles_result=resp_res,
        output_dir=out_path,
        summary_markdown=summary_md,
        summary_json=summary_json,
    )


if __name__ == "__main__":
    import sys
    bench = sys.argv[1] if len(sys.argv) > 1 else "assets/data/assembly_metrics.csv"
    surv = sys.argv[2] if len(sys.argv) > 2 else "microbial-qc-survey.csv"
    out = sys.argv[3] if len(sys.argv) > 3 else "results/robustness"
    generate_robustness_evidence(benchmark_path=bench, survey_path=surv, output_dir=out)
    print(f"Robustness evidence generated successfully in {out}")



