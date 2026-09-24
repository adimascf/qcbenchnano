from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import io
from typing import Dict, List, Optional, Set
import pandas as pd
from qc_scoring.criteria import calculate_cohort_criteria, CohortCriteriaResult
from qc_scoring.models import Scenario
from qc_scoring.preferences import GateConfig, WeightsConfig, validate_weights
from qc_scoring.validation import filter_scenario_data

SCORING_VERSION = "1.0"
SOURCE_DATA_BASELINE_COMMIT = "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"


@dataclass
class Provenance:
    scoring_version: str
    source_data_hash: str
    source_data_commit: Optional[str]
    scenario_model: str
    scenario_depth: str
    weights: Dict[str, float]
    gates: Dict[str, bool]
    generated_at: str
    preset: Optional[str] = None
    replicon_denominator: float = 3.0


@dataclass
class CombinationRecommendation:
    combo: str
    rank: Optional[int]
    is_eligible: bool
    ineligible_reason: Optional[str]
    is_near_tie: bool
    overall_score: float
    display_score: Optional[float]
    cohort: CohortCriteriaResult

    warnings: List[str] = field(default_factory=list)
    is_insufficient_data: bool = False

    # Delegation properties for convenience
    @property
    def score_contiguity(self) -> float:
        return self.cohort.score_contiguity

    @property
    def score_accuracy(self) -> float:
        return self.cohort.score_accuracy

    @property
    def score_residual(self) -> float:
        return self.cohort.score_residual

    @property
    def score_replicon(self) -> float:
        return self.cohort.score_replicon

    @property
    def mean_auNGA_ratio(self) -> float:
        return self.cohort.mean_auNGA_ratio

    @property
    def mean_error_rate(self) -> float:
        return self.cohort.mean_error_rate

    @property
    def mean_mismatches(self) -> float:
        return self.cohort.mean_mismatches

    @property
    def mean_indels(self) -> float:
        return self.cohort.mean_indels

    @property
    def residual_total_hits(self) -> int:
        return self.cohort.residual_total_hits

    @property
    def residual_clean_isolates(self) -> int:
        return self.cohort.residual_clean_isolates

    @property
    def residual_affected_isolates(self) -> int:
        return self.cohort.residual_affected_isolates

    @property
    def replicon_total_missed(self) -> int:
        return self.cohort.replicon_total_missed

    @property
    def replicon_full_missed(self) -> int:
        return self.cohort.replicon_full_missed

    @property
    def replicon_partial_missed(self) -> int:
        return self.cohort.replicon_partial_missed

    @property
    def replicon_affected_isolates(self) -> int:
        return self.cohort.replicon_affected_isolates

    @property
    def isolate_contiguity_scores(self) -> Dict[str, float]:
        return self.cohort.isolate_contiguity_scores

    @property
    def isolate_accuracy_scores(self) -> Dict[str, float]:
        return self.cohort.isolate_accuracy_scores

    @property
    def mean_duplication_ratio(self) -> float:
        return self.cohort.mean_duplication_ratio

    @property
    def total_misassemblies(self) -> int:
        return self.cohort.total_misassemblies

    @property
    def isolate_duplication_ratios(self) -> Dict[str, float]:
        return self.cohort.isolate_duplication_ratios

    @property
    def isolate_misassemblies(self) -> Dict[str, int]:
        return self.cohort.isolate_misassemblies

    @property
    def near_tie_label(self) -> Optional[str]:
        return "similar overall scores" if self.is_near_tie else None


@dataclass
class ScoringResult:
    scoring_version: str
    scenario: Scenario
    weights: WeightsConfig
    gates: GateConfig
    recommendations: List[CombinationRecommendation]
    provenance: Provenance
    preset: Optional[str] = None

    @property
    def eligible_count(self) -> int:
        return sum(1 for r in self.recommendations if r.is_eligible)

    @property
    def ineligible_count(self) -> int:
        return sum(1 for r in self.recommendations if not r.is_eligible)

    @property
    def leading_recommendation(self) -> Optional[CombinationRecommendation]:
        eligible = [r for r in self.recommendations if r.is_eligible]
        return eligible[0] if eligible else None


def _calculate_content_hash(df: pd.DataFrame) -> str:
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    return hashlib.sha256(csv_bytes).hexdigest()


def score_benchmark(
    df: pd.DataFrame,
    scenario: Scenario,
    weights: WeightsConfig,
    gates: Optional[GateConfig] = None,
    source_data_commit: Optional[str] = SOURCE_DATA_BASELINE_COMMIT,
    preset: Optional[str] = None,
    replicon_denominator: float = 3.0,
    expected_samples: Optional[Set[str]] = None,
) -> ScoringResult:
    if gates is None:
        gates = GateConfig()

    validate_weights(weights)
    filtered_df, incomplete_reasons = filter_scenario_data(
        df, scenario, expected_samples=expected_samples
    )

    if preset is None:
        from qc_scoring.preferences import PRESETS
        for p_name, p_def in PRESETS.items():
            if (
                abs(weights.accuracy - p_def.weights.accuracy) < 1e-4
                and abs(weights.contiguity - p_def.weights.contiguity) < 1e-4
                and abs(weights.residual - p_def.weights.residual) < 1e-4
                and abs(weights.replicon - p_def.weights.replicon) < 1e-4
                and gates.complete_recovery == p_def.gates.complete_recovery
                and gates.zero_residual_hits == p_def.gates.zero_residual_hits
            ):
                preset = p_name.value
                break
        if preset is None:
            preset = "custom"

    eligible_recs: List[CombinationRecommendation] = []
    excluded_recs: List[CombinationRecommendation] = []
    insufficient_recs: List[CombinationRecommendation] = []

    # Process complete combinations
    grouped = filtered_df.groupby("combo", observed=True)
    for combo, combo_df in grouped:
        combo_name = str(combo)
        cohort: CohortCriteriaResult = calculate_cohort_criteria(
            combo_df, replicon_denominator=replicon_denominator
        )

        # Preference alignment score (unrounded full precision)
        overall = (
            (weights.accuracy * cohort.score_accuracy)
            + (weights.contiguity * cohort.score_contiguity)
            + (weights.residual * cohort.score_residual)
            + (weights.replicon * cohort.score_replicon)
        ) / 100.0

        # Evaluate eligibility gates
        is_eligible = True
        fail_reasons: List[str] = []

        if gates.complete_recovery and not cohort.all_replicons_complete:
            is_eligible = False
            fail_reasons.append("fails complete-recovery gate (<95% coverage on at least one replicon)")

        if gates.zero_residual_hits and not cohort.zero_residual_hits:
            is_eligible = False
            fail_reasons.append(
                f"fails zero-residual-hits gate ({cohort.residual_total_hits} residual hit(s) detected)"
            )

        # Generate concrete warnings
        warnings: List[str] = []
        if cohort.residual_total_hits > 0:
            warnings.append(
                f"{cohort.residual_total_hits} residual adapter/barcode hit(s) across "
                f"{cohort.residual_affected_isolates} isolate(s)"
            )

        if cohort.replicon_total_missed > 0:
            warnings.append(
                f"{cohort.replicon_total_missed} missed or severely incomplete replicon(s) across "
                f"{cohort.replicon_affected_isolates} isolate(s) "
                f"({cohort.replicon_full_missed} full, {cohort.replicon_partial_missed} partial)"
            )

        # Isolate variability warnings (>= 25 points between mean and min isolate)
        min_acc = min(cohort.isolate_accuracy_scores.values()) if cohort.isolate_accuracy_scores else cohort.score_accuracy
        if (cohort.score_accuracy - min_acc) >= 25.0:
            warnings.append(
                f"Variable sequence accuracy across isolates: lowest isolate score is {min_acc:.1f} "
                f"(cohort mean {cohort.score_accuracy:.1f})"
            )

        min_cont = min(cohort.isolate_contiguity_scores.values()) if cohort.isolate_contiguity_scores else cohort.score_contiguity
        if (cohort.score_contiguity - min_cont) >= 25.0:
            warnings.append(
                f"Variable contiguity across isolates: lowest isolate score is {min_cont:.1f} "
                f"(cohort mean {cohort.score_contiguity:.1f})"
            )

        rec = CombinationRecommendation(
            combo=combo_name,
            rank=None,
            is_eligible=is_eligible,
            ineligible_reason="; ".join(fail_reasons) if fail_reasons else None,
            is_near_tie=False,
            overall_score=overall,
            display_score=round(overall, 1),
            cohort=cohort,
            warnings=warnings,
        )

        if is_eligible:
            eligible_recs.append(rec)
        else:
            excluded_recs.append(rec)

    # Process combinations with insufficient data
    empty_cohort = CohortCriteriaResult(
        score_contiguity=0.0,
        score_accuracy=0.0,
        score_residual=0.0,
        score_replicon=0.0,
        mean_auNGA_ratio=0.0,
        mean_error_rate=0.0,
        mean_mismatches=0.0,
        mean_indels=0.0,
        residual_total_hits=0,
        residual_clean_isolates=0,
        residual_affected_isolates=0,
        replicon_total_missed=0,
        replicon_full_missed=0,
        replicon_partial_missed=0,
        replicon_affected_isolates=0,
        all_replicons_complete=False,
        zero_residual_hits=False,
    )
    for combo_name, reason in incomplete_reasons.items():
        rec = CombinationRecommendation(
            combo=combo_name,
            rank=None,
            is_eligible=False,
            ineligible_reason=reason,
            is_near_tie=False,
            overall_score=0.0,
            display_score=None,
            cohort=empty_cohort,
            warnings=["Insufficient benchmark data"],
            is_insufficient_data=True,
        )
        insufficient_recs.append(rec)

    # Rank eligible combinations: higher overall_score first, stable alphabetical combo sort
    eligible_recs.sort(key=lambda r: (-r.overall_score, r.combo))

    # Competition ranking for exact ties (1, 1, 3, etc.)
    for idx, rec in enumerate(eligible_recs):
        if idx == 0:
            rec.rank = 1
        else:
            prev = eligible_recs[idx - 1]
            if abs(rec.overall_score - prev.overall_score) < 1e-9:
                rec.rank = prev.rank
            else:
                rec.rank = idx + 1

    # Near-tie calculation (< 1.0 point difference between adjacent eligible combinations)
    for idx in range(len(eligible_recs) - 1):
        diff = abs(eligible_recs[idx].overall_score - eligible_recs[idx + 1].overall_score)
        if diff < 1.0:
            eligible_recs[idx].is_near_tie = True
            eligible_recs[idx + 1].is_near_tie = True

    # Order excluded combinations by overall score descending (then combo name)
    excluded_recs.sort(key=lambda r: (-r.overall_score, r.combo))

    # Order insufficient combinations alphabetically by combo name
    insufficient_recs.sort(key=lambda r: r.combo)

    all_recommendations = eligible_recs + excluded_recs + insufficient_recs

    now_iso = datetime.now(timezone.utc).isoformat()
    content_hash = _calculate_content_hash(df)

    provenance = Provenance(
        scoring_version=SCORING_VERSION,
        source_data_hash=content_hash,
        source_data_commit=source_data_commit,
        scenario_model=scenario.model,
        scenario_depth=scenario.depth,
        preset=preset,
        weights=weights.as_dict(),
        gates=gates.as_dict(),
        generated_at=now_iso,
        replicon_denominator=replicon_denominator,
    )

    return ScoringResult(
        scoring_version=SCORING_VERSION,
        scenario=scenario,
        weights=weights,
        gates=gates,
        recommendations=all_recommendations,
        provenance=provenance,
        preset=preset,
    )


def recommendations_to_dataframe(result: ScoringResult) -> pd.DataFrame:
    rows = []
    for r in result.recommendations:
        rows.append({
            "rank": r.rank if r.rank is not None else "",
            "combo": r.combo,
            "overall_score": r.overall_score,
            "display_score": r.display_score if r.display_score is not None else "",
            "is_eligible": r.is_eligible,
            "is_insufficient_data": r.is_insufficient_data,
            "ineligible_reason": r.ineligible_reason or "",
            "is_near_tie": r.is_near_tie,
            "near_tie_label": r.near_tie_label or "",
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
            "warnings": "; ".join(r.warnings),
            "preset": result.preset or "custom",
            "scenario_model": result.scenario.model,
            "scenario_depth": result.scenario.depth,
            "weight_accuracy": result.weights.accuracy,
            "weight_contiguity": result.weights.contiguity,
            "weight_residual": result.weights.residual,
            "weight_replicon": result.weights.replicon,
            "gate_complete_recovery": result.gates.complete_recovery,
            "gate_zero_residual_hits": result.gates.zero_residual_hits,
            "scoring_version": result.scoring_version,
            "source_data_commit": result.provenance.source_data_commit or "",
            "source_data_hash": result.provenance.source_data_hash,
            "generated_at": result.provenance.generated_at,
        })
    return pd.DataFrame(rows)


def export_recommendations_csv(result: ScoringResult) -> str:
    df = recommendations_to_dataframe(result)
    buffer = io.StringIO()
    df.to_csv(buffer, index=False)
    return buffer.getvalue()
