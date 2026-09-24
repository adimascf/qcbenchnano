from dataclasses import dataclass, field
from typing import Dict, List
import numpy as np
import pandas as pd
from qc_scoring.parser import validate_and_parse_replicon_row

def score_isolate_contiguity(auNGA_ratio: float) -> float:
    return float(max(0.0, 1.0 - abs(auNGA_ratio - 1.0)) * 100.0)

def score_isolate_accuracy(mismatches_per_100kbp: float, indels_per_100kbp: float) -> float:
    event_rate = mismatches_per_100kbp + indels_per_100kbp
    ratio = 1.0 - (event_rate / 10.0)
    clipped = min(1.0, max(0.0, ratio))
    return float(clipped * 100.0)

def score_isolate_residual_clean(contamination_count: int | float) -> float:
    return 100.0 if contamination_count == 0 else 0.0

def score_replicon_recovery(total_missed_cohort: int | float, denominator: float = 3.0) -> float:
    if denominator <= 0:
        raise ValueError(f"Replicon calibration denominator must be positive, got {denominator}")
    ratio = 1.0 - (total_missed_cohort / denominator)
    return float(max(0.0, ratio) * 100.0)


@dataclass(frozen=True)
class CohortCriteriaResult:
    # 4 fixed 0-100 criterion scores
    score_contiguity: float
    score_accuracy: float
    score_residual: float
    score_replicon: float

    # Raw evidence summaries
    mean_auNGA_ratio: float
    mean_error_rate: float
    mean_mismatches: float
    mean_indels: float

    residual_total_hits: int
    residual_clean_isolates: int
    residual_affected_isolates: int

    replicon_total_missed: int
    replicon_full_missed: int
    replicon_partial_missed: int
    replicon_affected_isolates: int

    # Gate evaluations
    all_replicons_complete: bool
    zero_residual_hits: bool

    # Supporting diagnostics
    mean_duplication_ratio: float = 1.0
    total_misassemblies: int = 0

    # Per-isolate details
    isolate_contiguity_scores: Dict[str, float] = field(default_factory=dict)
    isolate_accuracy_scores: Dict[str, float] = field(default_factory=dict)
    isolate_duplication_ratios: Dict[str, float] = field(default_factory=dict)
    isolate_misassemblies: Dict[str, int] = field(default_factory=dict)

def calculate_cohort_criteria(combo_df: pd.DataFrame, replicon_denominator: float = 3.0) -> CohortCriteriaResult:
    n_isolates = len(combo_df)
    if n_isolates == 0:
        raise ValueError("Cannot calculate cohort criteria on empty DataFrame")

    contiguity_scores: Dict[str, float] = {}
    accuracy_scores: Dict[str, float] = {}
    duplication_ratios: Dict[str, float] = {}
    misassemblies_counts: Dict[str, int] = {}
    clean_isolates_count = 0
    total_residual_hits = 0
    affected_residual_count = 0

    full_missed_sum = 0
    partial_missed_sum = 0
    total_missed_sum = 0
    affected_replicon_count = 0

    all_complete_recovery = True

    mismatches_sum = 0.0
    indels_sum = 0.0
    aunga_sum = 0.0
    duplication_sum = 0.0
    misassemblies_sum = 0

    for _, row in combo_df.iterrows():
        sample_name = str(row["sample"])

        # Contiguity
        aunga_ratio = float(row["auNGA_ratio"])
        aunga_sum += aunga_ratio
        cont_score = score_isolate_contiguity(aunga_ratio)
        contiguity_scores[sample_name] = cont_score

        # Accuracy
        mismatches = float(row["Mismatches per 100kbp"])
        indels = float(row["Indels per 100kbp"])
        mismatches_sum += mismatches
        indels_sum += indels
        acc_score = score_isolate_accuracy(mismatches, indels)
        accuracy_scores[sample_name] = acc_score

        # Supporting diagnostics
        dup_val = float(row["Duplication_ratio"]) if "Duplication_ratio" in row and pd.notna(row["Duplication_ratio"]) else 1.0
        duplication_sum += dup_val
        duplication_ratios[sample_name] = dup_val

        mis_val = int(row["misassemblies"]) if "misassemblies" in row and pd.notna(row["misassemblies"]) else 0
        misassemblies_sum += mis_val
        misassemblies_counts[sample_name] = mis_val

        # Residual hits
        contamination_count = int(row["contamination_count"])
        total_residual_hits += contamination_count
        if contamination_count == 0:
            clean_isolates_count += 1
        else:
            affected_residual_count += 1

        # Replicons
        full_missed = int(row["full_missed"])
        partial_missed = int(row["partial_missed"])
        total_missed = int(row["total_missed"])
        full_missed_sum += full_missed
        partial_missed_sum += partial_missed
        total_missed_sum += total_missed
        if total_missed > 0:
            affected_replicon_count += 1

        # Parse replicons and verify completeness gate (>= 95% for all replicons)
        _, is_complete = validate_and_parse_replicon_row(dict(row))
        if not is_complete:
            all_complete_recovery = False

    # Cohort scores
    score_contiguity = float(np.mean(list(contiguity_scores.values())))
    score_accuracy = float(np.mean(list(accuracy_scores.values())))
    score_residual = float((clean_isolates_count / n_isolates) * 100.0)
    score_replicon = score_replicon_recovery(total_missed_sum, denominator=replicon_denominator)

    return CohortCriteriaResult(
        score_contiguity=score_contiguity,
        score_accuracy=score_accuracy,
        score_residual=score_residual,
        score_replicon=score_replicon,
        mean_auNGA_ratio=aunga_sum / n_isolates,
        mean_error_rate=(mismatches_sum + indels_sum) / n_isolates,
        mean_mismatches=mismatches_sum / n_isolates,
        mean_indels=indels_sum / n_isolates,
        residual_total_hits=total_residual_hits,
        residual_clean_isolates=clean_isolates_count,
        residual_affected_isolates=affected_residual_count,
        replicon_total_missed=total_missed_sum,
        replicon_full_missed=full_missed_sum,
        replicon_partial_missed=partial_missed_sum,
        replicon_affected_isolates=affected_replicon_count,
        all_replicons_complete=all_complete_recovery,
        zero_residual_hits=(total_residual_hits == 0),
        mean_duplication_ratio=duplication_sum / n_isolates,
        total_misassemblies=misassemblies_sum,
        isolate_contiguity_scores=contiguity_scores,
        isolate_accuracy_scores=accuracy_scores,
        isolate_duplication_ratios=duplication_ratios,
        isolate_misassemblies=misassemblies_counts,
    )
