from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from qc_scoring.models import Scenario

SUPPORTED_MODELS: Set[str] = {"hac", "sup"}
SUPPORTED_DEPTHS: Set[str] = {"20x", "100x"}

EXPECTED_COMBINATIONS: Set[str] = {
    "chopper-barbell",
    "chopper-dorado",
    "chopper-porechop_abi",
    "chopper-untrimmed",
    "fastplong-all",
    "filtlong-barbell",
    "filtlong-dorado",
    "filtlong-porechop_abi",
    "filtlong-untrimmed",
    "seqkit-barbell",
    "seqkit-dorado",
    "seqkit-porechop_abi",
    "seqkit-untrimmed",
    "unprocessed-barbell",
    "unprocessed-dorado",
    "unprocessed-porechop_abi",
    "unprocessed-untrimmed",
}

EXPECTED_SAMPLES: Set[str] = {
    "AJ292__202310",
    "AMtb_1__202402",
    "ATCC_10708__202309",
    "ATCC_17802__202309",
    "ATCC_19119__202309",
    "ATCC_25922__202309",
    "ATCC_33560__202309",
    "ATCC_35221__202309",
    "ATCC_35897__202309",
    "ATCC_BAA-679__202309",
    "BPH2947__202310",
    "MMC234__202311",
    "RDH275__202311",
}

REQUIRED_COLUMNS: List[str] = [
    "combo",
    "depth",
    "sample",
    "model",
    "Mismatches per 100kbp",
    "Indels per 100kbp",
    "auNGA_ratio",
    "contamination_count",
    "full_missed",
    "partial_missed",
    "total_missed",
    "all_contigs_coverage",
]

NUMERIC_COLUMNS: List[str] = [
    "Mismatches per 100kbp",
    "Indels per 100kbp",
    "auNGA_ratio",
    "contamination_count",
    "full_missed",
    "partial_missed",
    "total_missed",
]


class ScoringError(Exception):
    """Base exception for scoring errors."""
    pass


class UnsupportedScenarioError(ScoringError):
    """Raised when an unsupported basecalling model or sequencing depth is requested."""
    pass


class InputValidationError(ScoringError):
    """Raised when input benchmark data is malformed or invalid."""
    pass


def validate_scenario(scenario: Scenario) -> None:
    if scenario.model not in SUPPORTED_MODELS:
        raise UnsupportedScenarioError(
            f"Unsupported basecalling model '{scenario.model}'. Supported models: {sorted(SUPPORTED_MODELS)}"
        )
    if scenario.depth not in SUPPORTED_DEPTHS:
        raise UnsupportedScenarioError(
            f"Unsupported sequencing depth '{scenario.depth}'. Supported depths: {sorted(SUPPORTED_DEPTHS)}"
        )


def validate_benchmark_dataframe(df: pd.DataFrame) -> None:
    missing_cols = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_cols:
        raise InputValidationError(f"Missing required columns: {missing_cols}")

    # Check duplicate observation keys
    key_cols = ["model", "depth", "combo", "sample"]
    dups = df[df.duplicated(subset=key_cols, keep=False)]
    if not dups.empty:
        raise InputValidationError(
            f"Duplicate observation keys found for {len(dups)} rows on {key_cols}"
        )

    # Check numeric columns for NaN / inf
    for col in NUMERIC_COLUMNS:
        vals = df[col]
        if vals.isna().any() or np.isinf(vals).any():
            raise InputValidationError(f"NaN or infinite values found in numeric column '{col}'")
        if (vals < 0).any():
            raise InputValidationError(f"Negative values found in column '{col}' where non-negative required")


def filter_scenario_data(
    df: pd.DataFrame, scenario: Scenario, expected_samples: Optional[Set[str]] = None
) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    Validates scenario and filters benchmark DataFrame to matching rows.
    Identifies combinations with missing isolates as incomplete (insufficient data).
    Returns (filtered_complete_df, incomplete_reasons_dict).
    """
    validate_scenario(scenario)
    validate_benchmark_dataframe(df)

    if expected_samples is None:
        target_samples = EXPECTED_SAMPLES
    else:
        target_samples = expected_samples

    mask = (df["model"] == scenario.model) & (df["depth"] == scenario.depth)
    scenario_df = df[mask].copy()

    incomplete: Dict[str, str] = {}
    valid_combos: List[str] = []

    # Check each combo for completeness across expected isolates
    combos_present = set(scenario_df["combo"].unique())
    all_combos = EXPECTED_COMBINATIONS.union(combos_present)

    for combo in all_combos:
        combo_rows = scenario_df[scenario_df["combo"] == combo]
        samples_present = set(combo_rows["sample"].unique())
        missing_samples = target_samples - samples_present
        if missing_samples:
            incomplete[combo] = (
                f"insufficient benchmark data: missing {len(missing_samples)} isolate(s): "
                f"{sorted(missing_samples)}"
            )
        else:
            valid_combos.append(combo)

    filtered_df = scenario_df[scenario_df["combo"].isin(valid_combos)].copy()
    return filtered_df, incomplete
