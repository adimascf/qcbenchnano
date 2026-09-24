import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.validation import (
    UnsupportedScenarioError,
    InputValidationError,
    validate_scenario,
    validate_benchmark_dataframe,
    filter_scenario_data,
    EXPECTED_COMBINATIONS,
    EXPECTED_SAMPLES,
)

def test_supported_scenarios():
    for model in ["hac", "sup"]:
        for depth in ["20x", "100x"]:
            sc = Scenario(model=model, depth=depth)
            assert sc.model == model
            assert sc.depth == depth
            validate_scenario(sc)

def test_unsupported_scenarios():
    with pytest.raises(UnsupportedScenarioError):
        validate_scenario(Scenario(model="all", depth="20x"))
    with pytest.raises(UnsupportedScenarioError):
        validate_scenario(Scenario(model="hac", depth="50x"))
    with pytest.raises(UnsupportedScenarioError):
        validate_scenario(Scenario(model="hac", depth="all"))
    with pytest.raises(UnsupportedScenarioError):
        validate_scenario(Scenario(model="fast", depth="100x"))

def test_dataframe_missing_required_columns():
    df = pd.DataFrame({
        "combo": ["c1"],
        "depth": ["20x"],
        "sample": ["s1"],
        "model": ["hac"],
    })
    with pytest.raises(InputValidationError, match="Missing required columns"):
        validate_benchmark_dataframe(df)

def test_dataframe_duplicate_keys():
    # Duplicate combo, sample, model, depth
    df = pd.DataFrame({
        "combo": ["c1", "c1"],
        "depth": ["20x", "20x"],
        "sample": ["s1", "s1"],
        "model": ["hac", "hac"],
        "Mismatches per 100kbp": [0.0, 0.0],
        "Indels per 100kbp": [0.0, 0.0],
        "auNGA_ratio": [1.0, 1.0],
        "contamination_count": [0, 0],
        "full_missed": [0, 0],
        "partial_missed": [0, 0],
        "total_missed": [0, 0],
        "all_contigs_coverage": ["chr (1000bp, 100.0% cov)", "chr (1000bp, 100.0% cov)"],
    })
    with pytest.raises(InputValidationError, match="Duplicate observation keys"):
        validate_benchmark_dataframe(df)

def test_dataframe_invalid_numerics():
    df = pd.DataFrame({
        "combo": ["c1"],
        "depth": ["20x"],
        "sample": ["s1"],
        "model": ["hac"],
        "Mismatches per 100kbp": [float("nan")],
        "Indels per 100kbp": [0.0],
        "auNGA_ratio": [1.0],
        "contamination_count": [0],
        "full_missed": [0],
        "partial_missed": [0],
        "total_missed": [0],
        "all_contigs_coverage": ["chr (1000bp, 100.0% cov)"],
    })
    with pytest.raises(InputValidationError, match="NaN or infinite values"):
        validate_benchmark_dataframe(df)

def test_filter_scenario_data():
    df = pd.read_csv("logbook/assembly_metrics.csv")
    sc = Scenario(model="hac", depth="100x")
    filtered, incomplete = filter_scenario_data(df, sc)
    assert len(filtered) == 17 * 13
    assert incomplete == {}
    assert set(filtered["combo"].unique()) == EXPECTED_COMBINATIONS
    assert set(filtered["sample"].unique()) == EXPECTED_SAMPLES
