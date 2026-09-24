import pytest
import numpy as np
import pandas as pd
from qc_scoring.criteria import (
    score_isolate_contiguity,
    score_isolate_accuracy,
    score_isolate_residual_clean,
    calculate_cohort_criteria,
    CohortCriteriaResult,
)

def test_contiguity_value_function_anchors():
    assert score_isolate_contiguity(1.0) == 100.0
    assert score_isolate_contiguity(0.9) == pytest.approx(90.0)
    assert score_isolate_contiguity(1.1) == pytest.approx(90.0)
    assert score_isolate_contiguity(0.0) == 0.0
    assert score_isolate_contiguity(-0.5) == 0.0
    assert score_isolate_contiguity(2.0) == 0.0
    assert score_isolate_contiguity(2.5) == 0.0

def test_accuracy_value_function_anchors():
    assert score_isolate_accuracy(0.0, 0.0) == 100.0
    assert score_isolate_accuracy(2.5, 2.5) == pytest.approx(50.0)
    assert score_isolate_accuracy(5.0, 5.0) == 0.0
    assert score_isolate_accuracy(6.0, 6.0) == 0.0  # clipped at 0
    assert score_isolate_accuracy(-1.0, 0.0) == 100.0  # clipped at 100 if negative rate

def test_residual_removal_isolate_binary():
    assert score_isolate_residual_clean(0) == 100.0
    assert score_isolate_residual_clean(1) == 0.0
    assert score_isolate_residual_clean(5) == 0.0

def test_replicon_value_function_anchors():
    from qc_scoring.criteria import score_replicon_recovery
    assert score_replicon_recovery(0) == 100.0
    assert score_replicon_recovery(1) == pytest.approx(100.0 * (2.0 / 3.0))
    assert score_replicon_recovery(2) == pytest.approx(100.0 * (1.0 / 3.0))
    assert score_replicon_recovery(3) == 0.0
    assert score_replicon_recovery(4) == 0.0

def test_order_of_operations_distinction():
    # If isolate 1 has auNGA=0.0 (score 0) and isolate 2 has auNGA=2.0 (score 0),
    # isolate scoring -> mean is (0 + 0)/2 = 0.
    # If one averaged raw auNGA first: (0 + 2)/2 = 1.0 -> score would be 100!
    df = pd.DataFrame({
        "sample": ["s1", "s2"],
        "auNGA_ratio": [0.0, 2.0],
        "Mismatches per 100kbp": [0.0, 10.0],
        "Indels per 100kbp": [0.0, 10.0],
        "contamination_count": [0, 1],
        "full_missed": [0, 0],
        "partial_missed": [0, 0],
        "total_missed": [0, 0],
        "all_contigs_coverage": ["chr (1000bp, 100% cov)", "chr (1000bp, 100% cov)"],
    })
    res = calculate_cohort_criteria(df)
    assert res.score_contiguity == 0.0  # Not 100!
    # Accuracy: s1 is (0+0)->100, s2 is (10+10=20)->0. Cohort mean should be 50.0.
    # If raw averaged: (10+10)/2 = 10 -> score 0!
    assert res.score_accuracy == 50.0  # Not 0!
