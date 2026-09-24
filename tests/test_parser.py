import pytest
from qc_scoring.parser import (
    parse_replicon_coverage,
    validate_and_parse_replicon_row,
    RepliconParsingError,
    RepliconDisagreementError,
    RepliconInfo,
)

def test_parse_single_replicon():
    text = "chromosome (5445132bp, 99.9995% cov)"
    reps = parse_replicon_coverage(text)
    assert len(reps) == 1
    assert reps[0] == RepliconInfo(name="chromosome", size_bp=5445132, coverage_pct=99.9995)

def test_parse_multiple_replicons():
    text = "chromosome (4752118bp, 99.9996% cov); plasmid (49586bp, 99.9919% cov)"
    reps = parse_replicon_coverage(text)
    assert len(reps) == 2
    assert reps[0].name == "chromosome"
    assert reps[0].size_bp == 4752118
    assert reps[0].coverage_pct == 99.9996
    assert reps[1].name == "plasmid"
    assert reps[1].size_bp == 49586
    assert reps[1].coverage_pct == 99.9919

def test_replicon_threshold_boundaries():
    # 0% coverage -> full missed
    r_zero = parse_replicon_coverage("plasmid (1000bp, 0.0% cov)")[0]
    assert r_zero.is_full_missed is True
    assert r_zero.is_partial_missed is False
    assert r_zero.is_below_50 is True
    assert r_zero.is_complete_recovery is False

    # 49.99% coverage -> just below 50% -> partial missed
    r_49 = parse_replicon_coverage("plasmid (1000bp, 49.99% cov)")[0]
    assert r_49.is_full_missed is False
    assert r_49.is_partial_missed is True
    assert r_49.is_below_50 is True
    assert r_49.is_complete_recovery is False

    # 50.0% coverage -> exactly 50% -> NOT missed (< 50% required), but < 95%
    r_50 = parse_replicon_coverage("plasmid (1000bp, 50.0% cov)")[0]
    assert r_50.is_full_missed is False
    assert r_50.is_partial_missed is False
    assert r_50.is_below_50 is False
    assert r_50.is_complete_recovery is False

    # 94.99% coverage -> just below 95% -> not complete recovery
    r_94 = parse_replicon_coverage("plasmid (1000bp, 94.99% cov)")[0]
    assert r_94.is_below_50 is False
    assert r_94.is_complete_recovery is False

    # 95.0% coverage -> exactly 95% -> complete recovery
    r_95 = parse_replicon_coverage("plasmid (1000bp, 95.0% cov)")[0]
    assert r_95.is_below_50 is False
    assert r_95.is_complete_recovery is True

    # 100.0% coverage -> complete recovery
    r_100 = parse_replicon_coverage("plasmid (1000bp, 100.0% cov)")[0]
    assert r_100.is_complete_recovery is True

def test_malformed_entries():
    with pytest.raises(RepliconParsingError):
        parse_replicon_coverage("malformed string without format")
    with pytest.raises(RepliconParsingError):
        parse_replicon_coverage("chr (1000, 99% cov)")  # missing bp
    with pytest.raises(RepliconParsingError):
        parse_replicon_coverage("chr (1000bp, 99%)")  # missing cov
    with pytest.raises(RepliconParsingError):
        parse_replicon_coverage("")  # empty text
    with pytest.raises(RepliconParsingError):
        parse_replicon_coverage("chr (1000bp, 105.0% cov)")  # > 100%
    with pytest.raises(RepliconParsingError):
        parse_replicon_coverage("chr (1000bp, -5.0% cov)")  # < 0%

def test_row_cross_check_success():
    row = {
        "all_contigs_coverage": "chr (5000bp, 100.0% cov); p1 (2000bp, 0.0% cov); p2 (1000bp, 30.0% cov)",
        "full_missed": 1,
        "partial_missed": 1,
        "total_missed": 2,
    }
    reps, is_complete = validate_and_parse_replicon_row(row)
    assert len(reps) == 3
    assert is_complete is False

def test_row_cross_check_disagreement():
    # Disagreement on full_missed
    row = {
        "all_contigs_coverage": "chr (5000bp, 100.0% cov); p1 (2000bp, 0.0% cov)",
        "full_missed": 0,  # contradicts p1 with 0.0% cov
        "partial_missed": 0,
        "total_missed": 0,
    }
    with pytest.raises(RepliconDisagreementError):
        validate_and_parse_replicon_row(row)

    # Disagreement on total_missed
    row2 = {
        "all_contigs_coverage": "chr (5000bp, 100.0% cov); p1 (2000bp, 30.0% cov)",
        "full_missed": 0,
        "partial_missed": 1,
        "total_missed": 2,  # should be 1
    }
    with pytest.raises(RepliconDisagreementError):
        validate_and_parse_replicon_row(row2)
