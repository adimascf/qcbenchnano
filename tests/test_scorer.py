import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import (
    WeightsConfig,
    GateConfig,
    PresetName,
    get_preset,
    validate_weights,
    InvalidWeightsError,
)
from qc_scoring.scorer import (
    score_benchmark,
    ScoringResult,
    CombinationRecommendation,
)
from qc_scoring.validation import EXPECTED_SAMPLES

def test_presets_weights_and_gates():
    # Community-balanced
    cb = get_preset(PresetName.COMMUNITY_BALANCED)
    assert cb.weights == WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0)
    assert cb.gates == GateConfig(complete_recovery=False, zero_residual_hits=False)

    # Complete-replicon-recovery
    cr = get_preset(PresetName.COMPLETE_REPLICON_RECOVERY)
    assert cr.weights == WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0)
    assert cr.gates == GateConfig(complete_recovery=True, zero_residual_hits=False)

    # Sequence-accurate-assembly
    sa = get_preset(PresetName.SEQUENCE_ACCURATE_ASSEMBLY)
    assert sa.weights == WeightsConfig(accuracy=50.0, contiguity=14.0, residual=12.0, replicon=24.0)
    assert sa.gates == GateConfig(complete_recovery=False, zero_residual_hits=False)

def test_weight_validation_rules():
    # Valid
    validate_weights(WeightsConfig(accuracy=25.0, contiguity=25.0, residual=25.0, replicon=25.0))
    validate_weights(WeightsConfig(accuracy=100.0, contiguity=0.0, residual=0.0, replicon=0.0))
    
    # Negative weight
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=-5.0, contiguity=35.0, residual=35.0, replicon=35.0))
        
    # All zero
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=0.0, contiguity=0.0, residual=0.0, replicon=0.0))
        
    # Total below 100
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=20.0, contiguity=20.0, residual=20.0, replicon=20.0))
        
    # Total above 100
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=30.0, contiguity=30.0, residual=30.0, replicon=30.0))

def test_all_ineligible_scenario():
    # Load benchmark data and force strict gates that cause every combination to fail
    df = pd.read_csv("logbook/assembly_metrics.csv")
    scenario = Scenario(model="hac", depth="20x")
    # Complete recovery gate + zero residual hits gate
    result = score_benchmark(
        df,
        scenario,
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(complete_recovery=True, zero_residual_hits=True)
    )
    # Check that even if combinations fail, gates were not relaxed
    assert result.scenario == scenario
    for rec in result.recommendations:
        if not rec.is_eligible:
            assert rec.rank is None
            assert rec.ineligible_reason is not None

def test_near_ties_and_exact_ties():
    # Construct scenario with exact tie and near tie using expected isolate names
    sample_list = sorted(EXPECTED_SAMPLES)
    records = []
    for c in ["chopper-barbell", "chopper-dorado", "chopper-untrimmed"]:
        for sample in sample_list:
            if c in ("chopper-barbell", "chopper-dorado"):
                aunga = 1.0
                err = 0.0
            else:  # chopper-untrimmed slightly worse
                aunga = 1.0
                err = 0.05  # rate 0.05 -> score 99.5 instead of 100.0 (diff 0.5 < 1.0)
            records.append({
                "combo": c,
                "depth": "100x",
                "sample": sample,
                "model": "hac",
                "Mismatches per 100kbp": err,
                "Indels per 100kbp": 0.0,
                "auNGA_ratio": aunga,
                "contamination_count": 0,
                "full_missed": 0,
                "partial_missed": 0,
                "total_missed": 0,
                "all_contigs_coverage": "chr (1000000bp, 100.0% cov)",
            })
    df_tie = pd.DataFrame(records)
    res = score_benchmark(
        df_tie,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=50.0, contiguity=50.0, residual=0.0, replicon=0.0),
        gates=GateConfig()
    )
    # Filter to the 3 evaluated combos in eligible list
    eval_recs = [r for r in res.recommendations if r.combo in ("chopper-barbell", "chopper-dorado", "chopper-untrimmed")]
    # chopper-barbell and chopper-dorado should have exact same rank 1
    assert eval_recs[0].combo == "chopper-barbell"
    assert eval_recs[0].rank == 1
    assert eval_recs[1].combo == "chopper-dorado"
    assert eval_recs[1].rank == 1
    
    # chopper-untrimmed should have rank 3
    assert eval_recs[2].combo == "chopper-untrimmed"
    assert eval_recs[2].rank == 3
    
    # Near tie flags: adjacent difference is 0.25 (< 1.0)
    assert eval_recs[1].is_near_tie is True
    assert eval_recs[2].is_near_tie is True

def test_warning_25_point_gap():
    # Test accuracy variability warning: mean - min >= 25.0
    sample_list = sorted(EXPECTED_SAMPLES)
    records = []
    for idx, sample in enumerate(sample_list):
        # 12 isolates score 100 (err=0), 1 isolate scores 70 (err=3.0) -> mean = 97.69
        # min = 70.0. Gap = 97.69 - 70.0 = 27.69 >= 25.0 -> should warn!
        err = 3.0 if idx == 0 else 0.0
        records.append({
            "combo": "chopper-barbell",
            "depth": "100x",
            "sample": sample,
            "model": "hac",
            "Mismatches per 100kbp": err,
            "Indels per 100kbp": 0.0,
            "auNGA_ratio": 1.0,
            "contamination_count": 0,
            "full_missed": 0,
            "partial_missed": 0,
            "total_missed": 0,
            "all_contigs_coverage": "chr (1000000bp, 100.0% cov)",
        })
    df_warn = pd.DataFrame(records)
    res = score_benchmark(
        df_warn,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=100.0, contiguity=0.0, residual=0.0, replicon=0.0),
        gates=GateConfig()
    )
    rec = [r for r in res.recommendations if r.combo == "chopper-barbell"][0]
    has_var_warn = any("Variable sequence accuracy" in w for w in rec.warnings)
    assert has_var_warn is True

def test_insufficient_data_missing_isolate():
    # Omit 1 isolate for chopper-barbell
    df = pd.read_csv("logbook/assembly_metrics.csv")
    mask = ~((df["combo"] == "chopper-barbell") & (df["sample"] == "AJ292__202310"))
    df_missing = df[mask].copy()

    res = score_benchmark(
        df_missing,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(),
    )
    rec = [r for r in res.recommendations if r.combo == "chopper-barbell"][0]
    assert rec.is_eligible is False
    assert rec.rank is None
    assert "insufficient benchmark data" in (rec.ineligible_reason or "")
    assert "AJ292__202310" in (rec.ineligible_reason or "")

def test_pinned_real_data_acceptance_fixture():
    df = pd.read_csv("logbook/assembly_metrics.csv")
    cb_preset = get_preset(PresetName.COMMUNITY_BALANCED)

    for model in ["hac", "sup"]:
        for depth in ["20x", "100x"]:
            sc = Scenario(model=model, depth=depth)
            result = score_benchmark(
                df,
                sc,
                weights=cb_preset.weights,
                gates=cb_preset.gates,
            )
            assert len(result.recommendations) == 17
            assert result.scoring_version == "1.0"
            assert result.scenario == sc
            assert result.provenance is not None
            # Every recommendation must have 4 criterion scores in [0, 100]
            for rec in result.recommendations:
                assert 0.0 <= rec.score_contiguity <= 100.0
                assert 0.0 <= rec.score_accuracy <= 100.0
                assert 0.0 <= rec.score_residual <= 100.0
                assert 0.0 <= rec.score_replicon <= 100.0
                if rec.is_eligible:
                    assert rec.rank is not None
                    assert 0.0 <= rec.overall_score <= 100.0
                    assert rec.display_score == round(rec.overall_score, 1)

def test_export_recommendations_csv():
    from qc_scoring.scorer import export_recommendations_csv, recommendations_to_dataframe
    df = pd.read_csv("logbook/assembly_metrics.csv")
    res = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(),
    )
    csv_text = export_recommendations_csv(res)
    assert len(csv_text) > 0
    lines = csv_text.strip().split("\n")
    assert len(lines) == 18  # 1 header + 17 rows
    export_df = recommendations_to_dataframe(res)
    assert len(export_df) == 17
    assert "rank" in export_df.columns
    assert "combo" in export_df.columns
    assert "overall_score" in export_df.columns
    assert "scoring_version" in export_df.columns
    assert "source_data_commit" in export_df.columns


def test_complete_replicon_recovery_preset_and_gate():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    preset = get_preset(PresetName.COMPLETE_REPLICON_RECOVERY)
    assert preset.weights == WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0)
    assert preset.gates == GateConfig(complete_recovery=True, zero_residual_hits=False)

    res = score_benchmark(df, Scenario(model="hac", depth="100x"), weights=preset.weights, gates=preset.gates)
    assert res.eligible_count == 7
    assert res.ineligible_count == 10
    leader = res.leading_recommendation
    assert leader is not None
    assert leader.combo == "seqkit-barbell"
    assert leader.rank == 1
    assert leader.display_score == 98.0

    # Ineligible combinations must receive rank=None and exact failure reasons
    for rec in res.recommendations:
        if not rec.is_eligible:
            assert rec.rank is None
            assert "fails complete-recovery gate (<95% coverage on at least one replicon)" in (rec.ineligible_reason or "")


def test_sequence_accurate_assembly_preset_and_ranking():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    preset = get_preset(PresetName.SEQUENCE_ACCURATE_ASSEMBLY)
    assert preset.weights == WeightsConfig(accuracy=50.0, contiguity=14.0, residual=12.0, replicon=24.0)
    assert preset.gates == GateConfig(complete_recovery=False, zero_residual_hits=False)

    res = score_benchmark(df, Scenario(model="hac", depth="100x"), weights=preset.weights, gates=preset.gates)
    assert res.eligible_count == 17
    assert res.ineligible_count == 0
    leader = res.leading_recommendation
    assert leader is not None
    assert leader.combo == "seqkit-barbell"
    assert leader.rank == 1
    assert leader.display_score == 97.6


def test_custom_decimal_percentages_and_edge_cases():
    # Valid decimal percentages totaling exactly 100.0
    w_dec = WeightsConfig(accuracy=33.3, contiguity=33.3, residual=16.7, replicon=16.7)
    validate_weights(w_dec)

    # Valid with single non-zero weight
    w_single = WeightsConfig(accuracy=100.0, contiguity=0.0, residual=0.0, replicon=0.0)
    validate_weights(w_single)

    # Valid with zeros in other fields
    w_zero = WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0)
    validate_weights(w_zero)

    # Invalid: total is 99.9 (fails exact 100)
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=33.3, contiguity=33.3, residual=16.6, replicon=16.7))

    # Invalid: total is 100.1
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=33.3, contiguity=33.4, residual=16.7, replicon=16.7))

    # Invalid: negative decimal
    with pytest.raises(InvalidWeightsError):
        validate_weights(WeightsConfig(accuracy=-0.5, contiguity=50.5, residual=25.0, replicon=25.0))


def test_multiple_exclusion_reasons_and_reasons_format():
    # Load benchmark data and apply both gates
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    res = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(complete_recovery=True, zero_residual_hits=True),
    )
    # Find unprocessed-dorado which has both residual hits and sub-95% replicons
    rec = [r for r in res.recommendations if r.combo == "unprocessed-dorado"][0]
    assert rec.is_eligible is False
    assert rec.rank is None
    reason = rec.ineligible_reason or ""
    assert "fails complete-recovery gate (<95% coverage on at least one replicon)" in reason
    assert "fails zero-residual-hits gate" in reason


def test_complete_recovery_sub_95_gate_edge_case():
    # Construct a case where total_missed == 0 (no replicon < 50%), but one replicon has 94.0% coverage (< 95%)
    sample_list = sorted(EXPECTED_SAMPLES)
    records = []
    for idx, sample in enumerate(sample_list):
        # 12 isolates have 100% cov, 1 isolate has 94.0% cov on plasmid
        cov_str = "chr (1000000bp, 100.0% cov); plas (50000bp, 94.0% cov)" if idx == 0 else "chr (1000000bp, 100.0% cov)"
        records.append({
            "combo": "seqkit-dorado",
            "depth": "100x",
            "sample": sample,
            "model": "hac",
            "Mismatches per 100kbp": 1.0,
            "Indels per 100kbp": 0.0,
            "auNGA_ratio": 1.0,
            "contamination_count": 0,
            "full_missed": 0,
            "partial_missed": 0,
            "total_missed": 0,
            "all_contigs_coverage": cov_str,
        })
    df_edge = pd.DataFrame(records)
    res = score_benchmark(
        df_edge,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0),
        gates=GateConfig(complete_recovery=True, zero_residual_hits=False),
    )
    rec = [r for r in res.recommendations if r.combo == "seqkit-dorado"][0]
    # Even though total_missed == 0, 94.0% is < 95.0%, so it must FAIL the complete recovery gate!
    assert rec.is_eligible is False
    assert rec.rank is None
    assert "fails complete-recovery gate (<95% coverage on at least one replicon)" in (rec.ineligible_reason or "")


def test_warning_boundaries_immediately_below_and_at_25_points():
    sample_list = sorted(EXPECTED_SAMPLES)

    # 1. Test sequence accuracy at exactly 25.0 points difference:
    # All 13 isolates: 12 have score 100.0 (err=0), 1 has score 72.91666666666667
    # Mean = (12 * 100 + S) / 13 = (1200 + S) / 13.
    # We want Mean - S = 25.0 => 1200 + S - 13S = 25 * 13 => 1200 - 12S = 325 => 12S = 875 => S = 72.91666666666667
    # Then Mean = 97.91666666666667, S = 72.91666666666667, Mean - S = 25.0!
    # With score_isolate_accuracy: S = 100 * (1 - err/10) => err = 10 * (1 - S/100) = 10 * (1 - 0.7291666666666667) = 2.708333333333333
    target_s_exact = 875.0 / 12.0
    err_exact = 10.0 * (1.0 - target_s_exact / 100.0)

    records_exact = []
    for idx, sample in enumerate(sample_list):
        records_exact.append({
            "combo": "filtlong-dorado",
            "depth": "100x",
            "sample": sample,
            "model": "hac",
            "Mismatches per 100kbp": err_exact if idx == 0 else 0.0,
            "Indels per 100kbp": 0.0,
            "auNGA_ratio": 1.0,
            "contamination_count": 0,
            "full_missed": 0,
            "partial_missed": 0,
            "total_missed": 0,
            "all_contigs_coverage": "chr (1000000bp, 100.0% cov)",
            "Duplication_ratio": 1.0,
            "misassemblies": 0,
        })

    res_exact = score_benchmark(
        pd.DataFrame(records_exact),
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
    )
    rec_exact = [r for r in res_exact.recommendations if r.combo == "filtlong-dorado"][0]
    acc_warns_exact = [w for w in rec_exact.warnings if "Variable sequence accuracy" in w]
    assert len(acc_warns_exact) == 1
    # Check that it reports both the lowest isolate score and cohort mean
    warn_text = acc_warns_exact[0]
    assert f"{target_s_exact:.1f}" in warn_text
    assert f"{rec_exact.score_accuracy:.1f}" in warn_text

    # 2. Test immediately below 25.0 points difference (e.g. difference is 24.9 points):
    # Mean - S = 24.9 => 1200 - 12S = 24.9 * 13 = 323.7 => 12S = 876.3 => S = 73.025
    target_s_below = 876.3 / 12.0
    err_below = 10.0 * (1.0 - target_s_below / 100.0)

    records_below = []
    for idx, sample in enumerate(sample_list):
        records_below.append({
            "combo": "filtlong-dorado",
            "depth": "100x",
            "sample": sample,
            "model": "hac",
            "Mismatches per 100kbp": err_below if idx == 0 else 0.0,
            "Indels per 100kbp": 0.0,
            "auNGA_ratio": 1.0,
            "contamination_count": 0,
            "full_missed": 0,
            "partial_missed": 0,
            "total_missed": 0,
            "all_contigs_coverage": "chr (1000000bp, 100.0% cov)",
            "Duplication_ratio": 1.0,
            "misassemblies": 0,
        })

    res_below = score_benchmark(
        pd.DataFrame(records_below),
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
    )
    rec_below = [r for r in res_below.recommendations if r.combo == "filtlong-dorado"][0]
    acc_warns_below = [w for w in rec_below.warnings if "Variable sequence accuracy" in w]
    assert len(acc_warns_below) == 0, f"Expected no warning for 24.9 pt gap, got {acc_warns_below}"

    # 3. Test contiguity variability boundary: exactly 25.0 points vs 24.9 points
    # Contiguity score = 100 * max(0, 1 - |r - 1|).
    # S = 100 * (1 - (r - 1)) = 100 * (2 - r) => r = 2 - S/100
    r_exact = 2.0 - (target_s_exact / 100.0)
    for r in records_exact:
        r["Mismatches per 100kbp"] = 0.0
    records_exact[0]["auNGA_ratio"] = r_exact

    res_cont_exact = score_benchmark(
        pd.DataFrame(records_exact),
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
    )
    rec_cont_exact = [r for r in res_cont_exact.recommendations if r.combo == "filtlong-dorado"][0]
    cont_warns_exact = [w for w in rec_cont_exact.warnings if "Variable contiguity" in w]
    assert len(cont_warns_exact) == 1
    assert f"{target_s_exact:.1f}" in cont_warns_exact[0]
    assert f"{rec_cont_exact.score_contiguity:.1f}" in cont_warns_exact[0]

    # Contiguity immediately below 25 points:
    r_below = 2.0 - (target_s_below / 100.0)
    for r in records_below:
        r["Mismatches per 100kbp"] = 0.0
    records_below[0]["auNGA_ratio"] = r_below

    res_cont_below = score_benchmark(
        pd.DataFrame(records_below),
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
    )
    rec_cont_below = [r for r in res_cont_below.recommendations if r.combo == "filtlong-dorado"][0]
    cont_warns_below = [w for w in rec_cont_below.warnings if "Variable contiguity" in w]
    assert len(cont_warns_below) == 0


def test_warnings_remain_visible_when_criterion_weight_is_zero():
    # In real data, test with zero weights for various criteria
    df = pd.read_csv("logbook/assembly_metrics.csv")

    # Zero residual weight: residual warnings still present
    res_zero_res = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=50.0, contiguity=50.0, residual=0.0, replicon=0.0),
    )
    # Check combos that have residual hits
    recs_with_hits = [r for r in res_zero_res.recommendations if r.residual_total_hits > 0]
    assert len(recs_with_hits) > 0
    for r in recs_with_hits:
        hit_warns = [w for w in r.warnings if "residual adapter/barcode hit(s)" in w]
        assert len(hit_warns) == 1, f"Expected residual hit warning despite 0 weight for {r.combo}"

    # Zero replicon weight: replicon warnings still present
    res_zero_rep = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=50.0, contiguity=50.0, residual=0.0, replicon=0.0),
    )
    recs_with_miss = [r for r in res_zero_rep.recommendations if r.replicon_total_missed > 0]
    assert len(recs_with_miss) > 0
    for r in recs_with_miss:
        miss_warns = [w for w in r.warnings if "missed or severely incomplete replicon(s)" in w]
        assert len(miss_warns) == 1, f"Expected replicon warning despite 0 weight for {r.combo}"

    # Zero accuracy weight: variable accuracy warnings still present
    res_zero_acc = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=0.0, contiguity=50.0, residual=25.0, replicon=25.0),
    )
    recs_var_acc = [
        r for r in res_zero_acc.recommendations
        if r.isolate_accuracy_scores and (r.score_accuracy - min(r.isolate_accuracy_scores.values())) >= 25.0
    ]
    assert len(recs_var_acc) > 0
    for r in recs_var_acc:
        acc_warns = [w for w in r.warnings if "Variable sequence accuracy" in w]
        assert len(acc_warns) == 1, f"Expected var acc warning despite 0 weight for {r.combo}"


def test_supporting_metrics_duplication_ratio_and_misassemblies():
    df = pd.read_csv("logbook/assembly_metrics.csv")
    res = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
    )
    for r in res.recommendations:
        if r.is_eligible:
            assert hasattr(r, "mean_duplication_ratio")
            assert hasattr(r, "total_misassemblies")
            assert r.mean_duplication_ratio > 0.0
            assert r.total_misassemblies >= 0


def test_export_recommendations_csv_contract():
    import io
    from qc_scoring.scorer import export_recommendations_csv, recommendations_to_dataframe
    from qc_scoring.validation import UnsupportedScenarioError

    df = pd.read_csv("assets/data/assembly_metrics.csv")
    scenario = Scenario(model="hac", depth="100x")
    weights = WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0)

    # 1. Community-balanced export test
    result = score_benchmark(df, scenario, weights=weights, preset="community_balanced")
    assert result.preset == "community_balanced"
    assert result.provenance.preset == "community_balanced"

    csv_text = export_recommendations_csv(result)
    export_df = pd.read_csv(io.StringIO(csv_text), keep_default_na=False)

    assert len(export_df) == 17, "Export must contain all 17 combinations"
    
    # Check all required columns
    required_cols = [
        "rank", "combo", "overall_score", "display_score", "is_eligible", "is_insufficient_data",
        "ineligible_reason", "is_near_tie", "near_tie_label",
        "score_contiguity", "score_accuracy", "score_residual", "score_replicon",
        "mean_auNGA_ratio", "mean_error_rate", "mean_mismatches", "mean_indels",
        "residual_total_hits", "residual_clean_isolates", "residual_affected_isolates",
        "replicon_total_missed", "replicon_full_missed", "replicon_partial_missed", "replicon_affected_isolates",
        "mean_duplication_ratio", "total_misassemblies",
        "warnings", "preset", "scenario_model", "scenario_depth",
        "weight_accuracy", "weight_contiguity", "weight_residual", "weight_replicon",
        "gate_complete_recovery", "gate_zero_residual_hits",
        "scoring_version", "source_data_commit", "source_data_hash", "generated_at"
    ]
    for col in required_cols:
        assert col in export_df.columns, f"Missing required column: {col}"

    # Order matches result.recommendations
    for idx, r in enumerate(result.recommendations):
        row = export_df.iloc[idx]
        assert row["combo"] == r.combo
        if r.is_eligible:
            assert int(row["rank"]) == r.rank
            assert float(row["display_score"]) == r.display_score
        else:
            assert row["rank"] == ""
        assert abs(float(row["overall_score"]) - r.overall_score) < 1e-6
        assert row["preset"] == "community_balanced"
        assert row["scenario_model"] == "hac"
        assert row["scenario_depth"] == "100x"
        assert str(row["scoring_version"]) == "1.0"
        assert row["source_data_commit"] == "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"
        # No mutable branch head URL
        for val in row.values:
            assert "raw.githubusercontent.com" not in str(val)
            assert "/main/" not in str(val)

    # 2. Gate exclusion ordering: ranked eligible first, excluded second with blank rank
    res_gated = score_benchmark(
        df, scenario, weights=weights, gates=GateConfig(complete_recovery=True), preset="complete_replicon_recovery"
    )
    assert res_gated.eligible_count == 7
    assert res_gated.ineligible_count == 10
    gated_csv = export_recommendations_csv(res_gated)
    gated_df = pd.read_csv(io.StringIO(gated_csv), keep_default_na=False)

    assert len(gated_df) == 17
    # First 7 rows must be eligible with integer ranks 1..7
    for i in range(7):
        assert str(gated_df.iloc[i]["is_eligible"]).lower() == "true"
        assert int(gated_df.iloc[i]["rank"]) >= 1
        assert gated_df.iloc[i]["ineligible_reason"] == ""
    # Next 10 rows must be excluded with blank rank and explicit reason
    for i in range(7, 17):
        assert str(gated_df.iloc[i]["is_eligible"]).lower() == "false"
        assert gated_df.iloc[i]["rank"] == ""
        assert "complete-recovery gate" in gated_df.iloc[i]["ineligible_reason"]

    # 3. Invalid weights and invalid scenario prevent producing a misleading export
    with pytest.raises(InvalidWeightsError):
        score_benchmark(df, scenario, weights=WeightsConfig(accuracy=50, contiguity=50, residual=10, replicon=0))

    with pytest.raises(UnsupportedScenarioError):
        score_benchmark(df, Scenario(model="unknown", depth="20x"), weights=weights)


def test_score_benchmark_replicon_denominator_parameterization():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    scenario = Scenario(model="hac", depth="100x")
    weights = WeightsConfig(accuracy=0.0, contiguity=0.0, residual=0.0, replicon=100.0)

    # Base run: denominator 3.0
    res_base = score_benchmark(df, scenario, weights=weights, replicon_denominator=3.0)
    # Denominator 2.0
    res_d2 = score_benchmark(df, scenario, weights=weights, replicon_denominator=2.0)
    # Denominator 4.0
    res_d4 = score_benchmark(df, scenario, weights=weights, replicon_denominator=4.0)

    # Check a combo that has non-zero replicon losses, e.g. chopper-untrimmed or similar
    # All non-replicon criteria must remain identical
    for r_base, r_d2, r_d4 in zip(res_base.recommendations, res_d2.recommendations, res_d4.recommendations):
        assert r_base.combo == r_d2.combo == r_d4.combo
        assert r_base.score_contiguity == r_d2.score_contiguity == r_d4.score_contiguity
        assert r_base.score_accuracy == r_d2.score_accuracy == r_d4.score_accuracy
        assert r_base.score_residual == r_d2.score_residual == r_d4.score_residual

        missed = r_base.replicon_total_missed
        expected_d2 = max(0.0, 1.0 - (missed / 2.0)) * 100.0
        expected_d3 = max(0.0, 1.0 - (missed / 3.0)) * 100.0
        expected_d4 = max(0.0, 1.0 - (missed / 4.0)) * 100.0

        assert abs(r_d2.score_replicon - expected_d2) < 1e-6
        assert abs(r_base.score_replicon - expected_d3) < 1e-6
        assert abs(r_d4.score_replicon - expected_d4) < 1e-6


def test_score_benchmark_expected_samples_parameterization():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    scenario = Scenario(model="hac", depth="100x")
    weights = WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0)

    # Leave one isolate out: omit 'AJ292__202310'
    omitted = "AJ292__202310"
    subset_samples = EXPECTED_SAMPLES - {omitted}
    assert len(subset_samples) == 12

    filtered_df = df[df["sample"] != omitted].copy()

    # Without expected_samples, score_benchmark should treat all combos as having missing isolates
    res_strict = score_benchmark(filtered_df, scenario, weights=weights)
    assert res_strict.eligible_count == 0
    assert all(r.is_insufficient_data for r in res_strict.recommendations)

    # With expected_samples, score_benchmark should accept the 12 isolates as complete
    res_loo = score_benchmark(filtered_df, scenario, weights=weights, expected_samples=subset_samples)
    assert res_loo.eligible_count == 17
    assert not any(r.is_insufficient_data for r in res_loo.recommendations)




