import pandas as pd
import pytest
from qc_scoring.models import Scenario
from qc_scoring.preferences import PresetName, get_preset
from qc_scoring.robustness import run_denominator_sensitivity


def test_run_denominator_sensitivity():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    scenarios = [
        Scenario(model="hac", depth="20x"),
        Scenario(model="hac", depth="100x"),
        Scenario(model="sup", depth="20x"),
        Scenario(model="sup", depth="100x"),
    ]
    denominators = [2.0, 3.0, 4.0]

    result = run_denominator_sensitivity(df, scenarios=scenarios, denominators=denominators)

    # 4 scenarios * 3 denominators * 17 combinations = 204 records
    assert len(result.records) == 204
    df_res = result.to_dataframe()
    assert len(df_res) == 204

    # Verify all 4 scenarios are present
    models = set(df_res["scenario_model"].unique())
    depths = set(df_res["scenario_depth"].unique())
    assert models == {"hac", "sup"}
    assert depths == {"20x", "100x"}

    # Verify all 3 denominators are present
    assert set(df_res["denominator"].unique()) == {2.0, 3.0, 4.0}

    # Verify all 17 combos are present for each scenario & denominator
    for sc in scenarios:
        for d in denominators:
            sub = df_res[(df_res["scenario_model"] == sc.model) & (df_res["scenario_depth"] == sc.depth) & (df_res["denominator"] == d)]
            assert len(sub) == 17
            assert set(sub["combo"].unique()) == {
                "chopper-barbell", "chopper-dorado", "chopper-porechop_abi", "chopper-untrimmed",
                "fastplong-all", "filtlong-barbell", "filtlong-dorado", "filtlong-porechop_abi",
                "filtlong-untrimmed", "seqkit-barbell", "seqkit-dorado", "seqkit-porechop_abi",
                "seqkit-untrimmed", "unprocessed-barbell", "unprocessed-dorado",
                "unprocessed-porechop_abi", "unprocessed-untrimmed",
            }
            # Every combo must have a rank
            assert sub["rank"].notna().all()
            # Top-ranked combo has rank == 1
            assert (sub["rank"] == 1).any()

    # Verify summary reports winner stability and rank shift without failing on changed winner
    assert len(result.scenario_summaries) == 4
    for summary in result.scenario_summaries:
        assert summary.scenario in scenarios
        assert set(summary.winner_by_denominator.keys()) == {2.0, 3.0, 4.0}
        assert isinstance(summary.winner_is_stable, bool)
        assert isinstance(summary.max_rank_shift, int)
        assert summary.scoring_version == "1.0"
        assert len(summary.source_data_hash) == 64


def test_run_leave_one_isolate_out():
    from qc_scoring.robustness import run_leave_one_isolate_out
    from qc_scoring.validation import EXPECTED_SAMPLES

    df = pd.read_csv("assets/data/assembly_metrics.csv")
    scenarios = [
        Scenario(model="hac", depth="20x"),
        Scenario(model="hac", depth="100x"),
        Scenario(model="sup", depth="20x"),
        Scenario(model="sup", depth="100x"),
    ]

    result = run_leave_one_isolate_out(df, scenarios=scenarios)

    # 4 scenarios * 13 isolates * 17 combinations = 884 records
    assert len(result.records) == 884
    df_res = result.to_dataframe()
    assert len(df_res) == 884

    # Verify all 4 scenarios are present
    assert set(df_res["scenario_model"].unique()) == {"hac", "sup"}
    assert set(df_res["scenario_depth"].unique()) == {"20x", "100x"}

    # Verify all 13 isolates are omitted in turn
    assert set(df_res["omitted_isolate"].unique()) == EXPECTED_SAMPLES

    # Verify all 17 combos for each scenario & omitted isolate
    for sc in scenarios:
        for isolate in EXPECTED_SAMPLES:
            sub = df_res[
                (df_res["scenario_model"] == sc.model)
                & (df_res["scenario_depth"] == sc.depth)
                & (df_res["omitted_isolate"] == isolate)
            ]
            assert len(sub) == 17
            assert sub["rank"].notna().all()
            assert (sub["rank"] == 1).any()

    # Verify summary reports winner stability and rank shift
    assert len(result.scenario_summaries) == 4
    for summary in result.scenario_summaries:
        assert summary.scenario in scenarios
        assert len(summary.winner_by_omitted_isolate) == 13
        assert isinstance(summary.winner_is_stable, bool)
        assert 0.0 <= summary.winner_stability_frequency <= 1.0
        assert isinstance(summary.max_rank_shift, int)
        assert len(summary.combo_rank_ranges) == 17
        assert summary.scoring_version == "1.0"
        assert len(summary.source_data_hash) == 64


def test_run_respondent_profiles():
    from qc_scoring.robustness import run_respondent_profiles

    df = pd.read_csv("assets/data/assembly_metrics.csv")
    survey_df = pd.read_csv("microbial-qc-survey.csv")
    assert len(survey_df) == 22, "Survey must have 22 observed respondent profiles"

    scenarios = [
        Scenario(model="hac", depth="20x"),
        Scenario(model="hac", depth="100x"),
        Scenario(model="sup", depth="20x"),
        Scenario(model="sup", depth="100x"),
    ]

    result = run_respondent_profiles(df, survey_df, scenarios=scenarios)

    # 4 scenarios * 22 respondents * 17 combinations = 1496 records
    assert len(result.records) == 1496
    df_res = result.to_dataframe()
    assert len(df_res) == 1496

    # Verify all 4 scenarios are present
    assert set(df_res["scenario_model"].unique()) == {"hac", "sup"}
    assert set(df_res["scenario_depth"].unique()) == {"20x", "100x"}

    # Verify all 22 respondents are represented
    assert set(df_res["respondent_id"].unique()) == set(range(1, 23))

    # Verify each scenario has summaries with winner frequencies and rank variation
    assert len(result.scenario_summaries) == 4
    for summary in result.scenario_summaries:
        assert summary.scenario in scenarios
        assert summary.total_respondents == 22
        # Sum of winner frequencies across combos must equal 1.0 (or count sum == 22)
        total_winner_counts = sum(summary.winner_counts.values())
        assert total_winner_counts == 22
        total_freq = sum(summary.winner_frequencies.values())
        assert abs(total_freq - 1.0) < 1e-6

        # Rank variation recorded for all 17 combos
        assert len(summary.combo_rank_variation) == 17
        for combo, var in summary.combo_rank_variation.items():
            assert "min_rank" in var
            assert "max_rank" in var
            assert "mean_rank" in var
            assert 1 <= var["min_rank"] <= var["max_rank"] <= 17

        # Verify descriptive recording with provenance
        assert summary.scoring_version == "1.0"
        assert len(summary.source_data_hash) == 64
        assert len(summary.survey_data_hash) == 64


def test_generate_robustness_evidence_and_acceptance_criteria(tmp_path):
    import json
    from qc_scoring.robustness import generate_robustness_evidence

    benchmark_csv = "assets/data/assembly_metrics.csv"
    survey_csv = "microbial-qc-survey.csv"
    out_dir = tmp_path / "robustness"

    bundle = generate_robustness_evidence(
        benchmark_path=benchmark_csv,
        survey_path=survey_csv,
        output_dir=out_dir,
    )

    # Verify machine-readable outputs exist
    denom_csv = out_dir / "denominator_sensitivity.csv"
    loo_csv = out_dir / "leave_one_out.csv"
    resp_csv = out_dir / "respondent_profiles.csv"
    summary_json = out_dir / "evidence_summary.json"
    summary_md = out_dir / "summary.md"

    assert denom_csv.exists()
    assert loo_csv.exists()
    assert resp_csv.exists()
    assert summary_json.exists()
    assert summary_md.exists()

    # Verify denominator sensitivity coverage
    df_denom = pd.read_csv(denom_csv)
    assert len(df_denom) == 204
    assert set(df_denom["denominator"].unique()) == {2.0, 3.0, 4.0}
    assert set(df_denom["scenario_model"].unique()) == {"hac", "sup"}
    assert set(df_denom["scenario_depth"].unique()) == {"20x", "100x"}

    # Verify leave-one-out coverage (all 13 cases across 4 scenarios)
    df_loo = pd.read_csv(loo_csv)
    assert len(df_loo) == 884
    assert len(df_loo["omitted_isolate"].unique()) == 13
    assert set(df_loo["scenario_model"].unique()) == {"hac", "sup"}
    assert set(df_loo["scenario_depth"].unique()) == {"20x", "100x"}

    # Verify respondent profile coverage (all 22 profiles across 4 scenarios)
    df_resp = pd.read_csv(resp_csv)
    assert len(df_resp) == 1496
    assert len(df_resp["respondent_id"].unique()) == 22
    assert set(df_resp["scenario_model"].unique()) == {"hac", "sup"}
    assert set(df_resp["scenario_depth"].unique()) == {"20x", "100x"}

    # Verify JSON structure
    with open(summary_json) as f:
        meta = json.load(f)
    assert meta["scoring_version"] == "1.0"
    assert "source_data_hash" in meta
    assert "survey_data_hash" in meta
    assert len(meta["denominator_sensitivity"]) == 4
    assert len(meta["leave_one_out"]) == 4
    assert len(meta["respondent_profiles"]) == 4

    # Verify human-readable summary contains required descriptive elements
    md_content = summary_md.read_text()
    assert "Scoring Version 1.0" in md_content or "scoring-v1" in md_content.lower()
    assert "22 observed responses" in md_content or "22 observed respondent" in md_content
    assert "wider population" in md_content  # explains no inference made to wider population
    assert "superseded" in md_content or "geometric" in md_content  # confirms not compared to superseded
    assert "25" in md_content  # notes 25-point threshold unchanged
    assert "denominator" in md_content.lower()
    assert "leave-one" in md_content.lower()



