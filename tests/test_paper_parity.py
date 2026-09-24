from pathlib import Path
import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import get_preset, PresetName
from qc_scoring.scorer import score_benchmark


def test_supported_paper_scenarios():
    from qc_scoring.paper import SUPPORTED_PAPER_SCENARIOS
    expected = [
        Scenario("hac", "20x"),
        Scenario("hac", "100x"),
        Scenario("sup", "20x"),
        Scenario("sup", "100x"),
    ]
    assert SUPPORTED_PAPER_SCENARIOS == expected


def test_generate_paper_composite_results_parity_and_fixtures():
    from qc_scoring.paper import generate_paper_composite_results, SUPPORTED_PAPER_SCENARIOS
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    results = generate_paper_composite_results(df)
    assert set(results.keys()) == set(SUPPORTED_PAPER_SCENARIOS)

    expected_leaders = {
        Scenario("hac", "20x"): {
            "combo": "seqkit-dorado",
            "display_score": 94.6,
            "score_accuracy": 80.9,
            "score_contiguity": 99.9,
            "score_residual": 100.0,
            "score_replicon": 100.0,
        },
        Scenario("hac", "100x"): {
            "combo": "seqkit-barbell",
            "display_score": 98.7,
            "score_accuracy": 95.3,
            "score_contiguity": 100.0,
            "score_residual": 100.0,
            "score_replicon": 100.0,
        },
        Scenario("sup", "20x"): {
            "combo": "chopper-porechop_abi",
            "display_score": 97.5,
            "score_accuracy": 91.1,
            "score_contiguity": 100.0,
            "score_residual": 100.0,
            "score_replicon": 100.0,
        },
        Scenario("sup", "100x"): {
            "combo": "chopper-porechop_abi",
            "display_score": 99.3,
            "score_accuracy": 97.5,
            "score_contiguity": 100.0,
            "score_residual": 100.0,
            "score_replicon": 100.0,
        },
    }

    preset = get_preset(PresetName.COMMUNITY_BALANCED)

    for sc, exp in expected_leaders.items():
        res = results[sc]
        assert len(res.recommendations) == 17
        leading = res.leading_recommendation
        assert leading is not None
        assert leading.rank == 1
        assert leading.combo == exp["combo"]
        assert leading.display_score == exp["display_score"]
        assert round(leading.score_accuracy, 1) == exp["score_accuracy"]
        assert round(leading.score_contiguity, 1) == exp["score_contiguity"]
        assert round(leading.score_residual, 1) == exp["score_residual"]
        assert round(leading.score_replicon, 1) == exp["score_replicon"]

        # Assert full precision parity against canonical scorer
        canonical_res = score_benchmark(df, sc, weights=preset.weights, gates=preset.gates)
        for idx in range(17):
            rec = res.recommendations[idx]
            can_rec = canonical_res.recommendations[idx]
            assert rec.combo == can_rec.combo
            assert rec.rank == can_rec.rank
            assert abs(rec.overall_score - can_rec.overall_score) < 1e-9
            assert rec.display_score == can_rec.display_score
            assert abs(rec.score_contiguity - can_rec.score_contiguity) < 1e-9
            assert abs(rec.score_accuracy - can_rec.score_accuracy) < 1e-9
            assert abs(rec.score_residual - can_rec.score_residual) < 1e-9
            assert abs(rec.score_replicon - can_rec.score_replicon) < 1e-9


def test_paper_composite_rejects_unsupported_or_aggregate():
    from qc_scoring.paper import generate_paper_composite_results
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    # If df has unsupported scenarios (e.g. 50x)
    df_invalid = df.copy()
    df_invalid["depth"] = "50x"
    with pytest.raises(ValueError, match="Unsupported scenario values detected"):
        generate_paper_composite_results(df_invalid)

    # If df is missing a required scenario
    df_missing = df.query("depth != '100x'").copy()
    with pytest.raises(ValueError, match="Missing required scenario"):
        generate_paper_composite_results(df_missing)


def test_generate_paper_summary_dataframe_schema_and_values():
    from qc_scoring.paper import generate_paper_summary_dataframe
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    summary_df = generate_paper_summary_dataframe(df)

    # 4 scenarios x 17 combinations = 68 rows
    assert len(summary_df) == 68

    required_columns = [
        "scenario_model",
        "scenario_depth",
        "combo",
        "rank",
        "overall_score",
        "display_score",
        "score_contiguity",
        "score_accuracy",
        "score_residual",
        "score_replicon",
        "mean_auNGA_ratio",
        "mean_error_rate",
        "mean_mismatches",
        "mean_indels",
        "residual_total_hits",
        "residual_clean_isolates",
        "residual_affected_isolates",
        "replicon_total_missed",
        "replicon_full_missed",
        "replicon_partial_missed",
        "replicon_affected_isolates",
        "mean_duplication_ratio",
        "total_misassemblies",
        "is_near_tie",
        "near_tie_label",
        "warnings",
        "is_eligible",
        "ineligible_reason",
        "analysis_description",
        "preset",
        "weight_accuracy",
        "weight_contiguity",
        "weight_residual",
        "weight_replicon",
        "gate_complete_recovery",
        "gate_zero_residual_hits",
        "scoring_version",
        "source_data_commit",
        "source_data_hash",
        "generated_at",
    ]

    for col in required_columns:
        assert col in summary_df.columns, f"Missing column: {col}"

    assert (summary_df["analysis_description"] == "survey-weighted decision analysis (preference alignment)").all()
    assert (summary_df["preset"] == "community_balanced").all()
    assert (summary_df["scoring_version"] == "1.0").all()
    assert (summary_df["weight_accuracy"] == 28.0).all()
    assert (summary_df["weight_contiguity"] == 20.0).all()
    assert (summary_df["weight_residual"] == 17.0).all()
    assert (summary_df["weight_replicon"] == 35.0).all()


def test_export_paper_summary_csv(tmp_path):
    from qc_scoring.paper import export_paper_summary_csv
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    out_file = tmp_path / "paper_summary.csv"
    csv_str = export_paper_summary_csv(df, out_path=out_file)

    assert out_file.exists()
    assert out_file.read_text() == csv_str
    reloaded = pd.read_csv(out_file)
    assert len(reloaded) == 68


def test_plot_paper_composite_figure(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    from qc_scoring.paper import generate_paper_composite_results, plot_paper_composite_figure
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    results = generate_paper_composite_results(df)

    out_hac = tmp_path / "fig_hac.png"
    fig_hac = plot_paper_composite_figure(results, model="hac", out_path=out_hac)
    assert out_hac.exists()
    assert out_hac.stat().st_size > 0

    out_sup = tmp_path / "fig_sup.png"
    fig_sup = plot_paper_composite_figure(results, model="sup", out_path=out_sup)
    assert out_sup.exists()
    assert out_sup.stat().st_size > 0


def test_workflow_script_cli_execution(tmp_path):
    import subprocess
    import sys

    out_summary = tmp_path / "summary.csv"
    out_scores = tmp_path / "scores.csv"
    out_sup = tmp_path / "sup.png"
    out_hac = tmp_path / "hac.png"

    cmd = [
        sys.executable,
        "workflow/scripts/plot_assembly_aggregate_score.py",
        "--master-csv",
        "assets/data/assembly_metrics.csv",
        "--summary-csv-perdepth",
        str(out_summary),
        "--scores-csv-perdepth",
        str(out_scores),
        "--fig-perdepth-sup",
        str(out_sup),
        "--fig-perdepth-hac",
        str(out_hac),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode == 0, f"Script failed: {proc.stderr}\n{proc.stdout}"
    assert out_summary.exists()
    assert len(pd.read_csv(out_summary)) == 68
    assert out_scores.exists()
    assert len(pd.read_csv(out_scores)) == 68
    assert out_sup.exists()
    assert out_hac.exists()


def test_workflow_script_refuses_global_or_aggregate(tmp_path):
    import subprocess
    import sys

    out_global = tmp_path / "global.csv"
    cmd = [
        sys.executable,
        "workflow/scripts/plot_assembly_aggregate_score.py",
        "--master-csv",
        "assets/data/assembly_metrics.csv",
        "--summary-csv-global",
        str(out_global),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    assert proc.returncode != 0, "Script must refuse to generate global/cross-depth outputs"
    assert "cross-depth" in proc.stderr.lower() or "aggregate" in proc.stderr.lower() or "not supported" in proc.stderr.lower() or "refuse" in proc.stderr.lower()


def test_workflow_script_old_formulation_removed():
    content = Path("workflow/scripts/plot_assembly_aggregate_score.py").read_text()
    assert "calculate_scores" not in content, "Old calculate_scores function must be deleted"
    assert "plot_global_figure" not in content, "Old plot_global_figure function must be deleted"
    assert "s_contig_safe" not in content, "Old geometric mean formulation must be deleted"
    assert "qc_scoring.paper" in content, "Must import qc_scoring.paper"

