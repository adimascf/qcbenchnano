import hashlib
import io
import re
import socket
import subprocess
import time
from pathlib import Path
import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import get_preset, PresetName, WeightsConfig, GateConfig
from qc_scoring.scorer import score_benchmark, export_recommendations_csv
from qc_scoring.config import load_dashboard_config


def _extract_app_code(qmd_path: Path = Path("tool_weighting.qmd")) -> str:
    content = qmd_path.read_text()
    match = re.search(r"```{shinylive-python}(.*?)```", content, re.DOTALL)
    assert match is not None, "shinylive-python block not found in tool_weighting.qmd"
    code = match.group(1)
    lines = [l for l in code.splitlines() if not l.strip().startswith("#|")]
    return "\n".join(lines)


# 1. Source contract checks for Issue #17 in tool_weighting.qmd
def test_dashboard_source_csv_export_contract():
    qmd_content = Path("tool_weighting.qmd").read_text()

    # Must contain download_recommendations_csv identifier
    assert "download_recommendations_csv" in qmd_content, (
        "tool_weighting.qmd must define download_recommendations_csv"
    )

    # Must contain generate_dashboard_csv_export helper
    assert "generate_dashboard_csv_export" in qmd_content, (
        "tool_weighting.qmd must define generate_dashboard_csv_export"
    )

    # Must have shiny-download-link in shortlist content
    assert "shiny-download-link" in qmd_content, (
        "tool_weighting.qmd must render a shiny-download-link for CSV download"
    )

    # Must not contain mutable branch head URL as evidence provenance
    assert "raw.githubusercontent.com" not in qmd_content or "main/logbook" not in qmd_content, (
        "Evidence provenance must not embed mutable branch head URLs"
    )


# 2. Complete parity check between dashboard CSV export and canonical scorer export
@pytest.mark.parametrize(
    "model,depth,preset_enum",
    [
        ("hac", "20x", PresetName.COMMUNITY_BALANCED),
        ("hac", "100x", PresetName.COMMUNITY_BALANCED),
        ("sup", "20x", PresetName.COMMUNITY_BALANCED),
        ("sup", "100x", PresetName.COMMUNITY_BALANCED),
        ("hac", "100x", PresetName.COMPLETE_REPLICON_RECOVERY),
        ("sup", "100x", PresetName.SEQUENCE_ACCURATE_ASSEMBLY),
    ]
)
def test_csv_export_parity_with_canonical_scorer(model, depth, preset_enum):
    clean_code = _extract_app_code()
    scope = {}
    exec(clean_code, scope)

    generate_csv_func = scope["generate_dashboard_csv_export"]
    app_config = scope["app_config"]

    preset = get_preset(preset_enum)
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    # 1. Canonical scorer export
    canonical_result = score_benchmark(
        df=df,
        scenario=Scenario(model=model, depth=depth),
        weights=preset.weights,
        gates=preset.gates,
        preset=preset_enum.value,
    )
    canonical_csv_str = export_recommendations_csv(canonical_result)
    canonical_df = pd.read_csv(io.StringIO(canonical_csv_str), keep_default_na=False)

    # 2. Dashboard helper export
    dash_csv_str = generate_csv_func(
        df=df,
        model=model,
        depth=depth,
        preset=preset_enum.value,
        weight_accuracy=preset.weights.accuracy,
        weight_contiguity=preset.weights.contiguity,
        weight_residual=preset.weights.residual,
        weight_replicon=preset.weights.replicon,
        gate_complete=preset.gates.complete_recovery,
        gate_zero_hits=preset.gates.zero_residual_hits,
        config=app_config,
    )
    dash_df = pd.read_csv(io.StringIO(dash_csv_str), keep_default_na=False)

    # Assert row count parity (all 17 combinations)
    assert len(dash_df) == 17
    assert len(canonical_df) == 17

    # Assert exact column set
    assert list(dash_df.columns) == list(canonical_df.columns)

    # Assert exact row ordering and value equivalence across all combinations
    for idx in range(17):
        dash_row = dash_df.iloc[idx]
        canon_row = canonical_df.iloc[idx]

        assert dash_row["combo"] == canon_row["combo"]
        assert str(dash_row["rank"]) == str(canon_row["rank"])
        assert str(dash_row["is_eligible"]) == str(canon_row["is_eligible"])
        assert str(dash_row["is_insufficient_data"]) == str(canon_row["is_insufficient_data"])
        assert dash_row["ineligible_reason"] == canon_row["ineligible_reason"]
        assert str(dash_row["is_near_tie"]) == str(canon_row["is_near_tie"])
        assert dash_row["near_tie_label"] == canon_row["near_tie_label"]

        if dash_row["display_score"] != "":
            assert float(dash_row["display_score"]) == pytest.approx(float(canon_row["display_score"]), abs=1e-4)
            assert float(dash_row["overall_score"]) == pytest.approx(float(canon_row["overall_score"]), abs=1e-6)

        assert float(dash_row["score_accuracy"]) == pytest.approx(float(canon_row["score_accuracy"]), abs=1e-4)
        assert float(dash_row["score_contiguity"]) == pytest.approx(float(canon_row["score_contiguity"]), abs=1e-4)
        assert float(dash_row["score_residual"]) == pytest.approx(float(canon_row["score_residual"]), abs=1e-4)
        assert float(dash_row["score_replicon"]) == pytest.approx(float(canon_row["score_replicon"]), abs=1e-4)

        assert int(dash_row["residual_total_hits"]) == int(canon_row["residual_total_hits"])
        assert int(dash_row["residual_clean_isolates"]) == int(canon_row["residual_clean_isolates"])
        assert int(dash_row["residual_affected_isolates"]) == int(canon_row["residual_affected_isolates"])
        assert int(dash_row["replicon_total_missed"]) == int(canon_row["replicon_total_missed"])
        assert int(dash_row["replicon_full_missed"]) == int(canon_row["replicon_full_missed"])
        assert int(dash_row["replicon_partial_missed"]) == int(canon_row["replicon_partial_missed"])
        assert int(dash_row["replicon_affected_isolates"]) == int(canon_row["replicon_affected_isolates"])

        assert dash_row["warnings"] == canon_row["warnings"]
        assert dash_row["preset"] == canon_row["preset"]
        assert dash_row["scenario_model"] == canon_row["scenario_model"]
        assert dash_row["scenario_depth"] == canon_row["scenario_depth"]
        assert str(dash_row["scoring_version"]) == "1.0"
        assert dash_row["source_data_commit"] == "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"
        assert dash_row["source_data_hash"] == canon_row["source_data_hash"]


# 3. Guard against misleading downloads on invalid weights or unselected scenario
def test_dashboard_csv_export_validation_guards():
    clean_code = _extract_app_code()
    scope = {}
    exec(clean_code, scope)

    generate_csv_func = scope["generate_dashboard_csv_export"]
    app_config = scope["app_config"]
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    # Unselected model or depth
    with pytest.raises(ValueError, match="Scenario must be selected"):
        generate_csv_func(
            df=df,
            model="",
            depth="100x",
            preset="community_balanced",
            weight_accuracy=28,
            weight_contiguity=20,
            weight_residual=17,
            weight_replicon=35,
            gate_complete=False,
            gate_zero_hits=False,
            config=app_config,
        )

    with pytest.raises(ValueError, match="Scenario must be selected"):
        generate_csv_func(
            df=df,
            model="hac",
            depth="",
            preset="community_balanced",
            weight_accuracy=28,
            weight_contiguity=20,
            weight_residual=17,
            weight_replicon=35,
            gate_complete=False,
            gate_zero_hits=False,
            config=app_config,
        )

    # Invalid weights sum (not 100%)
    with pytest.raises(ValueError, match="total exactly 100%"):
        generate_csv_func(
            df=df,
            model="hac",
            depth="100x",
            preset="custom",
            weight_accuracy=50,
            weight_contiguity=20,
            weight_residual=20,
            weight_replicon=20,  # Sum = 110%
            gate_complete=False,
            gate_zero_hits=False,
            config=app_config,
        )

    # Negative weights
    with pytest.raises(ValueError, match="non-negative"):
        generate_csv_func(
            df=df,
            model="hac",
            depth="100x",
            preset="custom",
            weight_accuracy=-10,
            weight_contiguity=40,
            weight_residual=40,
            weight_replicon=30,
            gate_complete=False,
            gate_zero_hits=False,
            config=app_config,
        )


# 4. Playwright End-to-End Browser Journey Test
def test_browser_journey_csv_export(tmp_path):
    from playwright.sync_api import sync_playwright

    qmd_path = Path("tool_weighting.qmd")
    app_code = _extract_app_code(qmd_path)
    app_file = tmp_path / "app.py"
    app_file.write_text(app_code)

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("", 0))
    port = sock.getsockname()[1]
    sock.close()

    proc = subprocess.Popen(
        ["shiny", "run", "--port", str(port), str(app_file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    try:
        url = f"http://127.0.0.1:{port}"
        started = False
        for _ in range(50):
            try:
                import urllib.request
                urllib.request.urlopen(url, timeout=1)
                started = True
                break
            except Exception:
                time.sleep(0.1)

        assert started, f"Shiny server failed to start: {proc.stderr.read()}"

        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True)
            page = browser.new_page()
            page.goto(url)

            # 1. Unselected scenario: download button is NOT visible in shortlist
            assert page.locator("#download_recommendations_csv").count() == 0, (
                "Download button should not be rendered when scenario is unset"
            )

            # 2. Select scenario: HAC at 100x
            page.select_option("#model_select", "hac")
            page.select_option("#depth_select", "100x")

            page.wait_for_selector("#download_recommendations_csv")
            page.wait_for_selector(".ranking-table")

            # 3. Download CSV for Community-balanced
            with page.expect_download() as download_info:
                page.click("#download_recommendations_csv")
            download = download_info.value
            csv_path = download.path()
            csv_content = Path(csv_path).read_text()
            df_cb = pd.read_csv(io.StringIO(csv_content), keep_default_na=False)

            # Verify all 17 combinations in exact displayed order
            assert len(df_cb) == 17

            # Verify all 40 required columns
            expected_fields = [
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
            for field in expected_fields:
                assert field in df_cb.columns, f"Downloaded CSV missing field: {field}"

            # Canonical records comparison across all 17 rows
            canonical_benchmark = score_benchmark(
                df=pd.read_csv("assets/data/assembly_metrics.csv"),
                scenario=Scenario(model="hac", depth="100x"),
                weights=get_preset(PresetName.COMMUNITY_BALANCED).weights,
                preset="community_balanced",
            )
            for idx, canon_rec in enumerate(canonical_benchmark.recommendations):
                row = df_cb.iloc[idx]
                assert row["combo"] == canon_rec.combo
                assert int(row["rank"]) == canon_rec.rank
                assert float(row["display_score"]) == pytest.approx(canon_rec.display_score, abs=1e-4)
                assert float(row["overall_score"]) == pytest.approx(canon_rec.overall_score, abs=1e-6)
                assert float(row["score_accuracy"]) == pytest.approx(canon_rec.score_accuracy, abs=1e-4)
                assert float(row["score_contiguity"]) == pytest.approx(canon_rec.score_contiguity, abs=1e-4)
                assert float(row["score_residual"]) == pytest.approx(canon_rec.score_residual, abs=1e-4)
                assert float(row["score_replicon"]) == pytest.approx(canon_rec.score_replicon, abs=1e-4)
                assert int(row["residual_total_hits"]) == canon_rec.residual_total_hits
                assert int(row["replicon_total_missed"]) == canon_rec.replicon_total_missed
                assert row["warnings"] == "; ".join(canon_rec.warnings)

            assert df_cb.iloc[0]["combo"] == "seqkit-barbell"
            assert int(df_cb.iloc[0]["rank"]) == 1
            assert float(df_cb.iloc[0]["display_score"]) == 98.7
            assert df_cb.iloc[0]["preset"] == "community_balanced"
            assert df_cb.iloc[0]["scenario_model"] == "hac"
            assert df_cb.iloc[0]["scenario_depth"] == "100x"
            assert str(df_cb.iloc[0]["scoring_version"]) == "1.0"
            assert df_cb.iloc[0]["source_data_commit"] == "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"
            assert "raw.githubusercontent.com" not in csv_content
            assert "/main/" not in csv_content

            # Verify displayed table rows match CSV rows exactly
            table_combos = page.locator(".ranking-table tbody tr td:nth-child(2)").all_inner_texts()
            for idx, combo_name in enumerate(table_combos):
                assert df_cb.iloc[idx]["combo"] == combo_name

            # 4. Preset change: Complete-replicon-recovery
            page.select_option("#preset_select", "complete_replicon_recovery")
            page.wait_for_timeout(500)

            with page.expect_download() as download_info_cr:
                page.click("#download_recommendations_csv")
            download_cr = download_info_cr.value
            csv_cr_content = Path(download_cr.path()).read_text()
            df_cr = pd.read_csv(io.StringIO(csv_cr_content), keep_default_na=False)

            assert len(df_cr) == 17
            # In HAC 100x with complete replicon gate: 7 eligible, 10 excluded
            assert (df_cr["is_eligible"] == True).sum() == 7
            assert (df_cr["is_eligible"] == False).sum() == 10

            # Eligible combinations appear first with ranks
            for i in range(7):
                assert int(df_cr.iloc[i]["rank"]) >= 1
                assert df_cr.iloc[i]["ineligible_reason"] == ""

            # Excluded follow with blank ranks and explicit reasons
            for i in range(7, 17):
                assert df_cr.iloc[i]["rank"] == ""
                assert "complete-recovery gate" in df_cr.iloc[i]["ineligible_reason"]

            # 5. Valid Custom Weighting & Gate change
            # Reset to Community-balanced, then modify custom weights
            page.click("#btn_reset")
            page.wait_for_timeout(400)
            page.fill("#w_acc", "40")
            page.fill("#w_cont", "20")
            page.fill("#w_res", "20")
            page.fill("#w_rep", "20")
            page.wait_for_selector("#weight_total_indicator:has-text('100%')")
            page.wait_for_timeout(500)

            # Preset selector should automatically switch to Custom
            assert page.input_value("#preset_select") == "custom"

            with page.expect_download() as download_info_custom:
                page.click("#download_recommendations_csv")
            download_custom = download_info_custom.value
            csv_custom_content = Path(download_custom.path()).read_text()
            df_custom = pd.read_csv(io.StringIO(csv_custom_content), keep_default_na=False)

            assert len(df_custom) == 17
            assert df_custom.iloc[0]["preset"] == "custom"
            assert float(df_custom.iloc[0]["weight_accuracy"]) == 40.0
            assert float(df_custom.iloc[0]["weight_contiguity"]) == 20.0
            assert float(df_custom.iloc[0]["weight_residual"]) == 20.0
            assert float(df_custom.iloc[0]["weight_replicon"]) == 20.0

            # 6. Reveal-all interaction
            page.click("#btn_reveal_all")
            page.wait_for_selector(".ranking-table tbody tr:nth-child(17)")

            with page.expect_download() as download_info_revealed:
                page.click("#download_recommendations_csv")
            download_revealed = download_info_revealed.value
            csv_revealed_content = Path(download_revealed.path()).read_text()
            df_revealed = pd.read_csv(io.StringIO(csv_revealed_content), keep_default_na=False)

            assert len(df_revealed) == 17
            revealed_table_combos = page.locator(".ranking-table tbody tr td:nth-child(2)").all_inner_texts()
            assert len(revealed_table_combos) == 17
            for idx, combo_name in enumerate(revealed_table_combos):
                assert df_revealed.iloc[idx]["combo"] == combo_name

            # 7. Invalid weights prevent misleading download
            page.fill("#w_acc", "90")  # Total becomes 150%
            page.wait_for_selector("#weight_total_indicator:has-text('above 100%')")
            page.wait_for_timeout(400)

            # Shortlist content shows paused message and download button is NOT shown
            assert page.locator("#download_recommendations_csv").count() == 0, (
                "Download button must not be available when weights are invalid"
            )
            assert "Ranking paused" in page.locator("#region-shortlist").text_content()

            browser.close()
    finally:
        proc.terminate()
        proc.wait()
