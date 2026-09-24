import re
import subprocess
import time
from pathlib import Path
import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import get_preset, PresetName, WeightsConfig, GateConfig
from qc_scoring.scorer import score_benchmark
from qc_scoring.config import load_dashboard_config
def _extract_app_code(qmd_path: Path = Path("tool_weighting.qmd")) -> str:
    qmd_content = qmd_path.read_text()
    match = re.search(r"```\{shinylive-python\}\n(.*?)\n```", qmd_content, re.DOTALL)
    assert match is not None, "shinylive-python code chunk not found in dashboard source"
    app_lines = [l for l in match.group(1).split("\n") if not l.strip().startswith("#|")]
    return "\n".join(app_lines)



# 1. Source contract checks for presets, custom weights, and gates in tool_weighting.qmd
def test_dashboard_source_preset_and_gates_specifications():
    qmd_content = Path("tool_weighting.qmd").read_text()

    # Presets dropdown must contain all specified presets
    assert '"community_balanced":' in qmd_content
    assert '"complete_replicon_recovery":' in qmd_content
    assert '"sequence_accurate_assembly":' in qmd_content
    assert '"custom":' in qmd_content

    # Reset action button must exist
    assert "btn_reset" in qmd_content

    # Numeric inputs must support decimals (step="0.1" or step="any")
    for field in ["w_acc", "w_cont", "w_res", "w_rep"]:
        match = re.search(rf'ui\.input_numeric\(\s*"{field}".*?step=([0-9.]+)', qmd_content)
        assert match is not None, f"Input {field} should define a step attribute"
        step_val = float(match.group(1))
        assert step_val < 1.0, f"Input {field} must permit decimal steps, got {step_val}"

    # Non-compensatory gates must be present
    assert "gate_complete" in qmd_content
    assert "gate_zero_hits" in qmd_content


# 2. Complete parity for all presets and scenarios between dashboard calculate_scenario_rankings and canonical scorer
def test_dashboard_evaluated_results_all_presets_and_gates_parity():
    clean_code = _extract_app_code()
    scope = {}
    exec(clean_code, scope)

    calc_func = scope["calculate_scenario_rankings"]
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    test_presets = [
        PresetName.COMMUNITY_BALANCED,
        PresetName.COMPLETE_REPLICON_RECOVERY,
        PresetName.SEQUENCE_ACCURATE_ASSEMBLY,
    ]

    for p_name in test_presets:
        preset = get_preset(p_name)
        for model in ["hac", "sup"]:
            for depth in ["20x", "100x"]:
                sc = Scenario(model=model, depth=depth)
                can_res = score_benchmark(df, sc, weights=preset.weights, gates=preset.gates)
                res = calc_func(
                    df,
                    model=model,
                    depth=depth,
                    weight_accuracy=preset.weights.accuracy,
                    weight_contiguity=preset.weights.contiguity,
                    weight_residual=preset.weights.residual,
                    weight_replicon=preset.weights.replicon,
                    gate_complete=preset.gates.complete_recovery,
                    gate_zero_hits=preset.gates.zero_residual_hits,
                )
                dash_eligible, dash_excluded, *rest = res
                dash_insufficient = rest[0] if rest else []

                can_eligible = [r for r in can_res.recommendations if r.is_eligible]
                can_excluded = [r for r in can_res.recommendations if not r.is_eligible]

                assert len(dash_eligible) == len(can_eligible), (
                    f"Eligible count mismatch for {p_name} {model} {depth}: "
                    f"dashboard={len(dash_eligible)}, canonical={len(can_eligible)}"
                )
                assert len(dash_excluded) == len(can_excluded), (
                    f"Excluded count mismatch for {p_name} {model} {depth}: "
                    f"dashboard={len(dash_excluded)}, canonical={len(can_excluded)}"
                )

                for idx in range(len(can_eligible)):
                    can_rec = can_eligible[idx]
                    dash_rec = dash_eligible[idx]
                    assert dash_rec["combo"] == can_rec.combo
                    assert dash_rec["rank"] == can_rec.rank
                    assert abs(dash_rec["overall_score"] - can_rec.overall_score) < 1e-5
                    assert dash_rec["display_score"] == can_rec.display_score


# 3. Edge case parity: All-ineligible scenario and missing isolate insufficient data
def test_dashboard_all_ineligible_and_insufficient_data():
    clean_code = _extract_app_code()
    scope = {}
    exec(clean_code, scope)
    calc_func = scope["calculate_scenario_rankings"]

    df = pd.read_csv("assets/data/assembly_metrics.csv")

    # Force all-ineligible by creating synthetic df where all combos fail gates
    df_failing = df.copy()
    df_failing["contamination_count"] = 5  # Every row has residual hits

    res = calc_func(
        df_failing,
        model="hac",
        depth="100x",
        weight_accuracy=28,
        weight_contiguity=20,
        weight_residual=17,
        weight_replicon=35,
        gate_complete=False,
        gate_zero_hits=True,
    )
    eligible, excluded, *rest = res
    assert len(eligible) == 0, "No combinations should be eligible when all fail zero-residual-hits gate"
    assert len(excluded) == 17, "All 17 combinations should appear in excluded list"
    for ex in excluded:
        assert ex.get("rank") is None
        assert any("fails zero-residual-hits gate" in r for r in ex["reasons"])

    # Test missing isolate insufficient data
    mask = ~((df["combo"] == "chopper-barbell") & (df["sample"] == "AJ292__202310"))
    df_missing = df[mask].copy()

    res_missing = calc_func(
        df_missing,
        model="hac",
        depth="100x",
        weight_accuracy=28,
        weight_contiguity=20,
        weight_residual=17,
        weight_replicon=35,
        gate_complete=False,
        gate_zero_hits=False,
    )
    eligible_missing, excluded_missing, *insufficient_missing_rest = res_missing
    insufficient_missing = insufficient_missing_rest[0] if insufficient_missing_rest else [r for r in excluded_missing if r.get("is_insufficient_data")]
    assert any(r["combo"] == "chopper-barbell" for r in insufficient_missing)
    chopper = [r for r in insufficient_missing if r["combo"] == "chopper-barbell"][0]
    assert chopper.get("rank") is None
    assert "insufficient benchmark data" in (chopper.get("ineligible_reason") or " ".join(chopper.get("reasons", [])))


# 4. Playwright Browser-Level Behavioral Journey for Presets, Custom Weights, and Gates
def test_browser_acceptance_presets_custom_and_gates_journey(tmp_path):
    from playwright.sync_api import sync_playwright
    import socket

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

            page.wait_for_selector("#region-scenario")

            # 1. Select scenario: HAC + 100x
            page.select_option("#model_select", "hac")
            page.select_option("#depth_select", "100x")
            page.wait_for_selector("#region-shortlist table tbody tr:first-child:has-text('seqkit-barbell')")

            # Initial state is Community-balanced: 28, 20, 17, 35, both gates unchecked
            assert page.input_value("#w_acc") in ("28", "28.0")
            assert page.input_value("#w_cont") in ("20", "20.0")
            assert page.input_value("#w_res") in ("17", "17.0")
            assert page.input_value("#w_rep") in ("35", "35.0")
            assert not page.is_checked("#gate_complete")
            assert not page.is_checked("#gate_zero_hits")

            # 2. Select Complete-replicon-recovery preset
            page.select_option("#preset_select", "complete_replicon_recovery")
            page.wait_for_function("document.querySelector('#w_acc').value == '43' || document.querySelector('#w_acc').value == '43.0'")
            assert page.input_value("#w_acc") in ("43", "43.0")
            assert page.input_value("#w_cont") in ("31", "31.0")
            assert page.input_value("#w_res") in ("26", "26.0")
            assert page.input_value("#w_rep") in ("0", "0.0")
            assert page.is_checked("#gate_complete")
            assert not page.is_checked("#gate_zero_hits")

            # In Complete-replicon-recovery (HAC 100x), exactly 7 combinations are eligible
            page.wait_for_selector("#region-shortlist table tbody tr:first-child:has-text('98.0')")
            page.click("#btn_reveal_all")
            page.wait_for_function("document.querySelectorAll('#region-shortlist table tbody tr').length === 7")
            assert page.locator("#region-shortlist table tbody tr").count() == 7

            # Excluded section lists the 10 failed combinations with reasons
            excluded_text = page.text_content("#region-shortlist details")
            assert "Excluded / Ineligible Combinations (10)" in excluded_text
            assert "fails complete-recovery gate" in excluded_text

            # 3. Select Sequence-accurate-assembly preset
            page.select_option("#preset_select", "sequence_accurate_assembly")
            page.wait_for_function("document.querySelector('#w_acc').value == '50' || document.querySelector('#w_acc').value == '50.0'")
            assert page.input_value("#w_acc") in ("50", "50.0")
            assert page.input_value("#w_cont") in ("14", "14.0")
            assert page.input_value("#w_res") in ("12", "12.0")
            assert page.input_value("#w_rep") in ("24", "24.0")
            assert not page.is_checked("#gate_complete")
            assert not page.is_checked("#gate_zero_hits")

            # In Sequence-accurate-assembly, all 17 combinations are eligible; leader score is 97.6
            page.wait_for_selector("#region-shortlist table tbody tr:first-child:has-text('97.6')")

            # 4. Edit a weight visibly changes preset to "custom" and pauses ranking when sum != 100
            page.fill("#w_acc", "45")
            # Preset should visibly switch to custom
            page.wait_for_function("document.querySelector('#preset_select').value === 'custom'")
            assert page.input_value("#preset_select") == "custom"

            # Sum is now 45 + 14 + 12 + 24 = 95. Indicator shows 5% remaining and ranking paused
            page.wait_for_selector("#region-priorities :has-text('remaining to allocate')")
            indicator_text = page.text_content("#region-priorities")
            assert "remaining to allocate" in indicator_text
            assert "Ranking paused" in indicator_text

            # Ranking table must pause / disappear
            page.wait_for_selector("#region-shortlist .alert-warning")
            warning_text = page.text_content("#region-shortlist")
            assert "Ranking paused" in warning_text
            assert page.locator("#region-shortlist table tbody tr").count() == 0

            # Restore total to 100 with decimals in custom mode: 45.0 + 19.0 + 12.0 + 24.0 = 100.0
            page.fill("#w_cont", "19")
            page.wait_for_function("document.querySelector('#region-priorities').textContent.includes('100%')")
            page.wait_for_selector("#region-shortlist table tbody tr:first-child")
            assert page.locator("#region-shortlist table tbody tr").count() > 0

            # 5. Edit a gate visibly changes preset to "custom"
            # First reset to Community-balanced
            page.click("#btn_reset")
            page.wait_for_function("document.querySelector('#preset_select').value === 'community_balanced'")
            assert page.input_value("#preset_select") == "community_balanced"
            assert page.input_value("#w_acc") in ("28", "28.0")
            assert not page.is_checked("#gate_complete")

            # Toggle a gate
            page.check("#gate_zero_hits")
            page.wait_for_function("document.querySelector('#preset_select').value === 'custom'")
            assert page.input_value("#preset_select") == "custom"

            # 6. Reset restores Community-balanced
            page.click("#btn_reset")
            page.wait_for_function("document.querySelector('#preset_select').value === 'community_balanced'")
            assert page.input_value("#preset_select") == "community_balanced"
            assert page.input_value("#w_acc") in ("28", "28.0")
            assert page.input_value("#w_cont") in ("20", "20.0")
            assert page.input_value("#w_res") in ("17", "17.0")
            assert page.input_value("#w_rep") in ("35", "35.0")
            assert not page.is_checked("#gate_complete")
            assert not page.is_checked("#gate_zero_hits")

            # 7. Multiple exclusion reasons
            page.check("#gate_complete")
            page.check("#gate_zero_hits")
            # Wait for excluded list to update to 12 combinations and assert multiple exclusion reasons appear
            page.wait_for_selector("#region-shortlist details summary:has-text('(12)')")
            details_content = page.text_content("#region-shortlist details")
            assert "fails complete-recovery gate" in details_content
            assert "fails zero-residual-hits gate" in details_content
            # Specifically check that unprocessed-dorado has both exclusion reasons
            dorado_li = page.locator("#region-shortlist details li:has-text('unprocessed-dorado')").text_content()
            assert "fails complete-recovery gate" in dorado_li
            assert "fails zero-residual-hits gate" in dorado_li

            browser.close()
    finally:
        proc.terminate()
        proc.wait()


def test_browser_acceptance_all_ineligible_state(tmp_path):
    from playwright.sync_api import sync_playwright
    import socket

    # Create app copy pointing to synthetic metrics file where all combos have contamination_count = 5
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    df_all_failing = df.copy()
    df_all_failing["contamination_count"] = 5
    failing_csv = tmp_path / "assembly_metrics.csv"
    df_all_failing.to_csv(failing_csv, index=False)

    qmd_path = Path("tool_weighting.qmd")
    app_code = _extract_app_code(qmd_path)
    # Replace metrics path to load the failing csv
    app_code = app_code.replace('"assets/data/assembly_metrics.csv"', f'"{failing_csv}"')

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

            # Select scenario
            page.select_option("#model_select", "hac")
            page.select_option("#depth_select", "100x")

            # Enable zero residual hits gate (which will cause all 17 combos to fail)
            page.check("#gate_zero_hits")

            # Wait for all-ineligible warning message
            page.wait_for_selector("#no-eligible-msg")
            msg = page.text_content("#no-eligible-msg")
            assert "No evaluated combination meets the active eligibility requirements" in msg
            assert "All gates have been preserved" in msg

            # Ensure shortlist table has 0 rows and no failed combination is promoted
            assert page.locator("#region-shortlist table tbody tr").count() == 0

            # Ensure gates remain checked
            assert page.is_checked("#gate_zero_hits")

            # Detail section reports that no combination meets requirements
            detail_msg = page.text_content("#region-detail")
            assert "No combination meets active eligibility requirements" in detail_msg

            # Excluded section shows all 17 failed combinations
            page.wait_for_selector("#region-shortlist details summary:has-text('(17)')")
            summary_text = page.text_content("#region-shortlist details summary")
            assert "(17)" in summary_text

            browser.close()
    finally:
        proc.terminate()
        proc.wait()

