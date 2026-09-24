import re
import subprocess
import time
from pathlib import Path
import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import get_preset, PresetName, WeightsConfig, GateConfig
from qc_scoring.scorer import score_benchmark


def _extract_app_code(qmd_path: Path = Path("tool_weighting.qmd")) -> str:
    qmd_content = qmd_path.read_text()
    match = re.search(r"```\{shinylive-python\}\n(.*?)\n```", qmd_content, re.DOTALL)
    assert match is not None, "shinylive-python code chunk not found in dashboard source"
    app_lines = [l for l in match.group(1).split("\n") if not l.strip().startswith("#|")]
    return "\n".join(app_lines)


# 1. Source contract checks for Issue #16 acceptance criteria in tool_weighting.qmd
def test_dashboard_source_evidence_and_warning_contracts():
    qmd_content = Path("tool_weighting.qmd").read_text()

    # Near-tie label must be 'similar overall scores' and avoid statistical equivalence claim
    assert "similar overall scores" in qmd_content
    assert "statistical" in qmd_content.lower() or "statistically" in qmd_content.lower()

    # Recommendation selection control must exist in Region 4
    assert "selected_combo" in qmd_content

    # Labelled per-isolate dot plots for accuracy and contiguity with cohort means
    assert "dot-plot" in qmd_content.lower() or "dot_plot" in qmd_content.lower() or "svg" in qmd_content.lower()
    assert "Cohort Mean" in qmd_content or "cohort mean" in qmd_content.lower()

    # Compact per-isolate event table for residual hits and replicon losses
    assert "Residual Adapter/Barcode Hits" in qmd_content or "residual" in qmd_content.lower()
    assert "Replicon" in qmd_content

    # Supporting diagnostics: duplication ratio and misassemblies without primary sliders
    assert "Duplication_ratio" in qmd_content or "duplication" in qmd_content.lower()
    assert "misassemblies" in qmd_content.lower()
    # Ensure no primary slider/input was added for duplication or misassembly
    assert "ui.input_numeric(\"w_dup\"" not in qmd_content
    assert "ui.input_numeric(\"w_mis\"" not in qmd_content

    # Mismatches and indels individually visible without primary sliders
    assert "mean_mismatches" in qmd_content
    assert "mean_indels" in qmd_content
    assert "ui.input_numeric(\"w_mismatch\"" not in qmd_content
    assert "ui.input_numeric(\"w_indel\"" not in qmd_content

    # No population confidence interval or radar chart is introduced
    assert "polar" not in qmd_content.lower()
    assert "errorbar" not in qmd_content.lower()
    assert "confidence_interval" not in qmd_content.lower()

    # Equivalent textual interpretation for plotted information
    assert "aria-label" in qmd_content or "text-equivalent" in qmd_content or "table" in qmd_content


# 2. Parity check for calculate_scenario_rankings warnings and near-ties
def test_dashboard_calculate_scenario_rankings_warnings_and_near_ties():
    clean_code = _extract_app_code()
    scope = {}
    exec(clean_code, scope)

    calc_func = scope["calculate_scenario_rankings"]
    df = pd.read_csv("assets/data/assembly_metrics.csv")

    eligible, excluded, insufficient = calc_func(
        df=df,
        model="hac",
        depth="100x",
        weight_accuracy=28.0,
        weight_contiguity=20.0,
        weight_residual=17.0,
        weight_replicon=35.0,
        gate_complete=False,
        gate_zero_hits=False,
    )

    assert len(eligible) == 17
    # Verify near-tie calculation matches canonical scorer
    canonical_res = score_benchmark(
        df,
        Scenario(model="hac", depth="100x"),
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
    )
    for el, can in zip(eligible, canonical_res.recommendations[:len(eligible)]):
        assert el["combo"] == can.combo
        assert el["is_near_tie"] == can.is_near_tie
        # Any residual hit warning must report total hits and affected isolates
        if el["t_hits"] > 0:
            warn = [w for w in el["warnings"] if "residual adapter/barcode hit(s)" in w]
            assert len(warn) == 1
            assert f"{el['t_hits']} residual adapter/barcode hit(s)" in warn[0]
            assert f"{el['affected_residual_isolates']} isolate(s)" in warn[0]
        # Any replicon warning must report total, full, partial, and affected isolates
        if el["t_miss"] > 0:
            warn = [w for w in el["warnings"] if "missed or severely incomplete replicon(s)" in w]
            assert len(warn) == 1
            assert f"{el['t_miss']} missed or severely incomplete replicon(s)" in warn[0]
            assert f"{el['affected_replicon_isolates']} isolate(s)" in warn[0]
            assert f"{el['full_missed']} full" in warn[0]
            assert f"{el['partial_missed']} partial" in warn[0]


# 3. Playwright browser journey testing Issue #16 acceptance criteria
def test_browser_journey_warnings_and_per_isolate_evidence(tmp_path):
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

            # Scenario selection begins unset
            page.select_option("#model_select", "hac")
            page.select_option("#depth_select", "100x")

            # 1. Verify 'similar overall scores' badge is displayed on near-tie combinations
            page.wait_for_selector(".ranking-table")
            near_tie_badges = page.locator(".near-tie-badge")
            assert near_tie_badges.count() > 0
            first_badge_text = near_tie_badges.first.text_content()
            assert "similar overall scores" in first_badge_text

            # 2. Check shortlist table headers and contents: mismatches and indels visible
            table_text = page.locator(".ranking-table").text_content()
            assert "mm" in table_text
            assert "indel" in table_text

            # 3. Check selected combination inspector in Region 4
            # Combo selector exists and has choices
            page.wait_for_selector("#selected_combo")
            selected_val = page.input_value("#selected_combo")
            assert selected_val != ""

            # Check that dot plots for accuracy and contiguity are rendered
            page.wait_for_selector("#accuracy-dot-plot")
            page.wait_for_selector("#contiguity-dot-plot")

            # Check that cohort mean lines and text are present in the dot plots
            acc_plot_svg = page.locator("#accuracy-dot-plot").inner_html()
            assert "Cohort Mean" in acc_plot_svg or "mean" in acc_plot_svg.lower()

            cont_plot_svg = page.locator("#contiguity-dot-plot").inner_html()
            assert "Cohort Mean" in cont_plot_svg or "mean" in cont_plot_svg.lower()

            # Check that dot plots have equivalent textual interpretation
            page.wait_for_selector("#text-interpretation-accuracy")
            page.wait_for_selector("#text-interpretation-contiguity")

            # 4. Check compact per-isolate sparse events table (residual hits and replicon losses)
            page.wait_for_selector("#sparse-events-table")
            sparse_table = page.locator("#sparse-events-table").text_content()
            assert "Residual Adapter/Barcode Hits" in sparse_table
            assert "Replicon Losses" in sparse_table

            # 5. Check supporting duplication ratio and misassemblies in detail
            detail_text = page.locator("#region-detail").text_content()
            assert "Duplication Ratio" in detail_text or "Duplication ratio" in detail_text
            assert "Misassemblies" in detail_text or "misassemblies" in detail_text

            # 6. Select a different recommendation from the dropdown (e.g. index 1)
            options = page.locator("#selected_combo option").all_inner_texts()
            assert len(options) >= 5
            page.select_option("#selected_combo", index=1)
            # Check that Region 4 updates to show the newly selected combination
            new_selected_val = page.input_value("#selected_combo")
            assert new_selected_val != selected_val
            page.wait_for_selector(f"#region-detail:has-text('{new_selected_val}')")

            # 7. Verify zero-weight warnings remain visible
            # Reveal all combinations to inspect the complete 17-combination cohort
            page.click("#btn_reveal_all")
            page.wait_for_selector(".ranking-table tbody tr:nth-child(17)")

            # Set residual weight to 0.0 and increase accuracy to 45.0
            page.fill("#w_res", "0")
            page.fill("#w_acc", "45")
            page.wait_for_selector("#weight_total_indicator:has-text('100%')")
            page.wait_for_timeout(500)

            # Combinations with residual hits (e.g. chopper-untrimmed, seqkit-untrimmed) should STILL have warnings visible in Diagnostics column
            diag_column = page.locator(".ranking-table tbody tr td:last-child").all_inner_texts()
            has_res_warning = any("residual adapter/barcode hit(s)" in cell for cell in diag_column)
            assert has_res_warning is True, "Residual hit warnings must remain visible when residual weight is 0"

            browser.close()
    finally:
        proc.terminate()
        proc.wait()
