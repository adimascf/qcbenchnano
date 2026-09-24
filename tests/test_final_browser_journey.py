import io
import re
import socket
import subprocess
import time
from pathlib import Path
import pytest
import pandas as pd
from playwright.sync_api import sync_playwright


def _extract_app_code(qmd_path: Path = Path("tool_weighting.qmd")) -> str:
    content = qmd_path.read_text()
    match = re.search(r"```\{shinylive-python\}(.*?)```", content, re.DOTALL)
    assert match is not None, "shinylive-python block not found in tool_weighting.qmd"
    code = match.group(1)
    lines = [l for l in code.splitlines() if not l.strip().startswith("#|")]
    return "\n".join(lines)


def test_final_browser_comprehensive_journey_and_accessibility(tmp_path):
    """
    Acceptance Criteria 10 & 11:
    - A final browser journey covers scenario selection, named and custom preferences,
      gates, no-eligible state, warnings, details, reveal-all, and ordered download.
    - Keyboard operation, explicit labels, text equivalents, and non-colour-only
      status presentation pass the project's chosen accessibility checks.
    """
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

            # -----------------------------------------------------------------
            # 1. ACCESSIBILITY & INITIAL STATE INSPECTION
            # -----------------------------------------------------------------
            # Verify 5 accessible landmark flow regions with role="region" and aria-labelledby
            for region_id in [
                "region-scenario",
                "region-priorities",
                "region-shortlist",
                "region-detail",
                "region-methodology",
            ]:
                loc = page.locator(f"#{region_id}")
                assert loc.count() == 1, f"Region {region_id} must exist"
                assert loc.get_attribute("role") == "region"
                assert loc.get_attribute("aria-labelledby") is not None

            # Verify explicit labels for primary interactive controls
            assert page.locator("label[for='model_select']").count() >= 1 or page.locator("#model_select").get_attribute("aria-label")
            assert page.locator("label[for='depth_select']").count() >= 1 or page.locator("#depth_select").get_attribute("aria-label")
            assert page.locator("label[for='preset_select']").count() >= 1 or page.locator("#preset_select").get_attribute("aria-label")

            # Check Methods & Robustness link is present without interrupting journey
            methods_link = page.locator("a[href='methods_and_robustness.html']")
            assert methods_link.count() >= 1, "Must contain link to methods_and_robustness.html"

            # Initial state: scenario is unset, prompting user, no download button
            assert page.locator("#download_recommendations_csv").count() == 0
            page.wait_for_selector("#shortlist_content .alert-info")
            shortlist_text = page.locator("#shortlist_content").text_content()
            assert "select both a basecalling model and sequencing depth" in shortlist_text

            # -----------------------------------------------------------------
            # 2. KEYBOARD OPERATION & SCENARIO SELECTION
            # -----------------------------------------------------------------
            # Keyboard Tab navigation and focus verification
            page.focus("#model_select")
            assert page.evaluate("document.activeElement.id") == "model_select"
            page.select_option("#model_select", "hac")

            page.keyboard.press("Tab")
            assert page.evaluate("document.activeElement.id") == "depth_select"
            page.select_option("#depth_select", "100x")

            # Recommendations load: leading combo is seqkit-barbell
            page.wait_for_selector(".ranking-table")
            first_row = page.locator(".ranking-table tbody tr").first
            assert "seqkit-barbell" in first_row.text_content()
            assert "98.7" in first_row.text_content()

            # CSV Download button is now rendered
            page.wait_for_selector("#download_recommendations_csv")

            # -----------------------------------------------------------------
            # 3. NAMED PREFERENCES & SLIDERS
            # -----------------------------------------------------------------
            # Preset 1: Community-balanced (default)
            assert page.input_value("#preset_select") == "community_balanced"
            assert page.input_value("#w_acc") in ("28", "28.0")
            assert page.input_value("#w_cont") in ("20", "20.0")
            assert page.input_value("#w_res") in ("17", "17.0")
            assert page.input_value("#w_rep") in ("35", "35.0")
            assert not page.is_checked("#gate_complete")
            assert not page.is_checked("#gate_zero_hits")

            # Preset 2: Sequence-accurate-assembly
            page.select_option("#preset_select", "sequence_accurate_assembly")
            page.wait_for_selector("#region-shortlist table tbody tr:first-child:has-text('97.6')")
            assert page.input_value("#w_acc") in ("50", "50.0")
            assert page.input_value("#w_cont") in ("14", "14.0")

            # -----------------------------------------------------------------
            # 4. CUSTOM PREFERENCES, VALIDATION & RESET
            # -----------------------------------------------------------------
            # Modify a weight -> switches visibly to custom and pauses ranking when sum != 100
            page.fill("#w_acc", "45")
            page.wait_for_function("document.querySelector('#preset_select').value === 'custom'")
            assert page.input_value("#preset_select") == "custom"

            # Sum != 100% pauses ranking and shows remaining/excess
            page.wait_for_selector("#region-shortlist .alert-warning")
            assert "Ranking paused" in page.text_content("#region-shortlist")

            # Restore total to 100 with decimals in custom mode: 45 + 19 + 12 + 24 = 100
            page.fill("#w_cont", "19")
            page.wait_for_function("document.querySelector('#region-priorities').textContent.includes('100%')")
            page.wait_for_selector(".ranking-table tbody tr:first-child")

            # Reset back to Community-balanced restores 28, 20, 17, 35
            page.click("#btn_reset")
            page.wait_for_function("document.querySelector('#preset_select').value === 'community_balanced'")
            page.wait_for_selector("#region-shortlist table tbody tr:first-child:has-text('98.7')")
            assert page.input_value("#preset_select") == "community_balanced"
            assert page.input_value("#w_acc") in ("28", "28.0")

            # Preset 3: Complete-replicon-recovery
            page.select_option("#preset_select", "complete_replicon_recovery")
            page.wait_for_selector("#region-shortlist details summary:has-text('Excluded')")
            assert page.is_checked("#gate_complete")

            # -----------------------------------------------------------------
            # 5. GATES & EXCLUSION SEMANTICS
            # -----------------------------------------------------------------
            # Reset to Community-balanced
            page.click("#btn_reset")
            page.wait_for_function("document.querySelector('#preset_select').value === 'community_balanced'")
            page.wait_for_function("document.querySelector('#w_acc').value == '28' || document.querySelector('#w_acc').value == '28.0'")

            # Enable Complete-replicon-recovery gate
            page.check("#gate_complete")
            page.wait_for_selector("#region-shortlist details summary:has-text('(10)')")
            excluded_summary = page.text_content("#region-shortlist details summary")
            assert "Excluded / Ineligible Combinations (10)" in excluded_summary

            # Non-colour-only exclusion status: explicit text reasons in excluded list
            excluded_details = page.text_content("#region-shortlist details")
            assert "fails complete-recovery gate" in excluded_details

            # -----------------------------------------------------------------
            # 6. ALL-INELIGIBLE / NO-ELIGIBLE STATE
            # -----------------------------------------------------------------
            # Enable both gates; in HAC 100x 5 pass.
            # To test the all-ineligible state, check how UI behaves when no combo passes.
            # (Unit and component tests verify zero-pass state; here we verify gate state preservation)
            page.check("#gate_zero_hits")
            page.wait_for_selector("#region-shortlist details summary:has-text('(12)')")
            # Requirements are preserved without relaxation
            assert page.is_checked("#gate_complete")
            assert page.is_checked("#gate_zero_hits")

            # Reset to restore Community-balanced and uncheck all gates
            page.click("#btn_reset")
            page.wait_for_selector("#btn_reveal_all:has-text('Reveal All 17 Combinations')")

            # -----------------------------------------------------------------
            # 7. REVEAL-ALL INTERACTION & FULL COHORT VIEW
            # -----------------------------------------------------------------
            # Initially top 5 displayed
            assert page.locator(".ranking-table tbody tr").count() == 5
            btn_reveal = page.locator("#btn_reveal_all")
            assert "Reveal All 17 Combinations" in btn_reveal.text_content()

            # Click to reveal all 17
            btn_reveal.click()
            page.wait_for_function("document.querySelectorAll('.ranking-table tbody tr').length === 17")
            assert page.locator(".ranking-table tbody tr").count() == 17
            assert "Show Top 5 Only" in btn_reveal.text_content()

            # -----------------------------------------------------------------
            # 8. WARNINGS & NON-COLOUR-ONLY PRESENTATION
            # -----------------------------------------------------------------
            # Verify non-colour-only presentation:
            # - Near-tie indicator uses explicit text: "similar overall scores"
            near_tie_elements = page.locator(".near-tie-badge")
            assert near_tie_elements.count() > 0
            assert "similar overall scores" in near_tie_elements.first.text_content()

            # - Clean indicator has checkmark + text: "✓ Clean"
            assert "✓ Clean" in page.locator(".ranking-table").text_content()

            # - Warnings have explicit text descriptions across the full 17 combinations
            diag_column_text = page.locator(".ranking-table tbody tr td:last-child").all_inner_texts()
            assert any("hit(s)" in t or "replicon(s)" in t or "Variable" in t for t in diag_column_text)

            # -----------------------------------------------------------------
            # 9. DETAILS INSPECTION (DOT PLOTS, SPARSE TABLE, TEXT EQUIVALENTS)
            # -----------------------------------------------------------------
            # Dropdown combination inspector in Region 4
            page.wait_for_selector("#selected_combo")
            page.select_option("#selected_combo", value="seqkit-barbell")
            page.wait_for_selector("#accuracy-dot-plot")
            page.wait_for_selector("#contiguity-dot-plot")

            # Dot plots have role="img", aria-label, and cohort mean indicators
            acc_svg = page.locator("#accuracy-dot-plot")
            assert acc_svg.get_attribute("role") == "img"
            assert "Sequence Accuracy" in (acc_svg.get_attribute("aria-label") or "")
            assert "Cohort Mean" in acc_svg.inner_html()

            # Text equivalents for visual charts
            page.wait_for_selector("#text-interpretation-accuracy")
            page.wait_for_selector("#text-interpretation-contiguity")

            # Sparse-event isolate table exists and has aria-label
            sparse_tbl = page.locator("#sparse-events-table")
            assert sparse_tbl.count() == 1
            assert sparse_tbl.get_attribute("aria-label") is not None

            # -----------------------------------------------------------------
            # 10. ORDERED REPRODUCIBLE CSV EXPORT
            # -----------------------------------------------------------------
            with page.expect_download() as download_info:
                page.click("#download_recommendations_csv")
            download = download_info.value
            csv_path = download.path()
            csv_content = Path(csv_path).read_text()
            df_exported = pd.read_csv(io.StringIO(csv_content), keep_default_na=False)

            # Verify 17 records preserving exact visible display order
            assert len(df_exported) == 17
            table_combos = page.locator(".ranking-table tbody tr td:nth-child(2)").all_inner_texts()
            for idx, combo_name in enumerate(table_combos):
                assert df_exported.iloc[idx]["combo"] == combo_name

            # Verify required columns and provenance
            assert "rank" in df_exported.columns
            assert "overall_score" in df_exported.columns
            assert "scoring_version" in df_exported.columns
            assert "source_data_commit" in df_exported.columns
            assert str(df_exported.iloc[0]["scoring_version"]) == "1.0"
            assert df_exported.iloc[0]["source_data_commit"] == "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"
            assert "raw.githubusercontent.com" not in csv_content

            browser.close()
    finally:
        proc.terminate()
        proc.wait()


def test_final_browser_journey_all_ineligible_state(tmp_path):
    """
    Acceptance Criterion 10:
    Direct browser verification of the all-ineligible / no-eligible state:
    - Verifies #no-eligible-msg alert is displayed
    - Verifies requirements are not relaxed
    - Verifies failed combinations are not promoted (table rows = 0)
    - Verifies details shows 'No combination meets active eligibility requirements'
    - Verifies excluded container shows all 17 failed combinations with reasons
    - Verifies ordered CSV download contains all 17 combinations with blank ranks
    """
    # Create app copy pointing to synthetic metrics file where all combos fail the zero-residual-hits gate
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    df_all_failing = df.copy()
    df_all_failing["contamination_count"] = 5
    failing_csv = tmp_path / "assembly_metrics.csv"
    df_all_failing.to_csv(failing_csv, index=False)

    qmd_path = Path("tool_weighting.qmd")
    app_code = _extract_app_code(qmd_path)
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

            # 1. Wait for all-ineligible warning alert
            page.wait_for_selector("#no-eligible-msg")
            msg = page.text_content("#no-eligible-msg")
            assert "No evaluated combination meets the active eligibility requirements" in msg
            assert "All gates have been preserved" in msg

            # 2. Shortlist table has 0 eligible rows and no failed combination is promoted
            assert page.locator("#region-shortlist table tbody tr").count() == 0

            # 3. Gates remain checked (not relaxed)
            assert page.is_checked("#gate_zero_hits")

            # 4. Detail section reports that no combination meets requirements
            detail_msg = page.text_content("#region-detail")
            assert "No combination meets active eligibility requirements" in detail_msg

            # 5. Excluded section shows all 17 failed combinations
            page.wait_for_selector("#region-shortlist details summary:has-text('(17)')")
            summary_text = page.text_content("#region-shortlist details summary")
            assert "(17)" in summary_text

            # 6. Download button remains functional in all-ineligible state
            with page.expect_download() as download_info:
                page.click("#download_recommendations_csv")
            download = download_info.value
            csv_path = download.path()
            csv_content = Path(csv_path).read_text()
            df_ineligible_csv = pd.read_csv(io.StringIO(csv_content), keep_default_na=False)

            assert len(df_ineligible_csv) == 17
            assert (df_ineligible_csv["is_eligible"] == False).all()
            assert (df_ineligible_csv["rank"] == "").all()

            browser.close()
    finally:
        proc.terminate()
        proc.wait()

