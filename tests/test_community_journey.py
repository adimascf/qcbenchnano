import re
import subprocess
import time
from pathlib import Path
import pytest
import pandas as pd
from qc_scoring.models import Scenario
from qc_scoring.preferences import get_preset, PresetName
from qc_scoring.scorer import score_benchmark
from qc_scoring.config import load_dashboard_config


# 1. Canonical Community-balanced fixtures test
def test_canonical_community_balanced_fixtures():
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    preset = get_preset(PresetName.COMMUNITY_BALANCED)

    expected_leaders = {
        ("hac", "20x"): {
            "combo": "seqkit-dorado",
            "display_score": 94.6,
            "score_accuracy": 80.9,
            "score_contiguity": 99.9,
            "score_residual": 100.0,
            "score_replicon": 100.0,
            "mean_error_rate": 1.91,
            "mean_auNGA_ratio": 0.999,
            "residual_total_hits": 0,
            "replicon_total_missed": 0,
        },
        ("hac", "100x"): {
            "combo": "seqkit-barbell",
            "display_score": 98.7,
            "score_accuracy": 95.3,
            "score_contiguity": 100.0,
            "score_residual": 100.0,
            "score_replicon": 100.0,
            "mean_error_rate": 0.47,
            "mean_auNGA_ratio": 1.000,
            "residual_total_hits": 0,
            "replicon_total_missed": 0,
        },
        ("sup", "20x"): {
            "combo": "chopper-porechop_abi",
            "display_score": 97.5,
            "score_accuracy": 91.1,
            "score_contiguity": 100.0,
            "score_residual": 100.0,
            "score_replicon": 100.0,
            "mean_error_rate": 0.89,
            "mean_auNGA_ratio": 1.000,
            "residual_total_hits": 0,
            "replicon_total_missed": 0,
        },
        ("sup", "100x"): {
            "combo": "chopper-porechop_abi",
            "display_score": 99.3,
            "score_accuracy": 97.5,
            "score_contiguity": 100.0,
            "score_residual": 100.0,
            "score_replicon": 100.0,
            "mean_error_rate": 0.25,
            "mean_auNGA_ratio": 1.000,
            "residual_total_hits": 0,
            "replicon_total_missed": 0,
        },
    }

    for (model, depth), exp in expected_leaders.items():
        sc = Scenario(model=model, depth=depth)
        res = score_benchmark(df, sc, weights=preset.weights, gates=preset.gates)
        assert len(res.recommendations) == 17
        assert res.eligible_count == 17
        leading = res.leading_recommendation
        assert leading is not None
        assert leading.rank == 1
        assert leading.combo == exp["combo"]
        assert leading.display_score == exp["display_score"]
        assert round(leading.score_accuracy, 1) == exp["score_accuracy"]
        assert round(leading.score_contiguity, 1) == exp["score_contiguity"]
        assert round(leading.score_residual, 1) == exp["score_residual"]
        assert round(leading.score_replicon, 1) == exp["score_replicon"]
        assert round(leading.mean_error_rate, 2) == exp["mean_error_rate"]
        assert round(leading.mean_auNGA_ratio, 3) == exp["mean_auNGA_ratio"]
        assert leading.residual_total_hits == exp["residual_total_hits"]
        assert leading.replicon_total_missed == exp["replicon_total_missed"]


# 2. Static source contract checks on tool_weighting.qmd
def test_dashboard_source_journey_specifications():
    qmd_content = Path("tool_weighting.qmd").read_text()

    # Unset scenario selection
    assert 'selected=""' in qmd_content, "Scenario selections must begin unset"
    assert '"hac":' in qmd_content and '"sup":' in qmd_content
    assert '"20x":' in qmd_content and '"100x":' in qmd_content

    # Absence of unsupported options
    assert "50x" not in qmd_content, "50x depth must be absent"
    assert "all-model" not in qmd_content and "all_models" not in qmd_content
    assert "all-depth" not in qmd_content and "all_depths" not in qmd_content

    # Community-balanced percentages: 28/20/17/35
    assert "28" in qmd_content and "20" in qmd_content and "17" in qmd_content and "35" in qmd_content

    # No excellent/good/poor grades
    for forbidden in ["excellent", "good", "poor", "grade"]:
        matches = re.findall(rf"\b{forbidden}\b", qmd_content, re.IGNORECASE)
        assert len(matches) == 0, f"Found forbidden quality grade term: '{forbidden}'"

    # Limitation to 13-isolate benchmark displayed
    assert "13" in qmd_content and "isolate" in qmd_content


def _extract_app_code(qmd_path: Path = Path("tool_weighting.qmd")) -> str:
    qmd_content = qmd_path.read_text()
    match = re.search(r"```\{shinylive-python\}\n(.*?)\n```", qmd_content, re.DOTALL)
    assert match is not None, "shinylive-python code chunk not found in dashboard source"
    app_lines = [l for l in match.group(1).split("\n") if not l.strip().startswith("#|")]
    return "\n".join(app_lines)


# 3. Parity between dashboard evaluated_results and canonical scorer
def test_dashboard_evaluated_results_parity():
    clean_code = _extract_app_code()
    scope = {}
    exec(clean_code, scope)

    calc_func = scope["calculate_scenario_rankings"]
    df = pd.read_csv("assets/data/assembly_metrics.csv")
    preset = get_preset(PresetName.COMMUNITY_BALANCED)

    for model in ["hac", "sup"]:
        for depth in ["20x", "100x"]:
            sc = Scenario(model=model, depth=depth)
            canonical_res = score_benchmark(df, sc, weights=preset.weights, gates=preset.gates)
            dash_eligible, dash_excluded, *dash_insufficient_rest = calc_func(
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
            assert len(dash_eligible) == 17
            assert len(dash_excluded) == 0

            for idx in range(17):
                can_rec = canonical_res.recommendations[idx]
                dash_rec = dash_eligible[idx]
                assert dash_rec["combo"] == can_rec.combo
                assert dash_rec["rank"] == can_rec.rank
                assert abs(dash_rec["overall_score"] - can_rec.overall_score) < 1e-6
                assert dash_rec["display_score"] == can_rec.display_score
                assert dash_rec["accuracy"] == round(can_rec.score_accuracy, 1)
                assert dash_rec["contiguity"] == round(can_rec.score_contiguity, 1)
                assert dash_rec["residual"] == round(can_rec.score_residual, 1)
                assert dash_rec["replicon"] == round(can_rec.score_replicon, 1)
                assert abs(dash_rec["mean_error_rate"] - can_rec.mean_error_rate) < 1e-4
                assert abs(dash_rec["mean_aunga"] - can_rec.mean_auNGA_ratio) < 1e-4
                assert dash_rec["t_hits"] == can_rec.residual_total_hits
                assert dash_rec["t_miss"] == can_rec.replicon_total_missed


# 4. Browser-level acceptance check with Playwright (synchronous)
def test_browser_acceptance_community_journey(tmp_path):
    from playwright.sync_api import sync_playwright
    import socket

    config = load_dashboard_config()
    repo_owner = config.repository_owner
    repo_name = config.repository_name
    assert repo_owner is not None and len(repo_owner) > 0
    assert repo_name is not None and len(repo_name) > 0

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

            # Verify configured provenance in page text (repository independence check)
            methodology_text = page.text_content("#region-methodology")
            assert f"{repo_owner}/{repo_name}" in methodology_text, "Configured repository identity must appear in methodology provenance"
            assert config.pinned_commit[:7] in methodology_text, "Pinned commit must appear in methodology provenance"

            # 1. Verify model and depth begin unset
            model_val = page.input_value("#model_select")
            depth_val = page.input_value("#depth_select")
            assert model_val == "", "Model select must start unset"
            assert depth_val == "", "Depth select must start unset"

            # 2. Verify recommendations are not displayed when scenario is unset
            page.wait_for_selector("#region-shortlist .alert")
            initial_shortlist = page.text_content("#region-shortlist")
            assert "awaiting scenario selection" in initial_shortlist.lower() or "select both" in initial_shortlist.lower()
            rows = page.locator("#region-shortlist table tbody tr").count()
            assert rows == 0, "No recommendation rows should appear before scenario selection"

            # 3. Verify unsupported 50x, all-model, all-depth choices are absent
            model_options = page.locator("#model_select option").all_text_contents()
            depth_options = page.locator("#depth_select option").all_text_contents()
            assert not any("50x" in opt.lower() for opt in depth_options), "50x must be absent"
            assert not any(opt.strip().lower() in ("all", "all models", "all-models", "all model") for opt in model_options), "All-model must be absent"
            assert not any(opt.strip().lower() in ("all", "all depths", "all-depths", "all depth") for opt in depth_options), "All-depth must be absent"

            canonical_leaders = {
                ("hac", "20x"): ("seqkit-dorado", "94.6"),
                ("hac", "100x"): ("seqkit-barbell", "98.7"),
                ("sup", "20x"): ("chopper-porechop_abi", "97.5"),
                ("sup", "100x"): ("chopper-porechop_abi", "99.3"),
            }

            for (model, depth), (expected_combo, expected_score) in canonical_leaders.items():
                page.select_option("#model_select", model)
                page.select_option("#depth_select", depth)

                # Wait specifically for the first row to update to the expected scenario leader and score
                page.wait_for_selector(f"#region-shortlist table tbody tr:first-child:has-text('{expected_combo}')")
                page.wait_for_selector(f"#region-shortlist table tbody tr:first-child:has-text('{expected_score}')")

                first_row_text = page.locator("#region-shortlist table tbody tr").first.text_content()
                assert expected_combo in first_row_text
                assert expected_score in first_row_text

                # Verify all four raw summaries exist in the leading row
                assert "err/100kbp" in first_row_text
                assert "auNGA:" in first_row_text
                assert "hit(s)" in first_row_text
                assert "missed" in first_row_text

                initial_count = page.locator("#region-shortlist table tbody tr").count()
                assert initial_count == 5, f"Expected 5 initial rows, got {initial_count}"

                page.click("#btn_reveal_all")
                page.wait_for_function("document.querySelectorAll('#region-shortlist table tbody tr').length === 17")
                all_count = page.locator("#region-shortlist table tbody tr").count()
                assert all_count == 17, f"Expected 17 rows after reveal, got {all_count}"

                page.click("#btn_reveal_all")
                page.wait_for_function("document.querySelectorAll('#region-shortlist table tbody tr').length === 5")

            body_text = page.text_content("body")
            for word in ["excellent", "good", "poor"]:
                assert word not in body_text.lower()

            assert "13" in body_text and "isolate" in body_text.lower()

            browser.close()
    finally:
        proc.terminate()
        proc.wait()
