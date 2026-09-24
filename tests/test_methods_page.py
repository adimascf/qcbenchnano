import re
import subprocess
from pathlib import Path
import pytest
import yaml


def test_methods_source_exists_and_in_navbar():
    # 1. Source file exists
    methods_qmd = Path("methods_and_robustness.qmd")
    assert methods_qmd.exists(), "methods_and_robustness.qmd must exist in the repo root"

    # 2. Included in _quarto.yml navbar
    quarto_yml = Path("_quarto.yml")
    assert quarto_yml.exists()
    with open(quarto_yml, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)

    navbar_left = config.get("website", {}).get("navbar", {}).get("left", [])
    navbar_hrefs = [item.get("href") for item in navbar_left if isinstance(item, dict)]
    assert "methods_and_robustness.qmd" in navbar_hrefs, (
        "methods_and_robustness.qmd must be in the _quarto.yml navbar"
    )


def test_dashboard_links_to_methods_page():
    # Dashboard must link to methods_and_robustness.html without interrupting decision journey
    dash_source = Path("tool_weighting.qmd").read_text()
    assert "methods_and_robustness.html" in dash_source, (
        "tool_weighting.qmd must clearly link to methods_and_robustness.html"
    )


def test_methods_page_required_sections_and_disclosures():
    content = Path("methods_and_robustness.qmd").read_text()

    # Benchmark scope
    assert "13 reference bacterial isolates" in content or "13" in content
    assert "17" in content  # 17 combinations
    assert "884 observations" in content or "884" in content
    assert "50x" in content  # discussion of unsupported 50x depth
    assert "aggregate" in content or "cross-depth" in content  # rejection of cross-depth / aggregate scenarios

    # Four value functions
    # 1. Contiguity and deliberate symmetric auNGA
    assert "auNGA" in content
    assert "symmetric" in content
    assert "1.0" in content

    # 2. Sequence accuracy and project-policy calibration vs universal biological thresholds
    assert "Mismatches" in content or "mismatch" in content
    assert "Indels" in content or "indel" in content
    assert "10" in content  # 10 events / 100kbp
    assert "policy" in content or "calibration" in content

    # 3. Residual hits, 90% identity & 90% region-coverage, and existing validation
    assert "90%" in content
    assert "coverage" in content
    assert "identity" in content
    assert "validation" in content

    # 4. Replicon recovery (<50% severe loss, 95% complete gate, limits of coverage)
    assert "50%" in content
    assert "severe" in content or "loss" in content
    assert "95%" in content
    assert "circularisation" in content or "circular" in content or "structural correctness" in content

    # Cohort and arithmetic aggregation
    assert "arithmetic" in content
    assert "min–max" in content or "min-max" in content  # obsolete min-max excluded
    assert "geometric" in content  # obsolete geometric excluded

    # Presets and gates
    assert "Community-balanced" in content or "community_balanced" in content
    assert "Complete-replicon-recovery" in content or "complete_replicon_recovery" in content
    assert "Sequence-accurate-assembly" in content or "sequence_accurate_assembly" in content
    assert "Custom" in content or "custom" in content

    # Warnings, near-tie, and missing-data policy
    assert "25" in content  # 25-point mean-to-minimum gap
    assert "near-tie" in content or "similar overall scores" in content
    assert "insufficient" in content or "missing" in content

    # Scoring version and provenance
    assert "1.0" in content  # scoring version
    assert "4d6b8cb" in content  # pinned commit
    assert "ee7a50a" in content  # benchmark data hash
    assert "3e7cfdd" in content  # survey data hash

    # Robustness outputs linked (repository only, not in paper/supplement)
    assert "denominator_sensitivity.csv" in content or "denominator" in content
    assert "leave_one_out.csv" in content or "leave-one-isolate-out" in content
    assert "respondent_profiles.csv" in content or "respondent" in content
    assert "paper" in content or "supplement" in content

    # Link to Specification #10
    assert "issues/10" in content or "#10" in content


def test_methods_quarto_render_smoke():
    cmd = ["quarto", "render", "methods_and_robustness.qmd", "--to", "html", "--no-execute"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Quarto render failed: {res.stderr}\n{res.stdout}"

    output_html = Path("methods_and_robustness.html")
    if not output_html.exists():
        output_html = Path("_site/methods_and_robustness.html")
    assert output_html.exists(), "Rendered output methods_and_robustness.html should exist"

    html_content = output_html.read_text()
    assert "Methods & Robustness" in html_content or "Methods &amp; Robustness" in html_content
    assert "4d6b8cb" in html_content
