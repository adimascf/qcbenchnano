import os
import subprocess
from pathlib import Path
import pytest
from qc_scoring.config import DashboardConfig, load_dashboard_config

def test_load_default_dashboard_config():
    cfg = load_dashboard_config()
    assert cfg.repository_owner == "adimascf"
    assert cfg.repository_name == "biox7021"
    assert "adimascf.github.io" in cfg.site_base_url
    assert cfg.pinned_commit == "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"
    assert cfg.source_data_path is not None
    assert cfg.scoring_version == "1.0"

def test_repository_rename_configuration_independence(tmp_path):
    # Verify that a repository rename can be represented by changing configuration
    # without touching scoring logic or equations
    custom_yaml = tmp_path / "custom_config.yaml"
    custom_yaml.write_text("""
repository:
  owner: "new-org"
  name: "new-biox-repo"
  base_url: "https://new-org.github.io/new-biox-repo/"
  default_branch: "main"
evidence:
  source_data_path: "assets/data/assembly_metrics.csv"
  pinned_commit: "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562"
  bundled_data_path: "assets/data/assembly_metrics.csv"
scoring:
  version: "1.0"
""")
    cfg = load_dashboard_config(custom_yaml)
    assert cfg.repository_owner == "new-org"
    assert cfg.repository_name == "new-biox-repo"
    assert cfg.site_base_url == "https://new-org.github.io/new-biox-repo/"
    assert cfg.provenance_text() == "Source: new-org/new-biox-repo @ 4d6b8cb"

def test_dashboard_source_exists_and_has_agreed_regions():
    # Verify human-editable source exists and contains all 5 flow regions
    dash_source = Path("tool_weighting.qmd")
    assert dash_source.exists(), "Dashboard source tool_weighting.qmd must exist"
    
    content = dash_source.read_text()
    
    # 5 flow regions required by spec:
    # scenario, priorities, shortlist, selected detail, concise methodology
    assert "region-scenario" in content or "section-scenario" in content or "Scenario Selection" in content
    assert "region-priorities" in content or "section-priorities" in content or "Priorities" in content
    assert "region-shortlist" in content or "section-shortlist" in content or "Recommendation Shortlist" in content
    assert "region-detail" in content or "section-detail" in content or "Combination Details" in content
    assert "region-methodology" in content or "section-methodology" in content or "Methodology" in content

    # Pinned evidence / no moving branch head url
    assert "raw.githubusercontent.com/adimascf/biox7021/refs/heads/main" not in content, \
        "Dashboard must not silently fetch from moving branch head"

def test_dashboard_quarto_render_smoke():
    # Smoke check: quarto render on tool_weighting.qmd succeeds and contains configured provenance
    cmd = ["quarto", "render", "tool_weighting.qmd", "--to", "html", "--no-execute"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"Quarto render failed: {res.stderr}"
    
    output_html = Path("tool_weighting.html")
    if not output_html.exists():
        output_html = Path("_site/tool_weighting.html")
    assert output_html.exists(), "Rendered output html should exist"
    
    html_content = output_html.read_text()
    cfg = load_dashboard_config()
    # Check that configured provenance (e.g. pinned commit or repo) is visible
    assert cfg.pinned_commit[:7] in html_content or cfg.pinned_commit in html_content


def test_dashboard_pyodide_standalone_import_and_config(monkeypatch):
    # Regression test for Pyodide / Shinylive browser environment where qc_scoring is not installed
    import re
    import sys

    dash_source = Path("tool_weighting.qmd")
    content = dash_source.read_text()
    match = re.search(r"```\{shinylive-python\}\n(.*?)\n```", content, re.DOTALL)
    assert match is not None, "shinylive-python code chunk not found"
    app_lines = [l for l in match.group(1).split("\n") if not l.strip().startswith("#|")]
    app_code = "\n".join(app_lines)

    # Block qc_scoring imports in sys.modules to simulate Pyodide WASM environment
    monkeypatch.setitem(sys.modules, "qc_scoring", None)
    monkeypatch.setitem(sys.modules, "qc_scoring.config", None)

    scope = {}
    exec(app_code, scope)

    assert "app_config" in scope, "app_config must be loaded even in standalone Pyodide runtime"
    config = scope["app_config"]
    assert config.repository_owner == "adimascf"
    assert config.repository_name == "biox7021"
    assert config.scoring_version == "1.0"
    assert "https://raw.githubusercontent.com/adimascf/biox7021/" in config.data_url()
    assert "app" in scope, "Shiny App object must be instantiated without errors"

