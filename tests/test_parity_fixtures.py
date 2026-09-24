import re
from pathlib import Path
import pytest
import pandas as pd

from qc_scoring.models import Scenario
from qc_scoring.preferences import get_preset, PresetName, WeightsConfig, GateConfig
from qc_scoring.scorer import score_benchmark
from qc_scoring.paper import generate_paper_composite_results, SUPPORTED_PAPER_SCENARIOS


def _extract_dashboard_calc_func():
    dash_source = Path("tool_weighting.qmd").read_text()
    match = re.search(r"```\{shinylive-python\}(.*?)```", dash_source, re.DOTALL)
    assert match is not None, "shinylive-python block not found in tool_weighting.qmd"
    code = match.group(1)
    lines = [l for l in code.splitlines() if not l.strip().startswith("#|")]
    scope = {}
    exec("\n".join(lines), scope)
    return scope["calculate_scenario_rankings"]


@pytest.fixture(scope="module")
def benchmark_df():
    return pd.read_csv("assets/data/assembly_metrics.csv")


@pytest.fixture(scope="module")
def dashboard_calc_func():
    return _extract_dashboard_calc_func()


# -----------------------------------------------------------------------------
# Canonical Fixtures Table: Known Canonical Leaders for All 4 Supported Scenarios
# under Community-balanced preferences
# -----------------------------------------------------------------------------
CANONICAL_COMMUNITY_LEADERS = {
    Scenario("hac", "20x"): {
        "leader": "seqkit-dorado",
        "rank": 1,
        "display_score": 94.6,
        "score_accuracy": 80.9,
        "score_contiguity": 99.9,
        "score_residual": 100.0,
        "score_replicon": 100.0,
    },
    Scenario("hac", "100x"): {
        "leader": "seqkit-barbell",
        "rank": 1,
        "display_score": 98.7,
        "score_accuracy": 95.3,
        "score_contiguity": 100.0,
        "score_residual": 100.0,
        "score_replicon": 100.0,
    },
    Scenario("sup", "20x"): {
        "leader": "chopper-porechop_abi",
        "rank": 1,
        "display_score": 97.5,
        "score_accuracy": 91.1,
        "score_contiguity": 100.0,
        "score_residual": 100.0,
        "score_replicon": 100.0,
    },
    Scenario("sup", "100x"): {
        "leader": "chopper-porechop_abi",
        "rank": 1,
        "display_score": 99.3,
        "score_accuracy": 97.5,
        "score_contiguity": 100.0,
        "score_residual": 100.0,
        "score_replicon": 100.0,
    },
}


@pytest.mark.parametrize("scenario", SUPPORTED_PAPER_SCENARIOS)
def test_paper_dashboard_canonical_criterion_score_parity(scenario, benchmark_df, dashboard_calc_func):
    """
    Acceptance Criterion 9:
    Canonical fixtures prove exact paper/dashboard criterion-score, overall-score,
    eligibility, warning, and rank parity for all supported scenarios.
    """
    preset = get_preset(PresetName.COMMUNITY_BALANCED)

    # 1. Canonical scorer
    canonical_result = score_benchmark(
        benchmark_df,
        scenario,
        weights=preset.weights,
        gates=preset.gates,
        preset="community_balanced",
    )

    # 2. Paper composite generator
    paper_results = generate_paper_composite_results(benchmark_df)
    paper_result = paper_results[scenario]

    # 3. Dashboard calculation engine
    dash_eligible, dash_excluded, *dash_rest = dashboard_calc_func(
        benchmark_df,
        model=scenario.model,
        depth=scenario.depth,
        weight_accuracy=preset.weights.accuracy,
        weight_contiguity=preset.weights.contiguity,
        weight_residual=preset.weights.residual,
        weight_replicon=preset.weights.replicon,
        gate_complete=preset.gates.complete_recovery,
        gate_zero_hits=preset.gates.zero_residual_hits,
    )

    # Verify Leader Matches Known Canonical Fixture
    expected_fixture = CANONICAL_COMMUNITY_LEADERS[scenario]
    canonical_leader = canonical_result.leading_recommendation
    assert canonical_leader is not None
    assert canonical_leader.combo == expected_fixture["leader"]
    assert canonical_leader.rank == expected_fixture["rank"]
    assert canonical_leader.display_score == expected_fixture["display_score"]
    assert round(canonical_leader.score_accuracy, 1) == expected_fixture["score_accuracy"]
    assert round(canonical_leader.score_contiguity, 1) == expected_fixture["score_contiguity"]
    assert round(canonical_leader.score_residual, 1) == expected_fixture["score_residual"]
    assert round(canonical_leader.score_replicon, 1) == expected_fixture["score_replicon"]

    # All 17 combinations are eligible under Community-balanced
    assert len(canonical_result.recommendations) == 17
    assert len(paper_result.recommendations) == 17
    assert len(dash_eligible) == 17
    assert len(dash_excluded) == 0

    # Cross-compare all 17 combinations across Canonical, Paper, and Dashboard
    for idx in range(17):
        can_rec = canonical_result.recommendations[idx]
        pap_rec = paper_result.recommendations[idx]
        dash_rec = dash_eligible[idx]

        # 1. Combo identity and rank parity
        assert can_rec.combo == pap_rec.combo == dash_rec["combo"]
        assert can_rec.rank == pap_rec.rank == dash_rec["rank"] == (idx + 1)

        # 2. Overall score parity
        assert abs(can_rec.overall_score - pap_rec.overall_score) < 1e-8
        assert abs(can_rec.overall_score - dash_rec["overall_score"]) < 1e-5
        assert can_rec.display_score == pap_rec.display_score == dash_rec["display_score"]

        # 3. Criterion score parity
        assert abs(can_rec.score_accuracy - pap_rec.score_accuracy) < 1e-8
        assert abs(can_rec.score_contiguity - pap_rec.score_contiguity) < 1e-8
        assert abs(can_rec.score_residual - pap_rec.score_residual) < 1e-8
        assert abs(can_rec.score_replicon - pap_rec.score_replicon) < 1e-8

        assert dash_rec["accuracy"] == pytest.approx(can_rec.score_accuracy, abs=0.1)
        assert dash_rec["contiguity"] == pytest.approx(can_rec.score_contiguity, abs=0.1)
        assert dash_rec["residual"] == pytest.approx(can_rec.score_residual, abs=0.1)
        assert dash_rec["replicon"] == pytest.approx(can_rec.score_replicon, abs=0.1)

        # 4. Eligibility parity
        assert can_rec.is_eligible is True
        assert pap_rec.is_eligible is True
        assert dash_rec["is_eligible"] is True

        # 5. Near-tie parity
        assert can_rec.is_near_tie == pap_rec.is_near_tie == dash_rec["is_near_tie"]

        # 6. Warnings parity
        can_warnings = sorted(can_rec.warnings)
        pap_warnings = sorted(pap_rec.warnings)
        dash_warnings = sorted(dash_rec["warnings"])
        assert can_warnings == pap_warnings == dash_warnings


@pytest.mark.parametrize("scenario", SUPPORTED_PAPER_SCENARIOS)
@pytest.mark.parametrize("preset_name", [
    PresetName.COMPLETE_REPLICON_RECOVERY,
    PresetName.SEQUENCE_ACCURATE_ASSEMBLY,
])
def test_dashboard_canonical_all_presets_and_gates_parity(scenario, preset_name, benchmark_df, dashboard_calc_func):
    """
    Test parity under alternative named presets with gates across all 4 scenarios.
    """
    preset = get_preset(preset_name)

    canonical_result = score_benchmark(
        benchmark_df,
        scenario,
        weights=preset.weights,
        gates=preset.gates,
        preset=preset_name.value,
    )

    dash_eligible, dash_excluded, *dash_rest = dashboard_calc_func(
        benchmark_df,
        model=scenario.model,
        depth=scenario.depth,
        weight_accuracy=preset.weights.accuracy,
        weight_contiguity=preset.weights.contiguity,
        weight_residual=preset.weights.residual,
        weight_replicon=preset.weights.replicon,
        gate_complete=preset.gates.complete_recovery,
        gate_zero_hits=preset.gates.zero_residual_hits,
    )

    can_eligible = [r for r in canonical_result.recommendations if r.is_eligible]
    can_excluded = [r for r in canonical_result.recommendations if not r.is_eligible]

    assert len(dash_eligible) == len(can_eligible)
    assert len(dash_excluded) == len(can_excluded)

    # Check eligible items
    for idx, (can_rec, dash_rec) in enumerate(zip(can_eligible, dash_eligible)):
        assert can_rec.combo == dash_rec["combo"]
        assert can_rec.rank == dash_rec["rank"]
        assert abs(can_rec.overall_score - dash_rec["overall_score"]) < 1e-5
        assert can_rec.display_score == dash_rec["display_score"]
        assert dash_rec["is_eligible"] is True
        assert sorted(can_rec.warnings) == sorted(dash_rec["warnings"])

    # Check excluded items
    for can_rec, dash_rec in zip(can_excluded, dash_excluded):
        assert can_rec.combo == dash_rec["combo"]
        assert can_rec.rank is None
        assert dash_rec["rank"] is None
        assert dash_rec["is_eligible"] is False
        assert can_rec.ineligible_reason is not None
        assert dash_rec["ineligible_reason"] is not None
