from dataclasses import dataclass
from enum import Enum
from typing import Dict, Union
from qc_scoring.validation import ScoringError


class InvalidWeightsError(ScoringError):
    """Raised when weight allocation is invalid."""
    pass


@dataclass(frozen=True)
class WeightsConfig:
    accuracy: float
    contiguity: float
    residual: float
    replicon: float

    def as_dict(self) -> Dict[str, float]:
        return {
            "accuracy": self.accuracy,
            "contiguity": self.contiguity,
            "residual": self.residual,
            "replicon": self.replicon,
        }


@dataclass(frozen=True)
class GateConfig:
    complete_recovery: bool = False
    zero_residual_hits: bool = False

    def as_dict(self) -> Dict[str, bool]:
        return {
            "complete_recovery": self.complete_recovery,
            "zero_residual_hits": self.zero_residual_hits,
        }


class PresetName(str, Enum):
    COMMUNITY_BALANCED = "community_balanced"
    COMPLETE_REPLICON_RECOVERY = "complete_replicon_recovery"
    SEQUENCE_ACCURATE_ASSEMBLY = "sequence_accurate_assembly"
    CUSTOM = "custom"


@dataclass(frozen=True)
class Preset:
    name: PresetName
    title: str
    description: str
    weights: WeightsConfig
    gates: GateConfig


PRESETS: Dict[PresetName, Preset] = {
    PresetName.COMMUNITY_BALANCED: Preset(
        name=PresetName.COMMUNITY_BALANCED,
        title="Community-balanced",
        description="Rounded priorities from 22 microbial genomics community survey respondents.",
        weights=WeightsConfig(accuracy=28.0, contiguity=20.0, residual=17.0, replicon=35.0),
        gates=GateConfig(complete_recovery=False, zero_residual_hits=False),
    ),
    PresetName.COMPLETE_REPLICON_RECOVERY: Preset(
        name=PresetName.COMPLETE_REPLICON_RECOVERY,
        title="Complete-replicon-recovery",
        description="Requires ≥95% coverage on all reference replicons, prioritizing sequence accuracy among eligible tools.",
        weights=WeightsConfig(accuracy=43.0, contiguity=31.0, residual=26.0, replicon=0.0),
        gates=GateConfig(complete_recovery=True, zero_residual_hits=False),
    ),
    PresetName.SEQUENCE_ACCURATE_ASSEMBLY: Preset(
        name=PresetName.SEQUENCE_ACCURATE_ASSEMBLY,
        title="Sequence-accurate-assembly",
        description="Emphasizes base-level sequence accuracy while penalizing errors heavily.",
        weights=WeightsConfig(accuracy=50.0, contiguity=14.0, residual=12.0, replicon=24.0),
        gates=GateConfig(complete_recovery=False, zero_residual_hits=False),
    ),
}


def validate_weights(weights: WeightsConfig) -> None:
    vals = [weights.accuracy, weights.contiguity, weights.residual, weights.replicon]
    if any(w < 0.0 for w in vals):
        raise InvalidWeightsError(f"Weights must be non-negative: {weights}")
    total = sum(vals)
    if abs(total - 100.0) > 1e-4:
        raise InvalidWeightsError(
            f"Weights must total exactly 100%. Provided sum: {total:.2f}% (difference {total - 100.0:+.2f}%)"
        )
    if not any(w > 0.0 for w in vals):
        raise InvalidWeightsError("At least one criterion weight must be greater than zero.")


def get_preset(name: Union[PresetName, str]) -> Preset:
    if isinstance(name, str):
        try:
            p_name = PresetName(name)
        except ValueError:
            raise ValueError(f"Unknown preset name: '{name}'. Supported: {[p.value for p in PresetName]}")
    else:
        p_name = name
    if p_name not in PRESETS:
        raise ValueError(f"Preset '{p_name}' does not have fixed defaults (e.g. custom).")
    return PRESETS[p_name]
