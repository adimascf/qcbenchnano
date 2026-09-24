from dataclasses import dataclass
from typing import Literal

SupportedModel = Literal["hac", "sup"]
SupportedDepth = Literal["20x", "100x"]

@dataclass(frozen=True)
class Scenario:
    model: str
    depth: str

    def __post_init__(self) -> None:
        # We allow creation, but validation will enforce supported values
        pass
