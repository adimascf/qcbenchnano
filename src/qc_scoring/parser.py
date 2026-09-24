from dataclasses import dataclass
import re
from typing import Any, List, Mapping, Tuple
from qc_scoring.validation import InputValidationError, ScoringError

REPLICON_PATTERN = re.compile(
    r"^\s*([^(]+?)\s*\(\s*(\d+)\s*bp\s*,\s*([\d.]+)\s*%\s*cov\s*\)\s*$"
)


class RepliconParsingError(ScoringError):
    """Raised when replicon coverage text cannot be parsed or values are out of bounds."""
    pass


class RepliconDisagreementError(InputValidationError):
    """Raised when parsed replicon coverage counts disagree with numeric count columns."""
    pass


@dataclass(frozen=True)
class RepliconInfo:
    name: str
    size_bp: int
    coverage_pct: float

    @property
    def is_full_missed(self) -> bool:
        return self.coverage_pct == 0.0

    @property
    def is_partial_missed(self) -> bool:
        return 0.0 < self.coverage_pct < 50.0

    @property
    def is_below_50(self) -> bool:
        return self.coverage_pct < 50.0

    @property
    def is_complete_recovery(self) -> bool:
        return self.coverage_pct >= 95.0


def parse_replicon_coverage(text: str) -> List[RepliconInfo]:
    if not isinstance(text, str) or not text.strip():
        raise RepliconParsingError("Replicon coverage text is empty or not a string")

    parts = [p.strip() for p in text.split(";") if p.strip()]
    if not parts:
        raise RepliconParsingError("No valid replicon segments found in coverage text")

    replicons: List[RepliconInfo] = []
    for part in parts:
        match = REPLICON_PATTERN.match(part)
        if not match:
            raise RepliconParsingError(
                f"Malformed replicon text segment '{part}'. Expected 'name (Nbp, X.X% cov)'"
            )

        name = match.group(1).strip()
        size_bp = int(match.group(2))
        coverage_pct = float(match.group(3))

        if coverage_pct < 0.0 or coverage_pct > 100.0:
            raise RepliconParsingError(
                f"Coverage percentage {coverage_pct}% is out of valid range [0.0, 100.0]"
            )

        replicons.append(RepliconInfo(name=name, size_bp=size_bp, coverage_pct=coverage_pct))

    return replicons


def validate_and_parse_replicon_row(row: Mapping[str, Any]) -> Tuple[List[RepliconInfo], bool]:
    cov_text = str(row.get("all_contigs_coverage", ""))
    replicons = parse_replicon_coverage(cov_text)

    full_missed = int(row.get("full_missed", 0))
    partial_missed = int(row.get("partial_missed", 0))
    total_missed = int(row.get("total_missed", 0))

    derived_full = sum(1 for r in replicons if r.is_full_missed)
    derived_partial = sum(1 for r in replicons if r.is_partial_missed)
    derived_total = derived_full + derived_partial

    if derived_full != full_missed:
        raise RepliconDisagreementError(
            f"Disagreement for full missed: coverage text has {derived_full}, column has {full_missed}"
        )
    if derived_partial != partial_missed:
        raise RepliconDisagreementError(
            f"Disagreement for partial missed: coverage text has {derived_partial}, column has {partial_missed}"
        )
    if derived_total != total_missed:
        raise RepliconDisagreementError(
            f"Disagreement for total missed: coverage text has {derived_total}, column has {total_missed}"
        )

    is_complete_recovery = all(r.is_complete_recovery for r in replicons)
    return replicons, is_complete_recovery
