from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Union
import yaml

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent.parent / "config" / "dashboard_config.yaml"


@dataclass(frozen=True)
class DashboardConfig:
    repository_owner: str
    repository_name: str
    site_base_url: str
    default_branch: str
    source_data_path: str
    pinned_commit: str
    bundled_data_path: str
    scoring_version: str

    def provenance_text(self) -> str:
        commit_short = self.pinned_commit[:7] if self.pinned_commit else "unknown"
        return f"Source: {self.repository_owner}/{self.repository_name} @ {commit_short}"

    def data_url(self) -> str:
        return f"https://raw.githubusercontent.com/{self.repository_owner}/{self.repository_name}/{self.pinned_commit}/{self.source_data_path}"


def load_dashboard_config(config_path: Optional[Union[str, Path]] = None) -> DashboardConfig:
    path = Path(config_path) if config_path else DEFAULT_CONFIG_PATH
    if not path.exists():
        raise FileNotFoundError(f"Dashboard configuration file not found: {path}")

    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}

    repo = data.get("repository", {})
    evidence = data.get("evidence", {})
    scoring = data.get("scoring", {})

    return DashboardConfig(
        repository_owner=str(repo.get("owner", "adimascf")),
        repository_name=str(repo.get("name", "biox7021")),
        site_base_url=str(repo.get("base_url", "https://adimascf.github.io/biox7021/")),
        default_branch=str(repo.get("default_branch", "main")),
        source_data_path=str(evidence.get("source_data_path", "logbook/assembly_metrics.csv")),
        pinned_commit=str(evidence.get("pinned_commit", "4d6b8cb1d482e5066f9ca575ebd3b67af4a32562")),
        bundled_data_path=str(evidence.get("bundled_data_path", "assets/data/assembly_metrics.csv")),
        scoring_version=str(scoring.get("version", "1.0")),
    )
