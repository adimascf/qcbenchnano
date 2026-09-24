# BIOX7021: Evaluating Quality Control and Trimming Tools for Long-Read Bacterial Sequencing

This repository contains the benchmark code, analysis workflows, and interactive recommendation system for the final research project (BIOX7021) in the Master of Bioinformatics at The University of Queensland.

## Recommendation System & Scoring Framework

- **Scoring Specification:** `1.0`
- **Specification Issue:** [GitHub Issue #10](https://github.com/adimascf/biox7021/issues/10)
- **Specification Document:** [`docs/qc-trimming-recommendation-system-specification.md`](docs/qc-trimming-recommendation-system-specification.md)
- **Interactive Dashboard:** [`tool_weighting.qmd`](tool_weighting.qmd) ([Rendered HTML](tool_weighting.html))
- **Methods & Robustness Documentation:** [`methods_and_robustness.qmd`](methods_and_robustness.qmd) ([Rendered HTML](methods_and_robustness.html))
- **Dashboard Source & Architecture:** [`docs/dashboard.md`](docs/dashboard.md)

### Key Features
- **Four Fixed Value Functions:** Evaluates reference-aware contiguity (symmetric auNGA), sequence accuracy (error event rate per 100kbp), residual adapter/barcode removal (% clean isolates), and replicon recovery (< 50% severe loss).
- **No Dynamic Min–Max or Geometric Scoring:** Replaces superseded dynamic scaling with four transparent, fixed 0–100 rulers.
- **Scenario Independence:** Evaluates HAC 20×, HAC 100×, SUP 20×, and SUP 100× within independent experimental contexts (rejects 50× and aggregate scenarios).
- **Presets & Custom Trade-offs:** Supports Community-balanced, Complete-replicon-recovery, Sequence-accurate-assembly, and Custom weight allocations.
- **Eligibility Gates:** Optional non-compensatory gates for complete replicon recovery ($\ge 95\%$) and zero residual hits.
- **Diagnostic Warnings:** Automated flags for isolate variability ($\ge 25$ pt gap), residual hits, and replicon losses.
- **Ordered CSV Export:** 40-column reproducible export preserving visible table order, metadata, and data provenance.
- **Robustness Evidence:** Pinned repository analyses for denominator sensitivity, leave-one-out stability, and observed respondent profiles in `results/robustness/`.

## Repository Structure

```
├── config/
│   └── dashboard_config.yaml  # Configured repository metadata & data paths
├── src/
│   └── qc_scoring/            # Canonical Python scoring library (models, scorer, criteria, paper, robustness)
├── assets/
│   └── data/                  # Pinned benchmark metrics (884 observations at 4d6b8cb)
├── results/
│   └── robustness/            # Repository-only robustness outputs (CSVs, JSON, summary.md)
├── docs/                      # Technical specification, dashboard architecture, and agent guides
├── tests/                     # Comprehensive unit, parity, and browser journey test suite
├── tool_weighting.qmd         # Human-editable Shinylive dashboard source
└── methods_and_robustness.qmd # Human-editable Methods & Robustness reference page
```

## Verification & Build

```bash
# Run full test suite (unit, parity, and browser journey tests)
pytest

# Typecheck scoring library
mypy src/qc_scoring

# Render documentation and dashboard with Quarto
quarto render
```
