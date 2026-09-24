# Scoring v1 Robustness Evidence Summary

## Provenance & Metadata
- **Scoring Version**: 1.0
- **Benchmark Source Commit**: `4d6b8cb1d482e5066f9ca575ebd3b67af4a32562`
- **Benchmark Source Data Hash (SHA-256)**: `ee7a50a277fef6153a1cc05f9cc2afd227e940efbd36a4c2b0abe9a57580901b`
- **Survey Data Hash (SHA-256)**: `3e7cfddea1d70f28dbf0ec0a605e0a0e4201d8a77c82fc3d93e97e30380c297f`
- **Supported Scenarios**: HAC 20×, HAC 100×, SUP 20×, SUP 100×

## Methodological Scope and Boundary Disclosures
- **Descriptive Non-Inferential Scope**: Respondent-specific results are explicitly descriptive of the 22 observed responses in the microbial QC survey. They make no inference that these 22 respondents represent a wider population of bioinformaticians or microbiologists.
- **No Re-survey**: No respondents were re-contacted and no new survey was required.
- **Superseded Formulations Excluded**: This analysis does not compare against the superseded dynamic min–max or geometric mean scoring systems.
- **Fixed Warning Threshold**: The accepted 25-point isolate-variability warning threshold is preserved unchanged.
- **Outcome Interpretation**: A changed winner under alternative denominators, sub-cohorts, or preferences is reported as an empirical robustness result, not treated as a failed test.

## 1. Replicon Calibration Denominator Sensitivity ($B \in \{2, 3, 4\}$)
Rankings were recomputed with replicon-loss denominators $B = 2.0, 3.0, 4.0$ while holding every other scoring-v1 decision fixed under Community-balanced weights (accuracy: 28%, contiguity: 20%, residual: 17%, replicon: 35%).

| Scenario | Winner ($B=2$) | Winner ($B=3$, Canonical) | Winner ($B=4$) | Winner Stable? | Max Rank Shift |
| :--- | :--- | :--- | :--- | :---: | :---: |
| HAC 20x | `seqkit-dorado` | `seqkit-dorado` | `seqkit-dorado` | Yes | 2 |
| HAC 100x | `seqkit-barbell` | `seqkit-barbell` | `seqkit-barbell` | Yes | 0 |
| SUP 20x | `chopper-porechop_abi` | `chopper-porechop_abi` | `chopper-porechop_abi` | Yes | 0 |
| SUP 100x | `chopper-porechop_abi` | `chopper-porechop_abi` | `chopper-porechop_abi` | Yes | 3 |

## 2. Leave-One-Isolate-Out Stability (13 Isolates)
Rankings were recomputed leaving out one of the 13 reference isolates at a time for all 4 supported scenarios (52 sub-cohort evaluations total).

| Scenario | Canonical Winner | LOO Winner Match Frequency | Winner Stable Across All 13? | Max Rank Shift Across Combinations |
| :--- | :--- | :---: | :---: | :---: |
| HAC 20x | `seqkit-dorado` | 13/13 (100.0%) | Yes | 11 |
| HAC 100x | `seqkit-barbell` | 10/13 (76.9%) | No | 16 |
| SUP 20x | `chopper-porechop_abi` | 13/13 (100.0%) | Yes | 7 |
| SUP 100x | `chopper-porechop_abi` | 13/13 (100.0%) | Yes | 12 |

## 3. Observed Respondent Weight Profiles (22 Respondents)
The 22 observed respondent preference profiles from `microbial-qc-survey.csv` were applied directly to the canonical cohort criterion scores.

| Scenario | Leading Combination(s) | Winner Frequency among Observed Respondents |
| :--- | :--- | :--- |
| HAC 20x | `seqkit-dorado`: 22/22 (100.0%) |
| HAC 100x | `seqkit-barbell`: 21/22 (95.5%), `unprocessed-barbell`: 1/22 (4.5%) |
| SUP 20x | `chopper-porechop_abi`: 21/22 (95.5%), `unprocessed-porechop_abi`: 1/22 (4.5%) |
| SUP 100x | `chopper-porechop_abi`: 21/22 (95.5%), `seqkit-porechop_abi`: 1/22 (4.5%) |

### Top Combinations Rank Variation Across Observed Profiles

| Scenario | Combination | Min Rank | Max Rank | Mean Rank |
| :--- | :--- | :---: | :---: | :---: |
| HAC 20x | `seqkit-dorado` | 1 | 1 | 1.00 |
| HAC 20x | `chopper-dorado` | 2 | 3 | 2.05 |
| HAC 20x | `unprocessed-porechop_abi` | 2 | 3 | 2.95 |
| HAC 20x | `unprocessed-dorado` | 4 | 5 | 4.14 |
| HAC 20x | `seqkit-barbell` | 5 | 10 | 6.09 |
| HAC 100x | `seqkit-barbell` | 1 | 3 | 1.09 |
| HAC 100x | `seqkit-dorado` | 2 | 5 | 2.14 |
| HAC 100x | `chopper-barbell` | 3 | 4 | 3.05 |
| HAC 100x | `chopper-dorado` | 4 | 6 | 4.09 |
| HAC 100x | `seqkit-porechop_abi` | 2 | 5 | 4.86 |
| SUP 20x | `chopper-porechop_abi` | 1 | 3 | 1.09 |
| SUP 20x | `chopper-dorado` | 2 | 4 | 2.09 |
| SUP 20x | `unprocessed-dorado` | 2 | 3 | 2.95 |
| SUP 20x | `unprocessed-porechop_abi` | 1 | 4 | 3.86 |
| SUP 20x | `fastplong-all` | 5 | 7 | 5.32 |
| SUP 100x | `chopper-porechop_abi` | 1 | 3 | 1.09 |
| SUP 100x | `fastplong-all` | 2 | 2 | 2.00 |
| SUP 100x | `unprocessed-dorado` | 3 | 4 | 3.05 |
| SUP 100x | `seqkit-porechop_abi` | 1 | 4 | 3.86 |
| SUP 100x | `chopper-barbell` | 5 | 5 | 5.00 |

## Artifact Manifest
- `denominator_sensitivity.csv`: 204 records covering all 4 scenarios × 3 denominators × 17 combinations.
- `leave_one_out.csv`: 884 records covering all 4 scenarios × 13 isolate omissions × 17 combinations.
- `respondent_profiles.csv`: 1,496 records covering all 4 scenarios × 22 observed survey respondents × 17 combinations.
- `evidence_summary.json`: Machine-readable structured provenance and stability summaries.
- `summary.md`: Human-readable summary report (this file).
