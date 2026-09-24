# QC–trimming recommendation system specification

**Scoring specification:** 1.0

**Parent Specification:** [GitHub Issue #10](https://github.com/adimascf/biox7021/issues/10)

**Status:** Proposed for coauthor review

**Benchmark scope:** 13 bacterial reference isolates, 17 QC–trimming combinations, HAC/SUP basecalling, and 20×/100× depth

**Source-data baseline:** 884 complete scenario–isolate–combination observations at project commit `4d6b8cb1d482e5066f9ca575ebd3b67af4a32562`

## Problem Statement

Researchers choosing Oxford Nanopore read quality-control and adapter-trimming methods face genuine trade-offs. A combination that improves sequence accuracy may lose small, low-copy plasmids; another may preserve all reference replicons but retain adapter or barcode sequence. No single combination is therefore best independently of the researcher's priorities.

The current analysis does not yet provide a methodologically stable decision aid. It rescales some metrics relative to whichever alternatives happen to be present and combines them geometrically. Consequently, the same raw result can receive a different score after filtering, a small difference can be magnified into the full scoring range, and the geometric zero workaround introduces an undocumented balance preference. The paper and published dashboard also implement different grouping and aggregation sequences, so they can disagree despite ostensibly using the same metrics and weights.

The project needs one transparent scoring specification that:

- preserves the biological meaning of four distinct assembly criteria;
- maps their indicators to comparable, fixed 0–100 preference scales;
- separates compensatory preferences from non-compensatory requirements;
- supports the rounded community-survey priorities without presenting them as universal values;
- exposes important weaknesses and isolate-level variation rather than hiding them in a mean;
- gives an individual researcher a simple way to apply different priorities; and
- produces identical results in the paper and dashboard from pinned benchmark evidence.

The recommendation scope is deliberately narrow. It covers only the 13 named bacterial reference isolates and 884 observations in the current benchmark. It must not imply that the resulting order is a universal ranking of QC or adapter-trimming software.

## Solution

Provide a versioned, repository-name-independent recommendation system built around one canonical scoring implementation. A user first selects the basecalling model and sequencing depth that match their experiment. The system then evaluates all 17 QC–trimming combinations in that scenario using four fixed criterion scores: reference-aware contiguity, sequence accuracy, residual adapter/barcode removal, and reference-replicon recovery.

The user may start from a clearly described preset or enter four percentages that total 100. A normalized weighted arithmetic mean produces a preference-alignment score. Eligibility gates are applied separately so that an unacceptable result cannot be compensated for by strong performance elsewhere. The dashboard initially presents the five highest-scoring eligible combinations, while preserving access to all evaluated and excluded combinations.

Every recommendation presents the raw evidence behind its normalized scores, identifies detected failures, and provides isolate-level details. The paper retains its raw metric-specific figures and separately reports the Community-balanced composite result for each supported scenario. A Methods & robustness page documents the equations, assumptions, provenance, and repository-only stability checks.

## User Stories

1. As a bacterial genomics researcher, I want to state whether my reads were HAC- or SUP-basecalled, so that I am not shown evidence from an irrelevant basecalling model.
2. As a bacterial genomics researcher, I want to state whether my data resemble 20× or 100× depth, so that the recommendation matches an evaluated sequencing-depth scenario.
3. As a dashboard user, I want both scenario selections to begin unset, so that I cannot accidentally interpret a default scenario as my own.
4. As a dashboard user, I want unsupported scenarios such as 50× to be absent, so that every offered choice corresponds to complete benchmark evidence.
5. As a researcher seeking a general starting point, I want a Community-balanced preset, so that I can use the rounded priorities submitted in the community survey.
6. As a researcher who cannot tolerate reference-replicon loss, I want a Complete-replicon-recovery preset, so that incomplete combinations are excluded before ranking.
7. As a researcher prioritising base-level fidelity, I want a Sequence-accurate-assembly preset, so that sequence errors dominate the recommendation without erasing other consequences.
8. As a researcher with different downstream needs, I want Custom priorities, so that the survey average does not dictate my decision.
9. As a custom-priority user, I want to set any criterion to zero, so that a consequence irrelevant to my use case has no effect on my overall score.
10. As a custom-priority user, I want my four percentages to total exactly 100, so that their relative influence is explicit and communicable.
11. As a custom-priority user, I want to see how far my allocation is below or above 100, so that I can correct an invalid entry.
12. As a custom-priority user, I want ranking to pause while my percentages are invalid, so that a partial allocation is not silently reinterpreted.
13. As a custom-priority user, I want optional complete-recovery and zero-residual-hit requirements, so that I can express non-compensatory constraints separately from preferences.
14. As a preset user, I want any manual edit to switch the state visibly to Custom, so that I know I am no longer using the published preset unchanged.
15. As a dashboard user, I want a reset-to-community action, so that I can easily return to the published reference allocation.
16. As a dashboard user, I want the recommendation list ordered by overall preference alignment, so that the primary result is simple to interpret.
17. As a dashboard user, I want the five highest-scoring eligible combinations shown first, so that I am not overwhelmed by 17 alternatives.
18. As a dashboard user, I want to reveal all evaluated combinations, so that the shortlist does not hide evidence.
19. As a dashboard user, I want excluded combinations retained in a collapsed section with reasons, so that absence is not mistaken for lack of evaluation.
20. As a dashboard user, I want the system to preserve my gates when no combination is eligible, so that it does not recommend an option that violates my requirements.
21. As a dashboard user, I want scores displayed to one decimal place but calculated at full precision, so that presentation is readable without changing the ranking.
22. As a dashboard user, I want adjacent results less than one point apart identified as similar overall scores, so that a small numerical difference is not overstated.
23. As a dashboard user, I want all four criterion scores visible even when one has zero weight, so that I can see a combination's unselected weaknesses.
24. As a dashboard user, I want the mean combined error-event rate shown beside the sequence-accuracy score, so that the normalized value remains interpretable.
25. As a dashboard user, I want mismatches and indels shown separately in supporting detail, so that an indel-heavy result is not concealed by their combined rate.
26. As a dashboard user, I want the mean auNGA ratio shown beside the contiguity score, so that I can see the departure from the reference target.
27. As a dashboard user, I want total residual adapter/barcode hits and affected isolates reported together, so that concentrated and widespread events are distinguishable.
28. As a dashboard user, I want full and partial reference-replicon losses reported separately, so that complete absence is distinguishable from severe incompleteness.
29. As a dashboard user, I want the absolute number of missed or severely incomplete replicons reported, so that the loss is not diluted by an isolate's total replicon count.
30. As a dashboard user, I want concrete warnings for detected hits, replicon losses, and isolate-specific weakness, so that a high overall score cannot conceal a material failure.
31. As a dashboard user, I want labelled per-isolate accuracy and contiguity results, so that I can identify which reference genome drives inconsistency.
32. As a dashboard user, I want sparse event details in a compact table, so that contamination and replicon failures are legible without a misleading continuous plot.
33. As a dashboard user, I want no qualitative labels such as good or excellent, so that the preference-alignment score is not mistaken for a validated quality grade.
34. As a dashboard user, I want an explanation of the benchmark's limited scope, so that I do not generalise its recommendations to all bacterial isolates or sequencing protocols.
35. As a dashboard user, I want a separate Methods & robustness page, so that I can inspect the equations and assumptions without cluttering the main decision workflow.
36. As a dashboard user, I want to download all 17 results in their displayed order, so that my exported record matches what I saw.
37. As a dashboard user, I want the export to include scenario, weights, gates, scores, ranks, warnings, and provenance, so that I can reproduce and explain the decision.
38. As a paper reader, I want raw metric-specific figures retained, so that the empirical observations remain visible independently of the preference model.
39. As a paper reader, I want each metric panel ordered appropriately for its own scenario, so that the clearest raw-metric presentation is not constrained by dashboard ordering.
40. As a paper reader, I want the Community-balanced composite results reported separately, so that I can distinguish empirical measurements from survey-weighted decision analysis.
41. As a paper reader, I want community weights described as rounded priorities from 22 respondents, so that they are not presented as objective or population-representative truths.
42. As a coauthor, I want paper and dashboard results generated through one scoring contract, so that the two outputs cannot silently disagree.
43. As a coauthor, I want the scoring rules versioned independently of cosmetic site changes, so that methodological changes are traceable.
44. As a coauthor, I want the source evidence pinned to a commit or content hash, so that rankings do not change when a repository branch advances.
45. As a maintainer, I want repository ownership and naming supplied as deployment configuration, so that a repository rename does not break data loading or provenance links.
46. As a maintainer, I want malformed or incomplete benchmark inputs rejected explicitly, so that an invalid combination cannot receive an apparently valid ranking.
47. As a maintainer, I want the human-readable replicon coverage field cross-checked against its count columns, so that parsing errors cannot silently alter gates or scores.
48. As a maintainer, I want the dashboard source version-controlled rather than retaining only compiled output, so that future changes can be reviewed and reproduced.
49. As an implementer, I want behavioural acceptance tests at the canonical scoring boundary, so that I can modify internals without changing agreed outcomes.
50. As an implementer, I want small canonical fixtures shared by the paper and dashboard checks, so that grouping, clipping, gates, and rounding remain identical.
51. As a project reviewer, I want repository-only robustness reports for the selected modelling decisions, so that assumptions can be inspected without expanding the main paper.
52. As a project reviewer, I want observed respondent-specific rankings described without inferential population claims, so that preference heterogeneity is represented honestly.

## Implementation Decisions

1. **Supported evidence domain.** Version 1.0 accepts exactly the evaluated scenario dimensions: `hac` or `sup` basecalling and `20x` or `100x` sequencing depth. The reference dataset comprises 17 QC–trimming combinations and 13 named bacterial reference isolates in each scenario, for 884 unique observations. The canonical observation key is scenario, isolate, and QC–trimming combination.

2. **Required scenario selection.** Basecalling model and depth have no initial defaults. Ranking begins only after both have been selected. After a valid scenario is selected, Community-balanced becomes the initial preference preset.

3. **Canonical scoring boundary.** One pure, tested Python scoring module owns input validation, replicon-coverage parsing, criterion value functions, cohort aggregation, presets, gates, warnings, ranks, near-tie status, provenance, and export-ready records. It accepts the benchmark table plus an explicit scenario and preference configuration and returns the complete ordered result set. Paper analyses call this boundary directly. The dashboard imports the same implementation where practical; otherwise it consumes canonical criterion outputs and must pass parity fixtures produced at the same boundary.

4. **Fixed value functions.** Criterion transformations are fixed across all scenarios, filters, presets, and displayed alternatives. They are never recalculated from a selected group's minimum or maximum. This gives a raw value the same interpretation in HAC/SUP and 20×/100× results and prevents rank reversal caused only by adding or removing alternatives.

5. **Reference-aware contiguity.** For isolate auNGA ratio \(r\), the criterion score is

   \[
   S_{\mathrm{contiguity}}(r)=100\max(0,1-|r-1|).
   \]

   Equal departures above and below the reference target of 1.0 receive equal penalties. This is a deliberate project value judgement, not a claim that fragmentation and duplication have identical biological mechanisms. Supporting diagnostics such as duplication ratio and misassemblies remain available for interpretation but do not create additional primary sliders.

6. **Sequence accuracy.** Let \(M\) and \(I\) be QUAST mismatch and indel events per 100 kbp of aligned assembly sequence. The combined error-event rate is \(E=M+I\), and

   \[
   S_{\mathrm{accuracy}}(M,I)=100\,\operatorname{clip}\left(1-\frac{M+I}{10},0,1\right).
   \]

   Zero events score 100 and ten or more combined events per 100 kbp score 0. Each reported event receives equal primary-score treatment because the project has no validated biological conversion factor between mismatch and indel consequences. The dashboard retains both components and warns through supporting detail when errors are indel-heavy. The denominator ten is a frozen benchmark-policy calibration, not a universal ONT or bacterial quality threshold.

7. **Residual adapter/barcode removal.** A benchmark isolate scores 100 when it has zero resolved residual adapter/barcode hits and 0 when it has one or more. The cohort criterion score is the mean of those binary isolate scores, equivalently the percentage of the 13 isolates with zero detected hits. The raw summaries additionally report total hits and affected isolates. The paper may visualise summed hits, but summed hits are not the normalized criterion score.

8. **Residual-hit terminology and evidence.** A scored event is a resolved alignment to a defined ONT adapter or barcode region passing the active workflow thresholds of at least 90% identity and at least 90% region coverage. It is called a residual adapter/barcode hit, not generic biological contamination. The specification relies on the project's existing detection validation and does not require a new validation study.

9. **Reference-replicon recovery.** A missed or severely incomplete reference replicon has less than 50% reference coverage. A fully missed replicon has 0% coverage; a partially missed replicon has greater than 0% but less than 50% coverage. Exactly 50% is not included. For one QC–trimming combination and scenario, let \(T\) be the absolute count of qualifying replicons summed across all 13 isolates. Then

   \[
   S_{\mathrm{replicon}}(T)=100\max\left(0,1-\frac{T}{3}\right).
   \]

   Totals of zero, one, two, and at least three therefore score 100, 66.7, 33.3, and 0 before display rounding. Each missed replicon has the same effect regardless of how many replicons its isolate contains. The dashboard also reports affected isolates plus full and partial losses.

10. **Existing replicon-data contract.** No upstream pipeline schema change is required. The canonical scorer uses the numeric full, partial, and total missed columns for the primary score and parses the existing semicolon-delimited all-replicon coverage field for complete-recovery gates. It validates the text structure and verifies that coverage-derived below-50% counts equal the supplied count columns. Malformed text, impossible coverage, or disagreement is an explicit data error.

11. **Cohort aggregation.** Contiguity and accuracy are scored per isolate and then averaged across the 13-isolate cohort. Residual-hit removal is the mean binary clean-isolate score. Replicon recovery is intentionally the scenario-level absolute-loss score rather than a mean fraction recovered. The four resulting values are cohort-level criterion scores.

12. **Weighted arithmetic aggregation.** For active criterion weights \(w_k\) and cohort criterion scores \(S_k\), the preference-alignment score is

   \[
   S_{\mathrm{overall}}=\frac{\sum_k w_kS_k}{\sum_k w_k}.
   \]

   Custom weights must be non-negative percentages totalling exactly 100, and at least one criterion must be non-zero. The arithmetic model makes trade-offs explicit. Non-compensatory requirements are represented by eligibility gates rather than by a geometric mean or epsilon/shift workaround.

13. **Community-balanced preset.** The preset uses the rounded survey allocations: 28% sequence accuracy, 20% reference-aware contiguity, 17% residual adapter/barcode removal, and 35% reference-replicon recovery. These sum to 100. The preset has no automatic gate and is described as the priorities of the surveyed group, not an objective optimum.

14. **Complete-replicon-recovery preset.** Every expected reference replicon must have at least 95% reference coverage. Coverage alone does not claim circularisation or structural validation. Among eligible combinations, ranking uses 43% sequence accuracy, 31% reference-aware contiguity, and 26% residual adapter/barcode removal. Replicon recovery is enforced by the gate and is not counted again as a ranking advantage.

15. **Sequence-accurate-assembly preset.** This preset uses 50% sequence accuracy, 14% reference-aware contiguity, 12% residual adapter/barcode removal, and 24% reference-replicon recovery. It has no additional automatic gate.

16. **Custom-priorities preset.** The user enters four percentages. Both optional advanced gates are off initially. The complete-recovery gate applies the same 95%-per-replicon rule as the corresponding preset. The zero-residual-hit gate requires zero detected hits across all 13 isolates. Editing any preset weight or gate changes the visible mode to Custom; selecting a named preset replaces all weights and gate states.

17. **Eligibility before ranking.** Gates are evaluated before scores are ordered. Failed combinations are not assigned a rank. If every combination fails, the system preserves the user's constraints, reports that no combination meets them, and lists the reasons. It never weakens a gate automatically.

18. **Missing and invalid evidence.** Each scenario–combination requires all 13 named isolates and every indicator needed by an active or displayed criterion. A combination with incomplete evidence is marked `insufficient benchmark data`, receives no rank, and appears with an explicit reason. Values are not imputed and weights are not redistributed. Invalid values, duplicate observation keys, unexpected scenario values, and mismatched replicon representations fail validation explicitly.

19. **Ranking precision and near ties.** Eligible combinations are ordered by the full-precision preference-alignment score. Scores are displayed to one decimal place. Adjacent scores differing by less than 1.0 are labelled `similar overall scores`; this is a communication qualifier, not a statistical-equivalence claim. Exact score ties share the same rank, and their stable display order is determined by the QC–trimming identifier solely for reproducibility.

20. **Warnings.** Every residual-hit event and every missed or severely incomplete replicon generates a concrete count-and-prevalence warning. Accuracy or contiguity is labelled variable when the lowest isolate score is at least 25 points below its cohort mean. Warnings do not create a second ranking method and remain visible even if their criterion has zero weight.

21. **Raw summaries.** Each result shows the four normalized criterion scores and the underlying evidence: combined errors per 100 kbp with mismatch and indel components; mean auNGA ratio; clean and affected isolate counts plus total residual hits; and total, full, and partial replicon losses plus affected isolates. Supporting duplication ratio and misassembly information may appear in detail without becoming primary criteria.

22. **Dashboard interaction.** The main experience is one guided page: required scenario selection; preset and percentage allocation; ranked shortlist; selected-combination details; and a compact methodology disclosure. The initial list contains the five highest-ranked eligible alternatives. Users may reveal all 17. Ineligible and insufficient-data alternatives appear afterward in a collapsed section with no rank and an explicit reason.

23. **Fixed display order.** The recommendation list is not re-sortable by individual metrics; rank always means the overall preference-alignment order. Criterion-specific comparisons are available in details. This also gives the export a single unambiguous order.

24. **Per-isolate detail.** Accuracy and contiguity use labelled per-isolate dot plots with their cohort mean. Residual adapter/barcode and replicon outcomes use a compact per-isolate event table. Population confidence intervals are not shown because these 13 named reference isolates are the benchmark cohort rather than a random sample from a defined bacterial population.

25. **Score language.** The overall number is called a preference-alignment score within the selected benchmark scenario. It receives no labels such as excellent, good, or poor and is never described as a probability, universal tool quality, or externally validated assembly-quality grade.

26. **CSV export.** The export contains all 17 combinations in exactly the current display order. Eligible combinations come first with ranks; excluded or insufficient-data combinations follow with blank ranks and reasons. It includes raw summaries, criterion scores, full-precision overall score, displayed score, near-tie status, warnings, selected scenario, preset, four weights, active gates, scoring version, source-data version or content hash, and generation date.

27. **Methods & robustness page.** A separate site page documents scope, terminology, equations, policy anchors, survey derivation, presets, gates, warnings, data provenance, limitations, and the repository-only robustness results. The decision page contains only a concise methods summary and link.

28. **Evidence provenance and repository independence.** Dashboard evidence is bundled or pinned to a specific source commit or content hash. Repository owner, repository name, raw-data location, and site base URL are build/deployment configuration rather than embedded scoring assumptions. A repository rename must not change scores or require editing the canonical methodology.

29. **Scoring version.** The initial accepted rules are version 1.0. A change to a value function, calibration denominator, cohort aggregation, preset, gate, warning threshold, input interpretation, or ranking semantics creates a new scoring version. Styling, explanatory copy, or other cosmetic site changes do not.

30. **Paper presentation.** Raw metric-specific figures remain primary empirical evidence and may independently order their x-axes within scenario to aid interpretation. The Community-balanced composite is reported separately for all four supported scenarios as a survey-weighted decision analysis, with the leading combinations and their four criterion scores. The paper does not describe the dashboard's customizable output as a single universal result.

31. **Survey heterogeneity analysis.** The existing 22 respondent weight allocations are each applied descriptively to the canonical cohort criterion scores. The repository reports how often each combination ranks first and how ranks vary across the observed profiles. This requires no new survey and makes no inferential claim about a wider population.

32. **Repository-only robustness checks.** Store reproducible outputs for replicon calibration denominators \(B=2,3,4\), leave-one-isolate-out rankings, the 22 observed respondent profiles, and paper/dashboard numerical parity. These checks need not appear in the paper or supplement. Do not compare against the superseded dynamic min–max/geometric system, and do not vary the accepted 25-point warning threshold.

33. **Version-controlled dashboard source.** Maintainable dashboard and Methods & robustness source must live in version control. Compiled site output alone is not an implementation source of truth. The implementation may retain Quarto/Shinylive if it satisfies the canonical scoring and export contracts; the specification does not mandate a framework migration.

34. **Accessible presentation.** Criterion information is not conveyed by colour alone. Controls have explicit labels, keyboard operation, meaningful validation messages, and readable values. Plots and tables have text equivalents sufficient to interpret warnings and raw evidence.

## Testing Decisions

1. **Primary test seam.** The highest-value scoring seam is the canonical scorer's public operation: validated benchmark observations plus scenario and preference configuration in, complete ordered recommendation records out. Tests assert externally observable values, eligibility, warnings, ranks, and provenance without coupling to helper functions or dataframe implementation details.

2. **Canonical real-data acceptance test.** The pinned benchmark contains 884 unique observations representing 17 combinations × 13 isolates × two basecalling models × two depths. An end-to-end test validates this contract and produces 17 ordered result records for each supported scenario under Community-balanced preferences.

3. **Fixed-score examples.** Through the public scorer boundary, tests cover contiguity ratios at 1.0, 0.9, 1.1, and values whose distance clips to zero; combined accuracy event rates at 0, intermediate values, 10, and above 10; zero versus non-zero residual hits; and replicon totals of 0, 1, 2, 3, and above 3.

4. **Order-of-operations tests.** Fixtures distinguish scoring isolates before cohort aggregation from transforming already-averaged raw values, particularly where clipping occurs. Replicon recovery separately verifies its intentional cohort-level absolute count.

5. **Replicon parser contract.** Tests cover multiple reference sequences, decimals, exactly 0%, exactly 50%, just below 50%, exactly 95%, and just below 95%; malformed entries; out-of-range coverage; missing fields; and disagreement with full, partial, or total count columns.

6. **Preset and weight tests.** Each named preset produces its documented weights and gates. Custom percentages accept zeros and decimals but reject negative values, all-zero values, totals below 100, and totals above 100. Editing a named configuration yields visible Custom state; reset restores Community-balanced.

7. **Gate tests.** Complete recovery excludes a combination when any expected replicon is below 95%. The zero-residual-hit gate excludes a combination when any of the 13 isolates has a detected hit. Failed combinations have no rank. An all-failed scenario returns a no-eligible-recommendation state without relaxing requirements.

8. **Missing-data tests.** Missing isolates, indicators, malformed numerics, duplicate keys, unsupported scenarios, and unknown combinations produce explicit validation or insufficient-data outcomes according to whether the defect affects the entire input contract or one alternative. No test should observe imputation or silent weight redistribution.

9. **Ranking and display tests.** Ranking uses unrounded values, display uses one decimal, exact ties share rank, stable ordering is reproducible, and adjacent differences below one point receive near-tie status. Excluded records follow eligible records and carry blank ranks.

10. **Warning tests.** Residual hits and replicon losses always generate concrete event/affected-isolate warnings. Accuracy and contiguity warnings occur at a mean-to-minimum gap of exactly 25 or more and do not occur below 25. A zero-weight criterion still generates its warning.

11. **Export acceptance test.** Downloaded records exactly preserve the currently displayed order and contain all required settings, evidence, scores, qualifiers, reasons, and provenance. A repository rename or configured site-base change does not alter scoring data or break the export.

12. **Paper parity test.** For Community-balanced settings, every paper composite value and rank matches the canonical scorer for the same pinned data and scenario. Metric-specific plot ordering remains a presentation concern and does not modify criterion scores.

13. **Dashboard parity test.** Canonical fixtures are exercised through the interactive application at its highest practical boundary. A browser-level smoke journey selects a scenario and preset, verifies the leading recommendation and raw summaries, changes to a valid custom allocation, activates a gate, opens detail, reveals all alternatives, and downloads the ordered CSV.

14. **Methods-page test.** The published page exposes scoring version, source-data provenance, all four value functions, preset definitions, gates, limitations, and links to repository-only robustness outputs. It must not display obsolete min–max or geometric formulations as active alternatives.

15. **Robustness checks as executable analyses.** Denominator, leave-one-isolate-out, and respondent-profile analyses are deterministic from pinned inputs and save machine-readable results. They are evidence artifacts rather than tests that fail merely because a winner changes; automated checks instead ensure that the analyses run, cover all expected cases, and report their outputs completely.

16. **Testing prior art.** The repository currently validates scientific outputs primarily through Snakemake rules and generated tables and has shell-based determinism checks, but no general Python unit-test suite for scoring. New tests should introduce the smallest conventional Python test harness necessary to exercise the pure canonical boundary; scientific plotting scripts should consume its outputs rather than becoming a second test seam.

## Out of Scope

- Re-running or changing basecalling, read filtering, adapter trimming, assembly, QUAST, residual-hit detection, or reference-coverage generation.
- Requiring Dimas to add a normalized one-row-per-replicon pipeline table; version 1.0 works with and validates the existing results CSV.
- Conducting a new residual adapter/barcode false-positive validation study; the existing project validation is accepted.
- Surveying or re-contacting respondents, increasing the survey sample, or claiming the 22 respondents statistically represent the wider community.
- Comparing the accepted system with the superseded dynamic min–max/geometric system in the paper, supplement, site, or required repository analyses.
- Supporting 50× depth or an aggregate `all models`/`all depths` scenario.
- Introducing separate user-controlled mismatch and indel weights.
- Replacing the accepted symmetric auNGA-ratio distance with an asymmetric curve.
- Using TOPSIS, PROMETHEE, VIKOR, geometric aggregation, or another competing ranking method in the primary interface.
- Adding population confidence intervals for the 13 fixed benchmark isolates.
- Generalising the value functions or rankings beyond the current bacterial isolate benchmark without a new scoped scoring version.
- Assigning qualitative quality grades to either criterion or overall scores.
- Allowing the primary recommendation table to be re-sorted by individual metrics.
- Treating 95% reference coverage as proof of circularisation, structural correctness, or a finished assembly.
- Implementing a shareable-settings URL in the initial release; it may be added later without altering the scoring model.
- Renaming the project repository as part of this work.

## Further Notes

- The four explicit handoff assumptions have been resolved as follows: dynamic min–max scaling is replaced by fixed rulers; symmetric auNGA distance is retained deliberately; mismatch and indel indicators are combined without additional sliders while remaining visible; geometric shifts and epsilon floors are removed in favour of arithmetic trade-offs plus gates; and isolate variation is exposed descriptively rather than hidden or assigned population confidence intervals.
- QUAST's reference-aware contiguity and error indicators are distinct measurements rather than universal quality thresholds. The project calibrations at a combined error-event rate of ten and a replicon-loss total of three must therefore be described as benchmark-policy choices.
- The existing result file contains numeric full, partial, and total missed counts as well as all-reference coverage encoded as semicolon-delimited text. Across the pinned 884 rows it contains 62 full misses and five partial misses, for 67 total severe losses across 51 observations.
- The current published dashboard offers unsupported aggregate scenario choices and 50× depth, performs filter-dependent min–max scaling after averaging raw metrics, provides two geometric formulas, and hard-codes a repository raw-data URL. These behaviours are replaced by this specification.
- The current paper scoring script separately performs sample-level min–max transformations and an epsilon-floor geometric mean. It must delegate to the canonical scorer rather than retain an alternate formula.
- The published site currently contains compiled dashboard HTML, while maintainable dashboard source is not present on the current default branch. Restoring version-controlled source is an implementation prerequisite.
- Primary methodological references include the [QUAST paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3624806/), the [MIMAG/MISAG reporting standard](https://pmc.ncbi.nlm.nih.gov/articles/PMC6436528/), and [NCBI prokaryotic genome submission guidance](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/). None supplies universal value functions for this benchmark's four criteria.
- Repository names and owners are mutable. Documentation may identify the source commit used for this specification, but runtime behavior and provenance generation must use configured repository metadata or content hashes rather than assume a permanent `owner/repository` string.
