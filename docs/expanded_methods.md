# Expanded Reproducibility Methods

This document supplies reproducibility detail for the associated Research Letter and the shared public repository. It is not a formal Radiology supplement because that article type does not permit online supplementary material.

## Study Design and Data Sources

We conducted a cross-sectional ecological analysis of public country-level cancer projections and latest-reported radiotherapy resource records. Cancer projections were obtained from the International Agency for Research on Cancer Global Cancer Observatory Cancer Tomorrow application programming interface for data versions 2022 and 2024. We extracted country, cancer-site, sex, measure, year, predicted count, and population fields. The analysis used both sexes combined and incident cases. Resource data were obtained from country tables in the International Atomic Energy Agency Directory of Radiotherapy Centres (DIRAC). Source data were accessed in July 2026.

Countries were linked with three-letter International Organization for Standardization codes. Global Cancer Observatory country records without a DIRAC match were retained as resource-unknown observations. They were not assigned zero equipment. A DIRAC-reported value of zero, where present, was retained as an observed value. Duplicate country records and one-to-many joins were not permitted. The primary data version contained 186 Global Cancer Observatory country records.

The primary complete-case analysis included 150 countries with matched nonmissing megavoltage-unit records; 36 country records remained resource-unknown. The common-baseline burden-version analysis included 149 countries present in both burden versions and in the fixed resource snapshot. These analytic universes were defined separately because a country could contribute to missing-resource bounds without contributing to a complete-case threshold calculation.

Because the study used public, aggregate country observations and no individual health information, institutional review board approval and informed consent were not applicable.

## Cancer Burden Measures

The primary burden proxy was the unweighted sum of incident cases for 14 cancer sites commonly involving radiotherapy in at least some treatment pathways: lip and oral cavity, oropharynx, nasopharynx, hypopharynx, esophagus, rectum, larynx, trachea/bronchus/lung, breast, cervix uteri, corpus uteri, prostate, bladder, and brain/central nervous system. The corresponding Global Cancer Observatory site codes were 1, 3, 4, 5, 6, 9, 14, 15, 20, 23, 24, 27, 30, and 31. The proxy was specified to provide a reproducible burden measure and was not interpreted as radiotherapy indication, utilization, or unmet treatment need.

For data version 2024, selected-site cases were aggregated for 2024, 2040, and 2050. Relative projected growth was calculated as the difference between 2050 and 2024 cases divided by 2024 cases. Absolute increase was calculated as 2050 cases minus 2024 cases. The all-cancer sensitivity analysis used the same calculations after summing all cancer sites.

For country *i*, the primary growth measure was `(selected-site cases in 2050_i - selected-site cases in 2024_i) / selected-site cases in 2024_i`. Complete-case burden shares used absolute case increases, not relative growth, in the numerator and denominator.

## Radiotherapy Resource Measures

The primary resource numerator was the latest-reported national count of megavoltage photon or electron therapy units in DIRAC. The primary density denominator was projected selected-site incident cases in 2050. Density was expressed as megavoltage units per 1000 projected cases. The record year supplied by DIRAC was retained for each matched country. Registered unit counts were not adjusted for equipment operability, downtime, staffing, throughput, geographic distribution, or referral access.

For country *i*, unit density was `latest-reported megavoltage units_i x 1000 / selected-site cases in 2050_i`. Country rows with zero registered units were retained with zero density; cases per unit were left undefined rather than divided by zero.

Within the primary complete-case set, high growth was defined as relative growth at or above the 75th percentile, and lower resource density as density at or below the 25th percentile. Quantiles used linear interpolation. The resulting thresholds were 119.2% growth and 0.506 megavoltage units per 1000 projected cases. A country met both thresholds only when both criteria were satisfied. These thresholds define a relative screening group and do not represent an equipment-adequacy standard.

## Common-Baseline Burden-Version Analysis

To separate burden-data revision from different projection baselines, the main version comparison used 2025 as the baseline and 2050 as the horizon in both Global Cancer Observatory releases. It was restricted to the same 149 countries with burden estimates in both versions and matched resource records. The same latest-reported DIRAC snapshot was applied to both releases. Version-specific 75th-percentile growth and 25th-percentile density thresholds were calculated within this common set. Countries were classified as retained, new in data version 2024, no longer classified, or not classified in either version.

The version-2022 thresholds were 112.7% growth and 0.502 units per 1000 cases; the version-2024 thresholds were 110.5% and 0.500, respectively. Classification overlap was summarized by retained, new, and no-longer-classified countries and the Jaccard index. The comparison was descriptive and did not treat model versions as repeated measurements.

For transparency, the repository also preserves the earlier operational comparison in which each release used its native baseline year, 2022 or 2024. That comparison was treated as a secondary analysis because it combines burden-version and baseline-year differences.

## Sensitivity Analyses

The all-cancer analysis replaced selected-site incidence with all-cancer incidence in both the growth and density measures. The baseline-density analysis retained the primary 2024-to-2050 growth measure but divided latest-reported units by selected-site cases in 2024 instead of 2050. A stricter analysis used the 80th percentile of growth and 20th percentile of 2050 unit density. Additional repository analyses used a 2040 horizon, radiotherapy-centre density, selected-site mortality growth, and exclusion of countries with fewer than 1000 selected-site cases at baseline.

Resource-record recency was examined by restricting the complete-case set to DIRAC records dated 2023 or later. Growth and density thresholds were recalculated in the restricted set. Primary countries with earlier records were reported as ineligible for this sensitivity analysis rather than classified as negative. This restriction evaluates whether the observed core persists among countries with more recent inventory reports; it does not update older inventories or establish that recorded units remain functional.

## Missing-Resource Bounds

We used deterministic identification bounds to show how missing resource records could alter the proportion of countries meeting both thresholds. For each World Health Organization region and globally, the lower bound divided the observed number meeting both thresholds by all Global Cancer Observatory country records. The upper bound added resource-unknown countries whose projected growth exceeded the primary growth threshold to the numerator. The upper bound therefore represents an extreme scenario in which every unmatched high-growth country also meets the resource-density criterion. These intervals are not confidence intervals and do not assume that missingness is random.

For a stratum containing `N` country records, `O` observed countries meeting both thresholds, and `H` high-growth countries without a resource match, the lower bound was `O/N` and the upper bound was `(O+H)/N`. Resource-unknown countries below the growth threshold could not meet both criteria and therefore did not widen the interval.

## Statistical Analysis and Reproducibility

Analyses described the complete available country records and did not use sampling-based hypothesis tests. We reported counts, proportions, quantile thresholds, classification overlap, and Spearman rank correlations for construct diagnostics. Sensitivity analyses were interpreted by retention of primary countries and by named classification changes, not by statistical significance. Figure 1 used country-level analytic data for the map and burden-version panel and region-level identification-bound data for the interval panel. Table 1 was generated directly from the corresponding CSV source file without manual transcription.

All processing, analysis, and figure generation were performed in Python 3.12 with NumPy 2.3.5, pandas 2.3.3, Matplotlib 3.11.0, and SciPy 1.18.0. From the public repository root, the derived-data workflow runs `python scripts/analyze_additional_analyses.py` followed by `python scripts/make_radiotherapy_resource_mismatch_figure.py`. Full rebuilding additionally requires source files obtained from the original providers and then runs `python scripts/build_radiotherapy_resource_mismatch.py` before the additional-analysis and figure scripts.

The public repository contains source-data documentation, derived analytic files, analysis scripts, figure source data, checksums, and a file manifest. Raw source exports are not redistributed. The repository is available at https://github.com/ammkkr/radiotherapy-resource-mismatch-globocan2024-dirac.
