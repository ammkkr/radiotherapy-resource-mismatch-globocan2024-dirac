# Expanded Reproducibility Methods v6

This document provides reproducibility details for the Radiology Research Letter. It is intended for the public code repository and is not a formal online supplement; Radiology Research Letters do not permit online supplementary material.

## Study Design

This was a retrospective ecological country-level analysis of public aggregate cancer projections and radiotherapy-resource records. No individual-level data were used. Institutional review board approval, informed consent, and Health Insurance Portability and Accountability Act compliance were not applicable.

The analysis was designed as an updateable prioritization screen. It did not estimate patient-level radiotherapy indication, optimal utilization, treatment receipt, machine throughput, or unmet need.

## Data Sources and Versioning

Cancer projections were obtained from the International Agency for Research on Cancer Global Cancer Observatory Cancer Tomorrow application programming interface. The primary analysis used data version 2024. A burden-version analysis repeated the screen using data version 2022.

The version comparison held the same latest-reported International Atomic Energy Agency Directory of Radiotherapy Centres resource snapshot fixed. It therefore evaluates classification sensitivity to a burden-data update and is not a longitudinal comparison of radiotherapy resources.

The primary resource count was the latest-reported number of megavoltage photon/electron therapy units in the Directory of Radiotherapy Centres country record. Record years differed among countries. Counts describe registered equipment and do not establish operability, uptime, staffing, treatment-planning capability, or patient access.

Natural Earth geometry was used only to draw country boundaries in Figure 1a.

## Country Matching and Analytic Universes

Countries were matched using three-letter country codes.

- The data-version-2024 country universe contained 186 country records.
- Resource records matched 150 countries; 36 had no resource match.
- The burden-version comparison contained 149 countries with matched resources and burden projections in both versions.
- Countries without a resource match remained resource-unknown. They were not assigned zero units or an imputed density.

## Selected-Site Incidence Proxy

The prespecified proxy summed incident cases from 14 cancer sites commonly involving radiotherapy in at least some treatment pathways: lip/oral cavity, oropharynx, nasopharynx, hypopharynx, esophagus, rectum, larynx, trachea/bronchus/lung, breast, cervix uteri, corpus uteri, prostate, bladder, and brain/central nervous system.

Sites were summed with equal weight. The resulting quantity is selected-site incidence, not modeled radiotherapy demand. It does not incorporate stage, indication, retreatment, fractionation, or site-specific optimal-utilization weights.

## Primary Metrics

For country \(i\), relative selected-site incidence growth was:

```text
growth_i = (cases_2050,i - cases_2024,i) / cases_2024,i
```

Latest-reported megavoltage-unit density was:

```text
density_2050,i = latest-reported MV units_i / cases_2050,i * 1000
```

The descriptive case pressure shown in Table 1 was:

```text
cases_per_unit_i = cases_2050,i / latest-reported MV units_i
```

## Primary Screen

Percentile thresholds were calculated among the 150 countries with matched resource records. A country met both primary thresholds when:

```text
growth_i >= complete-case 75th percentile of growth
and
density_2050,i <= complete-case 25th percentile of density_2050
```

The thresholds were 119.2% growth and 0.506 megavoltage units per 1000 projected 2050 selected-site cases. Twenty-three countries met both thresholds. The quartiles are relative screening rules within the observed country set, not clinical or engineering adequacy standards.

## Burden-Version Stability Analysis

The screen was independently recalibrated within each burden version because the question was whether applying the same prespecified rule to an updated modeled burden would alter country classification. The same latest-reported resource count was joined to both versions.

Among 149 common matched countries:

- data version 2022 identified 25 countries;
- data version 2024 identified 23 countries;
- 22 were retained in both versions;
- Guatemala was newly classified with version 2024;
- Nigeria, Togo, and Zimbabwe were no longer classified with version 2024;
- the Jaccard similarity was 0.846.

Figure 1b plots each country's within-version growth percentile to show both rank stability and threshold crossings without treating the two burden versions as calendar-year observations.

## Denominator and Threshold Sensitivity

The primary growth and density metrics share projected 2050 case counts. Two analyses assessed whether this mathematical coupling drove classification.

First, density was recalculated using baseline selected-site incidence:

```text
density_2024,i = latest-reported MV units_i / cases_2024,i * 1000
```

The 25th-percentile threshold was 0.860 units per 1000 baseline cases. Twenty-one of the 23 primary countries remained classified; Guatemala and the Syrian Arab Republic did not.

Second, a stricter definition combined the 80th percentile of growth with the 20th percentile of primary density. Sixteen of the 23 primary countries remained classified.

Additional prespecified robustness analyses used all-cancer incidence, a 2040 projection horizon, radiotherapy-centre density instead of unit density, mortality growth instead of incidence growth, exclusion of countries with fewer than 1000 baseline selected-site cases, and alternative percentile thresholds. Table 1 reports, for each primary country, retention across seven alternatives: baseline-density denominator, strict 80th/20th percentiles, all-cancer burden, 2040 horizon, centre density, mortality growth, and exclusion of small baseline counts.

## Missing-Resource Identification Bounds

Missing resource status was summarized without imputing equipment counts. For a stratum with \(N\) Global Cancer Observatory country records, \(O\) observed countries meeting both thresholds, and \(H\) high-growth countries without a resource match:

```text
lower bound = O / N
upper bound = (O + H) / N
```

The lower bound treats every resource-unknown country as not meeting the resource criterion. The upper bound treats every high-growth resource-unknown country as meeting it. Resource-unknown countries below the growth threshold cannot meet both criteria and therefore do not widen the interval.

Globally, the observed count was 23 of 186 and 17 high-growth countries were resource-unknown, yielding 12.4%-21.5%. In the African Region, the corresponding values were 18 observed and 10 high-growth resource-unknown countries among 47 records, yielding 38.3%-59.6%. These are deterministic partial-identification bounds, not confidence intervals.

## Construct Diagnostics

Spearman rank correlations described relationships among country-level measures. Selected-site growth correlated strongly with all-cancer growth (rho = 0.997) and population growth (rho = 0.903). Its correlation with unit density was weaker when density used baseline rather than projected cases (rho = -0.378 vs -0.699), supporting the denominator sensitivity analysis. Growth ranks across burden versions had rho = 0.922.

No P values or confidence intervals were calculated. The analyses describe the complete available country records rather than a probability sample, and the inferential target was classification stability rather than a null hypothesis.

## Figure and Table Construction

Figure 1a maps latest-reported units per 1000 projected 2050 selected-site cases. Red outlines identify countries meeting both primary thresholds; hatching identifies resource-unknown countries.

Figure 1b compares within-version growth percentiles among the 149 common matched countries while holding the resource snapshot fixed. Colors identify retained, newly classified, and no-longer-classified countries. Point area represents version-2024 projected 2050 selected-site cases.

Figure 1c presents World Health Organization region lower and upper identification bounds among all country records. The plotted intervals do not represent sampling uncertainty.

Table 1 gives country characteristics, burden-version classification, baseline-density sensitivity, and retention across seven alternative analyses for the 23 countries meeting both primary thresholds.

## Reproduction

From the repository root, run:

```bash
python scripts/data_fetch/fetch_gco_2022_cancer_tomorrow_predictions.py --chunk-size 50
python scripts/analysis/build_radiotherapy_resource_mismatch.py
python scripts/analysis/analyze_radiology_v6_additional_analyses.py
python scripts/figures/make_radiotherapy_resource_mismatch_figure.py
python scripts/tables/build_radiology_revision_tables.py
```

The repository distributes derived analytic data, display source data, code, a sanitized build log, a manifest, and SHA-256 checksums. Raw Global Cancer Observatory and Directory of Radiotherapy Centres source exports are not redistributed.
