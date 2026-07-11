# Expanded Reproducibility Methods

This file documents reproducibility details for the Radiology Research Letter. It is intended for the public code repository and is not a formal Radiology online supplement.

## Study Design

Retrospective ecological country-level analysis of public aggregate cancer projection and radiotherapy-resource data. No individual-level data were used. IRB approval, informed consent, and HIPAA compliance were not applicable.

## Data Sources

Cancer projections were extracted from the International Agency for Research on Cancer Global Cancer Observatory Cancer Tomorrow API data version 2024. The formal GLOBOCAN 2024 estimates are cited in the manuscript as Global Cancer Statistics 2024.

Radiotherapy-resource variables were extracted from the IAEA Directory of Radiotherapy Centres (DIRAC) country-level records. The primary resource count was the latest-reported number of megavoltage photon/electron therapy units. DIRAC records have heterogeneous update years, so the manuscript consistently labels the resource metric as latest-reported.

Natural Earth country geometry was used for the world map.

## Country Matching

Countries were matched by ISO3 code. The analytic country universe contained 186 GLOBOCAN country records. DIRAC records were matched for 150 countries. Thirty-six country records had no DIRAC resource match and were treated as unknown resource observations without assigned measured resource density.

## Selected-Site Incidence

Selected-site incidence was a prespecified pragmatic incidence proxy, not a model of patient-level radiotherapy indication. The set included 14 cancer sites commonly involving radiotherapy in at least some treatment pathways:

- lip/oral cavity;
- oropharynx;
- nasopharynx;
- hypopharynx;
- esophagus;
- rectum;
- larynx;
- trachea, bronchus, and lung;
- breast;
- cervix uteri;
- corpus uteri;
- prostate;
- bladder;
- brain/central nervous system.

All selected sites were summed with equal weight because the analysis was designed as a simple country-level screen, not as a site-specific optimal radiotherapy utilization model.

## Primary Metrics

Relative selected-site case growth was calculated as:

```text
(projected 2050 selected-site incident cases - projected 2024 selected-site incident cases)
/ projected 2024 selected-site incident cases
```

Latest-reported MV-unit density was calculated as:

```text
latest-reported DIRAC MV photon/electron therapy units
/ projected 2050 selected-site incident cases * 1000
```

Selected-site cases per MV unit, used in Figure 1b and Table 1 sorting, was calculated as:

```text
projected 2050 selected-site incident cases
/ latest-reported DIRAC MV photon/electron therapy units
```

## Primary Screen

Analyses used complete cases with matched DIRAC data for the primary resource-density screen. High growth was defined as at or above the complete-case 75th percentile of relative selected-site case growth. Lower-quartile unit density was defined as at or below the complete-case 25th percentile of latest-reported MV units per 1000 projected 2050 selected-site cases.

Countries meeting both criteria were classified as acceleration-risk countries.

Primary thresholds:

- high-growth threshold: 119.2% relative selected-site case growth;
- lower-quartile unit-density threshold: 0.506 MV units per 1000 projected 2050 selected-site cases;
- complete-case countries: 150;
- acceleration-risk countries: 23.

## Missingness Checks

DIRAC-missing countries were not assigned an MV-unit value. Missingness was summarized by WHO region, income group, and HDI group. Among 36 countries without a DIRAC match, 17 exceeded the high-growth threshold and therefore had unknown resource status in additional high-growth settings.

## Sensitivity Analyses

Repository companion sensitivity analyses varied:

- threshold definition, including tertile thresholds;
- burden definition, replacing selected-site incidence with all-cancer incidence;
- projection horizon, using 2040 rather than 2050;
- resource metric, using radiotherapy centers rather than MV units;
- mortality growth rather than incidence growth;
- exclusion of countries with fewer than 1000 baseline selected-site cases.

All 23 primary countries were retained when all-cancer incidence replaced the selected-site definition. Other sensitivity definitions retained 22 to 23 of the 23 primary countries.

## Figure Construction

Figure 1a maps latest-reported MV units per 1000 projected 2050 selected-site cases. Red outlines mark acceleration-risk countries. Hatching marks countries without DIRAC resource matches.

Figure 1b ranks countries by projected 2050 selected-site cases per latest-reported MV unit. Red points meet both acceleration-risk criteria; gray points are other high-pressure countries by case-per-unit ranking.

Figure 1c reports acceleration-risk countries as a percentage of matched DIRAC countries within each WHO region. Labels show acceleration-risk countries over matched countries, and point area encodes matched country counts.

## Reproduction

The main analysis script is:

```bash
python scripts/build_radiotherapy_resource_mismatch.py
```

The figure script is:

```bash
python scripts/make_radiotherapy_resource_mismatch_figure.py
```

The repository includes processed analytic data, figure source data, table source data, a sanitized build log, a manifest, and SHA-256 checksums.
