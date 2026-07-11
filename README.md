# Cancer Tomorrow Burden Updates and Latest-Reported Radiotherapy Capacity

This repository contains derived country-level data, display source data, and
analysis code for a study linking IARC Global Cancer Observatory Cancer Tomorrow
projections with latest-reported IAEA Directory of Radiotherapy Centres records.
The primary burden analysis used API data version 2024; the stability analysis
repeated the screen with data version 2022 while holding the same resource
snapshot fixed.

## Study Questions

1. Which countries have both upper-quartile projected selected-site incidence
   growth and lower-quartile latest-reported megavoltage-unit density?
2. Does updating the modeled cancer-burden version change the countries identified?
3. Do classifications persist when the burden definition or resource-record
   recency changes?
4. How widely can missing resource observations change global and regional conclusions?

The two-threshold result is an early-warning screen for capacity verification.
It is not a model of radiotherapy indication, optimal utilization, functional
machine capacity, or patient-level unmet need.

## Key Derived Results

- 23 of 150 matched countries met both primary thresholds.
- With a common 2025 baseline and fixed resource snapshot, each burden version
  identified 24 countries; 23 were retained, Papua New Guinea entered, and
  Zimbabwe exited in version 2024.
- The all-cancer definition retained all 23 primary countries.
- A baseline-case density denominator retained 21 of 23 countries, and stricter
  80th-growth/20th-density percentiles retained 16.
- Restriction to 125 countries with DIRAC records dated 2023 or later retained
  all 16 eligible primary countries; 7 primary countries had older records.
- With 17 high-growth countries lacking resource matches, the global
  identification interval was 12.4%-21.5%; the African Region interval was
  38.3%-59.6%. These are missing-data bounds, not confidence intervals.

## Repository Contents

- `data/radiotherapy_resource_mismatch_country_v6_sensitivity.csv`: integrated
  country-level dataset used by the final figure and stability analyses.
- `data/figure_1_country_source_data.csv` and
  `data/figure_1_region_source_data.csv`: Figure 1 source data.
- `data/table_15_radiology_research_letter_country_characteristics.csv`: final
  manuscript Table 1 source.
- `data/table_19_burden_version_reclassification.csv`: version transitions.
- `data/table_20_denominator_threshold_stability.csv`: denominator and
  threshold sensitivity.
- `data/table_21_missing_resource_identification_bounds.csv`: global, region,
  income, and human-development bounds.
- `data/table_22_construct_diagnostics.csv`: rank-correlation diagnostics.
- `data/table_23_country_stability.csv`: seven-analysis support by primary country.
- `data/table_24_common_2025_baseline_burden_version_reclassification.csv`:
  burden-version transitions using the same 2025 baseline.
- `data/table_25_dirac_record_recency_sensitivity.csv` and
  `data/table_26_dirac_record_recency_country_detail.csv`: resource-record
  recency analysis and country detail.
- `data/table_27_country_stability_updated.csv`: updated country-level Table 1 source.
- `data/table_28_stability_summary.csv`: compact policy-facing stability table.
- `figures/`: Figure 1 in SVG, PDF, and PNG formats.
- `docs/expanded_methods.md`: expanded reproducibility methods; this is
  repository documentation, not a formal journal supplement.
- `scripts/`: data-fetch, analysis, and figure-generation code.
- `MANIFEST.csv` and `checksums_sha256.txt`: public-file inventory and checksums.

## Reproduce Figure 1

```bash
pip install -r requirements.txt
python scripts/make_radiotherapy_resource_mismatch_figure.py
```

The figure script uses the included derived dataset, common-baseline version
analysis, identification-bound table, Natural Earth geometry, and sanitized logs.

## Rebuild From Source

Raw Cancer Tomorrow and Directory of Radiotherapy Centres source exports are
not redistributed. The fetch and analysis scripts document the complete
workflow. Users rebuilding the analysis should obtain or fetch source data from
the original providers, preserve the documented directory structure, and then
run the primary and additional analysis scripts. The main burden-version
comparison uses a common 2025 baseline, a 2050 horizon, the same 149 matched
countries, and the same resource snapshot for both versions.

Full rebuilding expects `data/source/gco_cancer_tomorrow_2022_predictions_long.csv`,
`data/source/gco_cancer_tomorrow_2024_predictions_long.csv`, and
`data/source/dirac_country_resources.csv`. The GCO fetch scripts additionally
require the provider population dictionaries and site-rule file documented in
their module constants. Figure 1 can be reproduced directly from the included
derived files without downloading raw source data.

## Citation

Please cite the associated manuscript, Global Cancer Observatory Cancer
Tomorrow, the IAEA Directory of Radiotherapy Centres, and Natural Earth as
applicable.
