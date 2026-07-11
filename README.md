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
3. How widely can missing resource observations change global and regional conclusions?

The two-threshold result is an early-warning screen for capacity verification.
It is not a model of radiotherapy indication, optimal utilization, functional
machine capacity, or patient-level unmet need.

## Key Derived Results

- 23 of 150 matched countries met both primary thresholds.
- 22 countries were retained across burden versions; Guatemala was newly
  classified, while Nigeria, Togo, and Zimbabwe were no longer classified.
- A baseline-case density denominator retained 21 of 23 countries, and stricter
  80th-growth/20th-density percentiles retained 16.
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

The figure script uses the included v6 derived dataset, identification-bound
table, Natural Earth geometry, and sanitized logs.

## Rebuild From Source

Raw Cancer Tomorrow and Directory of Radiotherapy Centres source exports are
not redistributed. The fetch and analysis scripts document the complete
workflow. Users rebuilding the analysis should obtain or fetch source data from
the original providers, preserve the documented directory structure, and then
run the primary and additional analysis scripts. The burden-version comparison
must use the same resource snapshot for both versions.

## Citation

Please cite the associated manuscript, Global Cancer Observatory Cancer
Tomorrow, the IAEA Directory of Radiotherapy Centres, and Natural Earth as
applicable.
