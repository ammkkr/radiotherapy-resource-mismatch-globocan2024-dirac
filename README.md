# Inventory Gaps and Radiotherapy Resource Comparisons

**Correction, 9 September 2026:** Use the [corrected analysis folder](radiotherapy_oncology_20260909/README.md)
for the current manuscript. Earlier files omitted rectal cancer from a stated
14-site set and combined overlapping all-cancer aggregate codes. The correction
documents the effect, completes the site set and supplies independently tested
reproduction code. Historical files below are retained for traceability and
must not be treated as the corrected results.

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
  80th-growth/20th-density percentiles retained 17.
- Restriction to 125 countries with DIRAC records dated 2023 or later retained
  all 16 eligible primary countries; 7 primary countries had older records.
- Sixteen primary countries were retained under all seven alternative
  specifications, and 21 were retained under at least six.
- With 18 high-growth countries lacking resource matches, the global
  fixed-cutoff interval was 12.4%-22.0%; the African Region interval was
  38.3%-59.6%. These bounds cover unmatched records only, holding observed
  counts and cutoffs fixed. They are not confidence intervals or access estimates.

Inventory completeness is not treated as a neutral quality issue. A
complete-case-only screen can reward stronger reporting systems, so
resource-unknown countries remain eligible for capacity verification.

## Historical Repository Contents

- `data/radiotherapy_resource_mismatch_country_v6_sensitivity.csv`: integrated
  country-level dataset used by the final figure and stability analyses.
- `data/figure_1_country_source_data.csv` and
  `data/figure_1_region_source_data.csv`: Figure 1 source data.
- `data/table_15_radiology_research_letter_country_characteristics.csv`: legacy
  country-characteristics table retained for reproducibility.
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
- `data/table_29_radiotherapy_oncology_stability.csv`: compact main-table source
  for the Radiotherapy and Oncology Short Communication.
- `data/table_30_radiotherapy_oncology_primary_summary.csv`: primary-group
  burden, unit-pressure, support, and record-recency summary.
- `data/radiotherapy_oncology_short_communication_summary.json`: machine-readable
  summary generated with the journal-specific analysis script.
- `figures/`: Figure 1 in SVG, PDF, and PNG formats.
- `docs/expanded_methods.md`: expanded reproducibility methods; this is
  repository documentation, not a formal journal supplement.
- `scripts/`: data-fetch, analysis, and figure-generation code.
- `MANIFEST.csv` and `checksums_sha256.txt`: public-file inventory and checksums.

Checksums describe repository blob bytes, with LF newlines for text files.
Windows checkouts may convert newlines to CRLF; compare the canonical bytes from
`git show HEAD:path/to/file` when verifying such a checkout. Generated TIFF files,
Python caches and local reproduction outputs are not part of the public inventory.

## Reproduce the Corrected Analysis

```bash
pip install -r requirements.txt
python radiotherapy_oncology_20260909/reproduce.py --data radiotherapy_oncology_20260909/data --output reproduced_ro_20260909
python radiotherapy_oncology_20260909/make_figure.py
python scripts/test_burden_integrity.py
```

## Historical Figure Reproduction

The commands below reproduce the superseded outputs, not the current manuscript.

```bash
pip install -r requirements.txt
python scripts/make_radiotherapy_resource_mismatch_figure.py
python scripts/analyze_radiotherapy_oncology_short_communication.py
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
