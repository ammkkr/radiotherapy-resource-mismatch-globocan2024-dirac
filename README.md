# Cancer Tomorrow Projections and Latest-Reported Radiotherapy MV-Unit Density

This repository contains the derived data, figure source data, and analysis code
for a country-level analysis linking IARC Global Cancer Observatory Cancer
Tomorrow projections with IAEA DIRAC radiotherapy-resource data. The API data
version used for the extracted Cancer Tomorrow projections was 2024.

## Study Question

Which countries with matched DIRAC data have both high projected growth in
selected cancer-site incidence by 2050 and lower-quartile latest-reported megavoltage
(MV) unit density?

## Repository Contents

- `data/radiotherapy_resource_mismatch_country_2024_2050.csv`: country-level
  analytic dataset.
- `data/figure_1_country_source_data.csv` and
  `data/figure_1_region_source_data.csv`: source data for Figure 1.
- `data/table_15_radiology_research_letter_country_characteristics.csv`: source
  data for the manuscript Table 1 country-characteristics table.
- `data/table_14_radiology_research_letter_robustness.csv`,
  `data/table_16_radiology_dirac_match_comparison.csv`, and
  `data/table_17_radiology_dirac_missing_scenario.csv`: repository companion
  analyses for sensitivity and DIRAC-missingness checks.
- `data/selected_cancer_site_rules.csv` and
  `data/table_11_selected_site_burden.csv`: selected-site definitions and
  global selected-site burden summaries.
- `data/table_18_radiology_unit_pressure_top_countries.csv`: source data for
  Figure 1b country pressure profiles.
- `figures/`: exported Figure 1 files in SVG, PDF, and PNG formats.
- `scripts/`: analysis and figure-generation scripts.
- `docs/expanded_methods.md`: expanded reproducibility methods for the public
  code repository; this file is not a formal journal supplement.
- `docs/qc_report.md`: merge and threshold quality-control summary.
- `MANIFEST.csv` and `checksums_sha256.txt`: file inventory and checksums.

## Reproduce The Figure

Install dependencies:

```bash
pip install -r requirements.txt
```

Then run:

```bash
python scripts/make_radiotherapy_resource_mismatch_figure.py
```

The figure script uses the included processed country dataset, Natural Earth
country geometry, and sanitized build log.

## Rebuilding From Source

The scripts document the full analytic workflow. Raw GCO Cancer Tomorrow and
DIRAC source exports are not redistributed in this repository. Users who want to rebuild the
country-level dataset from raw inputs should obtain source data from the Global
Cancer Observatory and IAEA DIRAC and place them in the expected project
structure before running `scripts/build_radiotherapy_resource_mismatch.py`.

## Citation

Please cite the manuscript associated with this repository, Global Cancer
Observatory Cancer Tomorrow, IAEA DIRAC, and Natural Earth as applicable.
