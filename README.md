# Radiotherapy Resource Mismatch Using GLOBOCAN 2024 and DIRAC

This repository contains the derived data, figure source data, and analysis code
for a country-level analysis linking GLOBOCAN 2024 cancer projections with IAEA
DIRAC radiotherapy-resource data.

## Study Question

Which countries have both high projected growth in radiotherapy-relevant cancer
burden by 2050 and low current megavoltage radiotherapy-resource density?

## Repository Contents

- `data/radiotherapy_resource_mismatch_country_2024_2050.csv`: country-level
  analytic dataset.
- `data/figure_1_country_source_data.csv` and
  `data/figure_1_region_source_data.csv`: source data for Figure 1.
- `data/table_14_radiology_research_letter_robustness.csv`: source data for
  Table 1 robustness analysis.
- `figures/`: exported Figure 1 files in SVG, PDF, and PNG formats.
- `scripts/`: analysis and figure-generation scripts.
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

The scripts document the full analytic workflow. Raw GLOBOCAN and DIRAC source
exports are not redistributed in this repository. Users who want to rebuild the
country-level dataset from raw inputs should obtain source data from the Global
Cancer Observatory and IAEA DIRAC and place them in the expected project
structure before running `scripts/build_radiotherapy_resource_mismatch.py`.

## Citation

Please cite the manuscript associated with this repository, GLOBOCAN / Global
Cancer Observatory, IAEA DIRAC, and Natural Earth as applicable.
