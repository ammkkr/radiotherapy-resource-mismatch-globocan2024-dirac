# Radiotherapy Resource Mismatch Dataset QC

Run finished UTC: 2026-07-11T14:16:05.621868+00:00

## Inputs

- GCO Cancer Tomorrow API data version 2024: `data\interim\gco_cancer_tomorrow_2024_predictions_long.csv`
- DIRAC country resources: `data\interim\dirac_country_resources.csv`

## Merge Summary

- GCO country records: 186
- DIRAC country records: 155
- Matched records: 150
- GCO records without DIRAC match: 36
- DIRAC records without GCO match: 5

## Two-Threshold Screen

- Complete-case records: 150
- High-growth threshold, relative selected-site case growth q75: 1.1921
- Lower-quartile unit-density threshold, MV units per 1000 selected-site 2050 cases q25: 0.5063
- Countries meeting both thresholds: 23

## Interpretation Notes

- Countries absent from the DIRAC country table are treated as missing resource observations without assigned resource-density values.
- The primary mismatch analysis is therefore a complete-case comparison among GLOBOCAN country records with a DIRAC ISO3 match.

## Top 10 Countries by Selected-Site Cases per MV Unit

| Rank | Country | ISO3 | WHO region | Income | Selected-site cases 2050 | Relative growth | MV units | Selected-site cases per MV unit |
|---:|---|---|---|---|---:|---:|---:|---:|
| 1 | Congo, Democratic Republic of | COD | AFRO | Low income | 91916 | 1.391 | 1 | 91916 |
| 2 | Yemen | YEM | EMRO | Low income | 33953 | 1.836 | 1 | 33953 |
| 3 | Ethiopia | ETH | AFRO | Low income | 99613 | 1.409 | 3 | 33204 |
| 4 | Zimbabwe | ZWE | AFRO | Lower middle income | 30640 | 1.053 | 1 | 30640 |
| 5 | Mozambique | MOZ | AFRO | Low income | 28686 | 1.462 | 1 | 28686 |
| 6 | Korea, Democratic People's Republic of | PRK | SEARO | Low income | 44118 | 0.341 | 2 | 22059 |
| 7 | Niger | NER | AFRO | Low income | 19785 | 1.670 | 1 | 19785 |
| 8 | Malawi | MWI | AFRO | Low income | 34319 | 1.668 | 2 | 17160 |
| 9 | Angola | AGO | AFRO | Lower middle income | 36804 | 1.436 | 3 | 12268 |
| 10 | Nigeria | NGA | AFRO | Lower middle income | 202738 | 1.089 | 17 | 11926 |

## Output Files

- country dataset: `data/radiotherapy_resource_mismatch_country_2024_2050.csv`
- top countries: `data/table_18_radiology_unit_pressure_top_countries.csv`
- region summary: `data/table_10_region_summary.csv`
- site burden: `data/table_11_selected_site_burden.csv`
- site rules: `data/selected_cancer_site_rules.csv`
- sanitized build log: `data/build_log_sanitized.json`