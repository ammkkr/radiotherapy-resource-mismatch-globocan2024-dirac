# Radiotherapy Resource Mismatch Dataset QC

Run finished UTC: 2026-07-08T22:14:04.232211+00:00

## Inputs

- GLOBOCAN Cancer Tomorrow: `data\interim\gco_cancer_tomorrow_2024_predictions_long.csv`
- DIRAC country resources: `data\interim\dirac_country_resources.csv`

## Merge Summary

- GCO country records: 186
- DIRAC country records: 155
- Matched records: 150
- GCO records without DIRAC match: 36
- DIRAC records without GCO match: 5

## Mismatch Thresholds

- Complete-case records: 150
- High-growth threshold, relative RT-relevant case growth q75: 1.1921
- Low-resource threshold, MV units per 1000 RT-relevant 2050 cases q25: 0.5063
- High-growth/low-resource countries: 23

## Interpretation Notes

- Countries absent from the DIRAC country table are treated as missing resource observations, not as zero-resource countries.
- The primary mismatch analysis is therefore a complete-case comparison among GLOBOCAN country records with a DIRAC ISO3 match.

## Top 10 Countries by RT Mismatch Score

| Rank | Country | ISO3 | WHO region | Income | RT cases 2050 | Relative growth | MV units | MV units per 1000 RT cases |
|---:|---|---|---|---|---:|---:|---:|---:|
| 1 | Ethiopia | ETH | AFRO | Low income | 99613 | 1.409 | 3 | 0.030 |
| 2 | Congo, Democratic Republic of | COD | AFRO | Low income | 91916 | 1.391 | 1 | 0.011 |
| 3 | Uganda | UGA | AFRO | Low income | 46945 | 2.027 | 5 | 0.107 |
| 4 | Yemen | YEM | EMRO | Low income | 33953 | 1.836 | 1 | 0.029 |
| 5 | Malawi | MWI | AFRO | Low income | 34319 | 1.668 | 2 | 0.058 |
| 6 | Zambia | ZMB | AFRO | Lower middle income | 28900 | 2.028 | 3 | 0.104 |
| 7 | Nigeria | NGA | AFRO | Lower middle income | 202738 | 1.089 | 17 | 0.084 |
| 8 | Angola | AGO | AFRO | Lower middle income | 36804 | 1.436 | 3 | 0.082 |
| 9 | Tanzania, United Republic of | TZA | AFRO | Lower middle income | 44604 | 1.637 | 11 | 0.247 |
| 10 | Mozambique | MOZ | AFRO | Low income | 28686 | 1.462 | 1 | 0.035 |

## Output Files

- country dataset: `data/radiotherapy_resource_mismatch_country_2024_2050.csv`
- top countries: `data/table_09_radiotherapy_mismatch_top_countries.csv`
- region summary: `data/table_10_radiotherapy_mismatch_region_summary.csv`
- site burden: `data/table_11_radiotherapy_relevant_site_burden.csv`
- site rules: `data/radiotherapy_relevant_site_rules.csv`
- sanitized build log: `data/build_log_sanitized.json`
