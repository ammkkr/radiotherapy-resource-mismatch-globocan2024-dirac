# Radiotherapy Resource Mismatch Figure Contract

Core conclusion: Countries with upper-quartile projected selected cancer-site incidence growth are concentrated in settings with sparse latest-reported megavoltage unit density.
Figure archetype: asymmetric mixed-modality figure.
Target journal/output: Research Letter main figure, double-column width, editable SVG/PDF plus high-resolution PNG.
Backend: Python/matplotlib only.
Final size: 7.65 x 6.08 inches before tight bounding-box export.

## Panel Map

- a: Hero world map; fill encodes latest-reported MV units per 1000 projected 2050 selected-site cancer cases, red outline marks countries crossing both high-growth and lower-quartile unit-density thresholds, hatch marks DIRAC-missing observations.
- b: Country pressure profile plot; y-axis is selected-site cases per MV-unit pressure rank, x-axis is MV density, point area encodes projected 2050 selected-site cases, right text gives relative growth and latest-reported MV units.
- c: WHO-region acceleration-risk proportion plot; x-axis is acceleration-risk countries as a share of matched countries, labels give acceleration-risk count over matched countries, and point area encodes the number of matched countries.

## Evidence Hierarchy

- Hero evidence: geographic co-localisation of projected selected-site incidence growth and sparse latest-reported MV-unit density.
- Validation evidence: ranked country profiles show countries with the highest selected-site cases per latest-reported MV unit.
- Regional synthesis: WHO-region proportions show where acceleration-risk countries are concentrated among matched DIRAC records.

## Thresholds and n

- Complete-case countries: 150.
- High-growth threshold: relative selected-site case growth q75 = 1.1921.
- Lower-quartile unit-density threshold: MV units per 1000 projected 2050 selected-site cases q25 = 0.5063.
- Acceleration-risk countries: 23.

## Reviewer-Risk Notes

- DIRAC-absent countries are treated as missing resource observations without assigned measured resource density.
- The figure visualises selected cancer-site incidence, not modelled radiotherapy demand or utilisation.
- MV units are latest-reported DIRAC country-table counts and are not projected to 2050.
- Panel c reports regional proportions of country-level classifications, avoiding classification of regional aggregate ratios.

## Exported Files

- svg: `figures\figure_1_radiotherapy_resource_mismatch_3panel.svg`
- pdf: `figures\figure_1_radiotherapy_resource_mismatch_3panel.pdf`
- png: `figures\figure_1_radiotherapy_resource_mismatch_3panel.png`
- country source data: `data\figure_1_country_source_data.csv`
- region source data: `data\figure_1_region_source_data.csv`