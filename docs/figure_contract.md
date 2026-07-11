# Radiotherapy Resource Mismatch Figure Contract

Core conclusion: Countries with rapidly growing selected cancer-site incidence are concentrated in settings with sparse current megavoltage unit density.
Figure archetype: asymmetric mixed-modality figure.
Target journal/output: Research Letter main figure, double-column width, editable SVG/PDF plus high-resolution PNG.
Backend: Python/matplotlib only.
Final size: 7.65 x 6.08 inches before tight bounding-box export.

## Panel Map

- a: Hero world map; fill encodes MV units per 1000 projected 2050 selected-site cancer cases, red outline marks countries crossing both high-growth and lower-quartile unit-density thresholds, hatch marks DIRAC-missing observations.
- b: Country pressure profile plot; y-axis is selected-site cases per MV-unit pressure rank, x-axis is MV density, point area encodes projected 2050 selected-site cases, right text gives relative growth and current MV units.
- c: WHO-region burden-resource quadrant; x-axis is aggregate relative selected-site case growth, y-axis is aggregate MV density, bubble area encodes projected 2050 selected-site cases, labels give acceleration-risk count over matched countries.

## Evidence Hierarchy

- Hero evidence: geographic co-localisation of projected selected-site incidence growth and low current MV-unit density.
- Validation evidence: ranked country profiles show countries with the highest selected-site cases per current MV unit.
- Regional synthesis: AFRO is visually separated by high relative growth, low MV density, and the largest count of acceleration-risk countries.

## Thresholds and n

- Complete-case countries: 150.
- High-growth threshold: relative selected-site case growth q75 = 1.1921.
- Lower-quartile unit-density threshold: MV units per 1000 projected 2050 selected-site cases q25 = 0.5063.
- Acceleration-risk countries: 23.

## Reviewer-Risk Notes

- DIRAC-absent countries are treated as missing resource observations, not zero-capacity countries.
- The figure visualises selected cancer-site incidence, not modelled radiotherapy demand or utilisation.
- MV units are current DIRAC country-table counts and are not projected to 2050.
- Panel c uses complete-case regional aggregation for resource denominators.

## Exported Files

- svg: `figures\figure_1_radiotherapy_resource_mismatch_3panel.svg`
- pdf: `figures\figure_1_radiotherapy_resource_mismatch_3panel.pdf`
- png: `figures\figure_1_radiotherapy_resource_mismatch_3panel.png`
- country source data: `data\figure_1_country_source_data.csv`
- region source data: `data\figure_1_region_source_data.csv`