# Radiotherapy Resource Mismatch Figure Contract

Core conclusion: Countries with rapidly growing radiotherapy-relevant cancer burden are concentrated in settings with sparse current megavoltage radiotherapy capacity.
Figure archetype: asymmetric mixed-modality figure.
Target journal/output: Research Letter main figure, double-column width, editable SVG/PDF plus high-resolution PNG.
Backend: Python/matplotlib only.
Final size: 7.65 x 6.08 inches before tight bounding-box export.

## Panel Map

- a: Hero world map; fill encodes MV units per 1000 projected 2050 RT-relevant cancer cases, red outline marks countries crossing both high-growth and low-resource thresholds, hatch marks DIRAC-missing observations.
- b: Top-scoring country profile plot; y-axis is mismatch rank, x-axis is MV density, point area encodes projected 2050 RT-relevant cases, right text gives relative growth and current MV units.
- c: WHO-region burden-resource quadrant; x-axis is aggregate relative RT-relevant case growth, y-axis is aggregate MV density, bubble area encodes projected 2050 cases, labels give high-growth/low-resource count over matched countries.

## Evidence Hierarchy

- Hero evidence: geographic co-localisation of projected RT-relevant burden and low current MV capacity.
- Validation evidence: ranked country profiles show the same countries carry low MV density and large projected growth.
- Regional synthesis: AFRO is visually separated by high relative growth, low MV density, and the largest count of high-growth/low-resource countries.

## Thresholds and n

- Complete-case countries: 150.
- High-growth threshold: relative RT-relevant case growth q75 = 1.1921.
- Low-resource threshold: MV units per 1000 projected 2050 RT-relevant cases q25 = 0.5063.
- High-growth/low-resource countries: 23.

## Reviewer-Risk Notes

- DIRAC-absent countries are treated as missing resource observations, not zero-capacity countries.
- The figure visualises radiotherapy-relevant incident cancers, not modelled radiotherapy demand or utilisation.
- MV units are current DIRAC country-table counts and are not projected to 2050.
- Panel c uses complete-case regional aggregation for resource denominators.

## Exported Files

- svg: `figures\figure_1_radiotherapy_resource_mismatch_3panel.svg`
- pdf: `figures\figure_1_radiotherapy_resource_mismatch_3panel.pdf`
- png: `figures\figure_1_radiotherapy_resource_mismatch_3panel.png`
- country source data: `data\figure_1_country_source_data.csv`
- region source data: `data\figure_1_region_source_data.csv`