# Radiotherapy Resource Mismatch Figure Contract

Core conclusion: With a common 2025 baseline and fixed resource snapshot, the GCO version 2024 update retained 23 of 24 countries while missing resource observations widened regional uncertainty.
Figure archetype: asymmetric mixed-modality figure.
Target journal/output: Short Communication main figure, double-column width, editable SVG/PDF plus high-resolution PNG.
Backend: Python/matplotlib only.
Final size: 7.65 x 6.18 inches before tight bounding-box export.

## Panel Map

- a: Hero world map; fill encodes latest-reported MV units per 1000 projected 2050 selected-site cancer cases, red outline marks countries meeting both thresholds, and hatching marks DIRAC-missing observations.
- b: Common-baseline burden-version comparison; axes are within-version 2025-2050 projected growth percentiles for 149 common matched countries, colors identify retained, new, and no-longer-screen-positive countries, and point area encodes version 2024 selected-site cases in 2050.
- c: WHO-region partial-identification bounds; the lower endpoint is the observed count divided by all regional GCO country records, and the upper endpoint additionally treats high-growth DIRAC-missing countries as meeting the resource criterion.

## Evidence Hierarchy

- Hero evidence: geographic co-localisation of projected selected-site incidence growth and sparse latest-reported MV-unit density.
- Update evidence: the burden-version comparison separates the GCO data revision from baseline-year and resource changes by using a common 2025 baseline and holding the same latest-reported DIRAC snapshot fixed.
- Missingness evidence: regional intervals show the range compatible with observed and high-growth resource-unknown country records.

## Thresholds and n

- Complete-case countries: 150.
- High-growth threshold: relative selected-site case growth q75 = 1.1921.
- Lower-quartile unit-density threshold: MV units per 1000 projected 2050 selected-site cases q25 = 0.5063.
- Countries meeting both primary thresholds: 23.
- High-growth countries with unknown resource status: 17.
- Countries retained across GCO burden versions: 23.
- New in version 2024: PNG.
- No longer screen-positive in version 2024: ZWE.

## Reviewer-Risk Notes

- DIRAC-absent countries are treated as missing resource observations without assigned measured resource density.
- The figure visualises selected cancer-site incidence, not modelled radiotherapy demand or utilisation.
- MV units are latest-reported DIRAC country-table counts and are not projected to 2050.
- Panel b is a common-2025-baseline burden-version comparison, not a longitudinal resource analysis; the latest-reported DIRAC snapshot is fixed.
- Panel c intervals are deterministic missing-resource bounds, not confidence intervals.

## Exported Files

- svg: `figures\figure_1_radiotherapy_resource_mismatch_3panel.svg`
- pdf: `figures\figure_1_radiotherapy_resource_mismatch_3panel.pdf`
- png: `figures\figure_1_radiotherapy_resource_mismatch_3panel.png`
- country source data: `data\figure_1_country_source_data.csv`
- region source data: `data\figure_1_region_source_data.csv`