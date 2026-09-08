# Country data dictionary

One row represents one of the 186 country or territory records in the
GLOBOCAN 2024 dictionary. Blank resource fields mean unmatched or unknown,
not zero. Counts use both sexes and all ages. All computations use unrounded
values; display tables round only at presentation.

| Field | Definition |
| --- | --- |
| country_iso3 | Source-supplied three-letter linkage identifier |
| population_label | Display name; encoding normalised for CIV, CUW, REU and TUR |
| source_population_label | Original archived name, including any source encoding artifacts |
| who_region, income_label, hdi_label | Categories from the archived GLOBOCAN dictionary; blank means unclassified |
| units | Recorded DIRAC HePhotonAndElectronBeamRt count in the fixed July 2026 extract |
| centres | Recorded DIRAC RTCentersWithRT count |
| last_update | Displayed DIRAC LastUpdate year; not a verified country census year |
| cases_YEAR | Sum of 14 selected-site incident case counts, GLOBOCAN 2024 release |
| deaths_YEAR | Sum of 14 selected-site mortality counts, GLOBOCAN 2024 release |
| all40_YEAR | All-cancer incidence excluding non-melanoma skin cancer, code 40 |
| all39_YEAR | All-cancer incidence including non-melanoma skin cancer, code 39; never added to code 40 |
| population_YEAR | Population projection from the code-40 incidence response |
| release2022_cases_YEAR | Sum of the same 14 sites in the 2022 release, for 2025 or 2050; blank if no corresponding country record |
| thirteen_site_cases_YEAR | Selected-site incident count omitting rectum, used only for an explicit sensitivity check |
| growth | cases_2050 / cases_2024 - 1; proportion, not percentage |
| density | 1000 * units / cases_2050 |
| baseline_density | 1000 * units / cases_2024 |
| primary | Growth at or above the matched-country 75th percentile AND density at or below the 25th percentile |
| support_7 | Number of the seven defined alternative specifications selecting that country; 0-7, not a probability |
| additional_units_to_cross_cutoff | Minimum nonnegative integer increase crossing the fixed density threshold; unknown if units unknown; not procurement need |

YEAR is 2024, 2025, 2030, 2040 or 2050 unless the field description specifies
otherwise. These are projected diagnoses or deaths, not unique radiotherapy
patients, courses or fractions. The 2024 values are the release's baseline
estimates; later years assume zero annual change in age-specific cancer rates.

`version_comparison.csv` contains old/new growth and density, flags, eligibility
and growth percentiles using the common 2025 baseline. `sensitivity_results.csv`
reports the denominator, total selected, eligible primary records, retained
primary records and country-code transitions. `missingness_bounds.csv` gives
equal-country proportions conditional on fixed primary cutoffs. `correlation_diagnostics.csv`
reports Spearman correlations with the exact country count for each comparison.
