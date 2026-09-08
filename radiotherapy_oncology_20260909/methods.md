# Supplementary Material

Title: Inventory gaps and stability of radiotherapy resource comparisons across 150 countries

## Supplementary methods

### Data sources and extraction

We analysed public aggregate records from the IARC Global Cancer Observatory Cancer Tomorrow API (GLOBOCAN releases 2022 and 2024) and IAEA Directory of Radiotherapy Centres (DIRAC). The original Cancer Tomorrow files and DIRAC country extract were archived in July 2026. The DIRAC archive was processed on 8 July 2026 UTC; this timestamp records local processing, not independent confirmation of the original download time. Resource counts were kept fixed throughout the analysis. They should be described as counts in this extract rather than current September 2026 inventories.

Cancer Tomorrow requests used the release-specific endpoint https://gco.iarc.who.int/gateway_prod/api/globocan/v3/{release}/data/prediction/{type}/0/{population_codes}/{cancer_codes}/?prediction_annual=1&ages_group=0_17. Type 0 denotes incidence, type 1 mortality, and sex 0 both sexes. The requests selected all ages and zero annual percentage change in cancer rates. Responses supplied estimates for the release baseline and projection years 2025, 2030, 2035, 2040, 2045, and 2050. The case-count field was cases_pred. The analysis checked uniqueness by release, country identifier, cancer code, measure, sex, and year and compared every archived both-sex API case count with the corresponding long-format record.

Rectal-cancer records (code 9) were retrieved on 9 September 2026 to complete the defined 14-site set in both releases. The supplementary requests simultaneously retrieved lung-cancer records (code 15) as an extraction check. All overlapping lung counts were identical to the July archive. This check supports consistency for the anchor site; it does not establish that every other API field remained unchanged. Supplemental JSON files, query identifiers, retrieval timestamps, and hashes were archived. Every country-measure-year cell in the completed selected-site dataset contained all 14 sites, including observed zero counts. No absent cancer-site record was silently converted to zero.

DIRAC country records were obtained from https://dirac.iaea.org/api/DataGridWebApi/GetCountriesAndRegions. The archive contained 155 distinct country codes, of which 150 matched the 186 GLOBOCAN 2024 country or territory records. Matching used the supplied three-letter codes without fuzzy country-name matching. The numerator was HePhotonAndElectronBeamRt; RTCentersWithRT supplied the centre sensitivity measure. LastUpdate was retained as a displayed country-level update year. Its meaning is not a synchronised census date for all facilities. None of the 150 matched records had an observed zero megavoltage count. The 36 unmatched records were retained as unknown, not assigned zero units or assumed to have no radiotherapy services.

Region, income, and HDI categories were copied from the archived GLOBOCAN population dictionary and were not relabelled to a later World Bank fiscal-year classification. One record had no supplied WHO-region label, five had no income label, and ten had no HDI label; these remain explicit unclassified strata. The regional figure excludes the unclassified record, while the global calculation includes it. The study contains no individual health data, recruitment, intervention, or patient follow-up.

### Selected-site case volume

The selected-site measure is the unweighted sum of lip and oral cavity (code 1; C00-06), oropharynx (3; C09-10), nasopharynx (4; C11), hypopharynx (5; C12-13), oesophagus (6; C15), rectum (9; C19-20), larynx (14; C32), trachea/bronchus/lung (15; C33-34), breast (20; C50), cervix uteri (23; C53), corpus uteri (24; C54), prostate (27; C61), bladder (30; C67), and brain/central nervous system (31; C70-72). The set contains non-overlapping site categories with established radiotherapy applications. It is not exhaustive and does not assign clinical radiotherapy utilisation weights. The site selection and all robustness analyses should be considered exploratory, not preregistered confirmatory analyses.

All-cancer sensitivity analyses used aggregate code 40, all cancers excluding non-melanoma skin cancer. They did not sum site-specific codes, add colorectal subcategories to their parent category, or combine overlapping aggregate codes 39 and 40. The country sums for 2024 were 19,485,880 for code 40 and 20,623,564 for code 39. The 14-site country sums were 10,846,289 in 2024 and 16,396,240 in 2050. Selected-site cases are incident diagnoses, not unique patients requiring radiotherapy, treatment courses, fractions, or observed unmet need.

### Primary classification and mathematical coupling

Let B0 and B1 denote selected-site incidence in 2024 and 2050, and M the recorded megavoltage count. Growth equals B1/B0 - 1, future-case density equals 1000M/B1, and baseline-case density equals 1000M/B0. Therefore future-case density equals baseline-case density divided by one plus growth. This identity creates mathematical coupling between the primary indicators. The baseline-density sensitivity removes the shared future-case term, but it does not make the indicators statistically independent because they still contain baseline burden and may share structural determinants.

Primary thresholds were the 75th growth percentile and 25th density percentile across 150 matched countries, using pandas linear interpolation. A country was selected when growth was greater than or equal to 1.1753085502 and density was less than or equal to 0.4710147743 units per 1000 projected cases. The resource numerator was held at its recorded value. These thresholds quantify relative positions in the observed reference set and are not adequacy standards. Country comparisons give equal weight to each country; neither country-selection proportions nor their bounds are patient-weighted.

### Comparisons between burden releases

Release comparisons used 149 countries with records in both burden releases and the fixed DIRAC extract. The 2025-to-2050 growth interval was identical for both releases. Within-release cutoffs were 112.7483% growth and 0.483129 units per 1000 cases for the 2022 release, and 110.4897% and 0.468450 for the 2024 release. Both selected 24 countries, with 23 in common. The Jaccard index equals intersection size divided by union size: 23/25 = 0.92. This comparison differs from the primary analysis in its country set and baseline year; its retention denominator is therefore 24, not 23.

A fixed-threshold analysis applied the 2022-release cutoffs to 2024-release values in the same common set. It selected 23 countries and retained 22 of the 24 selected using 2022-release values. This separates changes in the indicator values from changes in the relative cutoffs. It does not isolate a causal effect of a database update or provide independent validation, because releases share data sources and projection assumptions. The growth correlation reported in the main text uses this common 2025-to-2050 interval and these 149 matched countries.

### Alternative specifications

Seven alternatives were evaluated separately. The all-cancer alternative used code 40 incidence for both growth and density. The baseline-density alternative retained 2024-to-2050 selected-site growth but used 2024 cases in the density denominator. Stricter thresholds were the 80th growth percentile and 20th density percentile. The 2040 alternative changed both the growth horizon and density denominator to 2040. The centre alternative replaced megavoltage counts with centres offering radiotherapy, scaled per 10,000 future cases. The mortality alternative used selected-site mortality growth with the primary incidence-based density measure. The small-burden restriction recalculated thresholds after excluding countries with fewer than 1000 selected-site incident cases in 2024. Support counts record selection across these seven related alternatives; they are not probabilities of true priority status.

Additional analyses were not included in the seven-specification count. Omitting rectal cancer reproduced a 13-site extraction and retained the same 23 primary countries. Record-recency analyses required LastUpdate of 2023, 2024, 2025, or 2026 or later and recalculated thresholds within each eligible set. Excluded primary countries were labelled ineligible rather than not selected. These restrictions change both the observed sample and the reference distribution; apparent persistence among eligible countries does not validate excluded records.

Deterministic inventory perturbations added one, two, or five units to every matched country, or multiplied all recorded unit counts by 0.8, 1.2, or 1.5. Burden estimates and both primary thresholds were held fixed. Recalculating a quantile after multiplying every numerator by a common factor would preserve ranks mechanically, so it was not used for these checks. Perturbation magnitudes were illustrative and were not calibrated from measurement-error data. They do not simulate funded expansion, estimate functional capacity, or define treatment adequacy. Country data also record the smallest integer increase that would cross the fixed low-density threshold; this is an arithmetic classification boundary, not the number of machines a country needs.

### Bounds for unmatched resource records

For each stratum, let N be all GLOBOCAN country records, O the observed countries meeting both primary thresholds, and H the unmatched countries above the fixed growth threshold. The lower proportion is O/N and the upper proportion is (O+H)/N. Unmatched countries below the growth threshold cannot meet both criteria. The growth and density cutoffs remain those estimated from the primary complete-case set; they are not recalculated after assigning hypothetical resource values. Thus the bounds refer to membership at fixed cutoffs, conditional on recorded values for matched countries, rather than the quartiles of a fully observed global resource distribution.

The global bounds are 23/186 to 41/186 (12.4%-22.0%), and African Region bounds are 18/47 to 28/47 (38.3%-59.6%). The bounds account only for unmatched inventories. They do not cover errors in observed unit counts, operation status, cancer estimates, demographic projections, or population classifications. They are deterministic extreme scenarios, not confidence intervals or evidence that all unmatched high-growth countries have low resource density.

### Statistical and computational reporting

Descriptive calculations included counts, proportions, linear quantiles, Spearman rank correlations, and Jaccard overlap. Sampling-based hypothesis tests were not used because the target was the available set of country records. This choice does not eliminate model uncertainty or measurement error. No causal interpretation, prospective accuracy claim, or patient-level inference was made. The source archive and derived data were checked for unique joins, full site coverage, count consistency, fixed resource numerators, and the density-growth identity.

Source identifiers, query parameters, hashes, software versions, derived country records, and the complete sensitivity table are available in the repository at https://github.com/ammkkr/radiotherapy-resource-mismatch-globocan2024-dirac. The figure uses an archived Natural Earth country map. Panel b point areas are a linear function of log10 projected 2050 cases; colours denote the joint classification in both releases. The dates of acquisition and processing are retained, and later retrievals are identified separately. Raw source exports are not redistributed.

## Supplementary Table S1

Supplementary Table S1. Country characteristics and stability among the 23 countries selected in the primary analysis

Note: Cases refer to the 14-site incidence measure. Update year is the displayed DIRAC LastUpdate field. Support counts refer to the seven alternatives defined above. MV, megavoltage; LIC, low income; LMIC, lower middle income; UMIC, upper middle income. Classifications follow the source dictionary. Display names were standardised where needed to correct encoding; source names and three-letter identifiers are retained in the country data.

## Supplementary Table S2

Supplementary Table S2. Country-record coverage and conditional bounds under unmatched DIRAC records

Note: All proportions give equal weight to country records. Lower and upper bounds use fixed primary cutoffs and hold observed resource counts fixed. Unclassified categories are retained explicitly. These are not confidence intervals or patient-level access estimates.

## Supplementary Table S3

Supplementary Table S3. Complete specification and inventory-perturbation results

Note: Unless specified otherwise, selection is compared with the 23 primary countries. Update-year restrictions exclude ineligible records and recalculate cutoffs. Inventory perturbations hold primary cutoffs fixed and use uncalibrated illustrative changes. Full country transition lists are supplied as machine-readable data. NMSC, non-melanoma skin cancer.

## Supplementary Table S4

Supplementary Table S4. Rank correlations and shared-denominator diagnostics

Note: Spearman correlations are descriptive. The common-release correlations use 149 matched countries and the same 2025-to-2050 interval. The 2050 and baseline density measures share underlying data and should not be interpreted as independent corroborating tests.

## Supplementary Table S5

Supplementary Table S5. Country records identified for inventory follow-up

Note: This table lists the seven primary countries with pre-2023 update years and the 18 unmatched high-growth countries. A place on this list indicates an information gap under the study definitions, not confirmed absence of radiotherapy, clinical treatment priority, or procurement need. None of the existing records has been independently validated by the authors against a national facility census.
