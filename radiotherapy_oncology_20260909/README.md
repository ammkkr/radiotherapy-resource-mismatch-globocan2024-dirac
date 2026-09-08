# Corrected country resource comparison

This folder contains the 9 September 2026 analysis supporting the Short
Communication on inventory gaps and stability of radiotherapy resource comparisons.
It supersedes the earlier 13-site calculations elsewhere in this repository.

## Corrections

The original general-site extraction did not include rectal cancer (code 9),
although the stated site set contained 14 sites. Rectal incidence and mortality
were retrieved for both burden releases on 9 September 2026. Lung counts were
retrieved in the same requests and matched the July archive. Complete site
coverage is now checked rather than treating an absent site as zero.

The earlier all-cancer fallback added overlapping aggregate codes 39 and 40.
The corrected all-cancer sensitivity uses code 40 only (excluding non-melanoma
skin cancer). In 2024 its country sum is 19,485,880, not 40,109,444. Historical
files remain available for traceability, but are not the current analysis.

The revised 14-site case totals are 10,846,289 in 2024 and 16,396,240 in 2050.
The same 23 primary countries remain selected. Their 2050 case total is 853,352,
with 103 recorded units. Stricter thresholds retain 17 rather than 16 countries;
21 rather than 20 are selected in at least six alternatives. Eighteen unmatched
countries exceed the growth threshold; the global fixed-cutoff bounds are now
12.4%-22.0%. The common-2025-baseline growth correlation is 0.933; the earlier
0.922 referred to native baselines and should not be cited for that comparison.

## Reproduction

From the repository root, after installing `requirements.txt`:

```bash
python radiotherapy_oncology_20260909/reproduce.py --data radiotherapy_oncology_20260909/data --output reproduced_ro_20260909
python radiotherapy_oncology_20260909/make_figure.py
python scripts/test_burden_integrity.py
```

The first command independently recomputes classifications, all sensitivity
rows, fixed-threshold release transitions, bounds and correlations from the
shared country aggregates, and checks them against the supplied outputs.
It does not download data, infer absent inventories, or substitute fabricated
case counts. The figure command uses the included Natural Earth geometry.

`data/table1.csv` and `data/table_s1.csv` through `data/table_s5.csv` are the
display data for the main and supplementary tables. `methods.md` documents the
formulas, inclusion rules and interpretation. `data_dictionary.md` defines
country-data fields. `source_checks.json` records source-level checks and the
calculation environment. `metadata.json` records analysis and extraction dates.

## Source-level audit

Original source JSON exports are not redistributed. `source_audit.py` is the
source-level recomputation code; `source_checksums.json` lists the archive paths
and hashes used in the analysis. To rerun it, obtain the source records under
the providers' terms and reproduce the paths listed in that file inside an
archive directory. Set `RADIOTHERAPY_SOURCE_ARCHIVE` to that directory and run
`python radiotherapy_oncology_20260909/source_audit.py`.

The audit also compares against the historical derived file at
`data/processed/radiotherapy_resource_mismatch_country_v6_sensitivity.csv` within
that archive. The same historical file is already shared in this repository's
top-level `data` directory. This comparison identifies the old errors; it is not
used to construct the corrected 14-site burden.

`supplemental_queries.json` contains the actual rectum/lung request URLs and
retrieval time. `fetch_rectum.py` documents those supplementary requests and uses
the same archive environment variable. It expects the country identifiers in
the archived long-format source files. A later API request may differ from the
archived data; an unchanged response is not assumed.

## Interpretation

These are relative country-record classifications, not radiotherapy demand,
machine adequacy standards, functional capacity, or patient-level unmet need.
The resource numerator is the fixed July 2026 DIRAC extract. Update years are
not verified national census dates. Bounds hold observed values and thresholds
fixed and cover unmatched records only. The analyses are exploratory and their
agreement is not independent validation. No national facility census or actual
treatment outcomes were independently collected for this study.
