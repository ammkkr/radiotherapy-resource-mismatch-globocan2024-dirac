from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "data" / "radiotherapy_resource_mismatch_country_v6_sensitivity.csv"
ADDITIONAL_LOG = ROOT / "data" / "additional_analyses_log.json"
TABLE_DIR = ROOT / "data"
REPORT_DIR = ROOT / "data"

STABILITY_TABLE = TABLE_DIR / "table_29_radiotherapy_oncology_stability.csv"
PRIMARY_SUMMARY_TABLE = TABLE_DIR / "table_30_radiotherapy_oncology_primary_summary.csv"
REPORT = REPORT_DIR / "radiotherapy_oncology_short_communication_summary.json"


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"1", "1.0", "true"})


def rounded_quantiles(series: pd.Series) -> tuple[int, int]:
    values = series.quantile([0.25, 0.75])
    return int(round(values.iloc[0])), int(round(values.iloc[1]))


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(DATASET)
    with ADDITIONAL_LOG.open("r", encoding="utf-8") as handle:
        log = json.load(handle)

    matched = as_bool(data["dirac_matched"])
    primary = as_bool(data["current_screen_positive"])
    recent = as_bool(data["recent_dirac_record"])
    primary_data = data.loc[primary].copy()

    if len(data) != 186 or int(matched.sum()) != 150 or int(primary.sum()) != 23:
        raise ValueError("Unexpected analytic universe; rebuild the upstream analysis first.")

    support = primary_data["sensitivity_support_count_7"].astype(int)
    cases_per_unit = primary_data["selected_site_cases_per_mv_unit_2050"].astype(float)
    growth_percent = 100 * primary_data["selected_site_relative_case_growth_2050"].astype(float)
    cases_iqr = rounded_quantiles(cases_per_unit)

    primary_metrics = {
        "countries": int(primary.sum()),
        "latest_reported_mv_units": int(round(primary_data["mv_therapy_units"].sum())),
        "selected_site_cases_2024": int(round(primary_data["selected_site_cases_2024"].sum())),
        "selected_site_cases_2050": int(round(primary_data["selected_site_cases_2050"].sum())),
        "additional_selected_site_cases": int(
            round(primary_data["selected_site_case_increase_2050"].sum())
        ),
        "median_growth_percent": round(float(growth_percent.median()), 1),
        "growth_iqr_percent": [
            round(float(value), 1) for value in growth_percent.quantile([0.25, 0.75])
        ],
        "median_2050_cases_per_mv_unit": int(round(cases_per_unit.median())),
        "cases_per_mv_unit_iqr": list(cases_iqr),
        "cases_per_mv_unit_range": [
            int(round(cases_per_unit.min())),
            int(round(cases_per_unit.max())),
        ],
        "support_all_7": int((support == 7).sum()),
        "support_at_least_6_of_7": int((support >= 6).sum()),
        "support_distribution": {
            str(index): int(value) for index, value in support.value_counts().sort_index().items()
        },
        "matched_records_before_2023": int((matched & ~recent).sum()),
        "primary_records_before_2023": int((primary & ~recent).sum()),
    }

    if primary_metrics["additional_selected_site_cases"] != 492196:
        raise ValueError("Primary-group case increase no longer matches the verified claim ledger.")
    if primary_metrics["support_all_7"] != 16 or primary_metrics["support_at_least_6_of_7"] != 20:
        raise ValueError("Cross-specification support counts changed unexpectedly.")

    primary_summary = pd.DataFrame(
        [
            ["Countries meeting both primary thresholds, No.", primary_metrics["countries"]],
            ["Latest-reported MV units, No.", primary_metrics["latest_reported_mv_units"]],
            ["Projected selected-site cases in 2024, No.", primary_metrics["selected_site_cases_2024"]],
            ["Projected selected-site cases in 2050, No.", primary_metrics["selected_site_cases_2050"]],
            ["Additional projected selected-site cases, No.", primary_metrics["additional_selected_site_cases"]],
            ["Projected growth, median (IQR), %", f"{primary_metrics['median_growth_percent']:.1f} ({primary_metrics['growth_iqr_percent'][0]:.1f}-{primary_metrics['growth_iqr_percent'][1]:.1f})"],
            ["Projected 2050 cases per recorded MV unit, median (IQR)", f"{primary_metrics['median_2050_cases_per_mv_unit']:,} ({primary_metrics['cases_per_mv_unit_iqr'][0]:,}-{primary_metrics['cases_per_mv_unit_iqr'][1]:,})"],
            ["Supported by all 7 alternative specifications, No.", primary_metrics["support_all_7"]],
            ["Supported by at least 6 of 7 alternative specifications, No.", primary_metrics["support_at_least_6_of_7"]],
            ["Primary countries with DIRAC records before 2023, No.", primary_metrics["primary_records_before_2023"]],
        ],
        columns=["Metric", "Value"],
    )
    primary_summary.to_csv(PRIMARY_SUMMARY_TABLE, index=False, encoding="utf-8-sig")

    version = log["version"]
    recency = log["recency"]
    all_cancer = log["all_cancer"]
    decoupled = log["decoupled"]

    stability_rows = [
        {
            "Analysis": "Primary selected-site screen",
            "Countries analyzed, No.": "150",
            "Definition": "Q75 growth; Q25 MV-unit density using 2050 cases",
            "Countries meeting both, No.": "23",
            "Primary countries retained": "Reference",
            "Classification changes or uncertainty": "None",
        },
        {
            "Analysis": "Burden versions, common 2025 baseline",
            "Countries analyzed, No.": "149",
            "Definition": "Version-specific Q75/Q25 thresholds; fixed DIRAC snapshot",
            "Countries meeting both, No.": "24 / 24",
            "Primary countries retained": "23/24 version-2022 countries",
            "Classification changes or uncertainty": f"New: {', '.join(version['new_in_2024'])}; no longer: {', '.join(version['no_longer_in_2024'])}; Jaccard {version['jaccard']:.2f}",
        },
        {
            "Analysis": "All-cancer burden",
            "Countries analyzed, No.": "150",
            "Definition": "All-cancer growth and 2050 all-cancer density",
            "Countries meeting both, No.": str(all_cancer["screen_positive_count"]),
            "Primary countries retained": f"{all_cancer['primary_retained_count']}/23",
            "Classification changes or uncertainty": f"Added: {', '.join(all_cancer['additional_iso3'])}",
        },
        {
            "Analysis": "2024-case density denominator",
            "Countries analyzed, No.": "150",
            "Definition": "Primary growth; Q25 MV-unit density using 2024 cases",
            "Countries meeting both, No.": str(decoupled["baseline_density_screen_positive_count"]),
            "Primary countries retained": f"{decoupled['baseline_density_primary_retained']}/23",
            "Classification changes or uncertainty": f"Not retained: {', '.join(decoupled['baseline_density_not_retained'])}",
        },
        {
            "Analysis": "Stricter relative thresholds",
            "Countries analyzed, No.": "150",
            "Definition": "Q80 growth; Q20 MV-unit density using 2050 cases",
            "Countries meeting both, No.": str(decoupled["strict_screen_positive_count"]),
            "Primary countries retained": f"{decoupled['strict_primary_retained']}/23",
            "Classification changes or uncertainty": f"Not retained: {', '.join(decoupled['strict_not_retained'])}",
        },
        {
            "Analysis": "Cross-specification support",
            "Countries analyzed, No.": "23 primary countries",
            "Definition": "Seven alternative burden, denominator, horizon, and resource definitions",
            "Countries meeting both, No.": "Not applicable",
            "Primary countries retained": "16/23 in all 7; 20/23 in at least 6",
            "Classification changes or uncertainty": "Support ranged from 4 to 7 specifications",
        },
        {
            "Analysis": "DIRAC records dated 2023 or later",
            "Countries analyzed, No.": str(recency["recent_record_countries"]),
            "Definition": "Q75/Q25 thresholds recalculated among recent records",
            "Countries meeting both, No.": str(recency["recent_screen_positive_count"]),
            "Primary countries retained": f"{recency['eligible_primary_retained_count']}/{recency['primary_recent_record_count']} eligible",
            "Classification changes or uncertainty": f"Seven primary records were older; added: {', '.join(recency['new_recent_screen_positive_iso3'])}",
        },
    ]
    pd.DataFrame(stability_rows).to_csv(STABILITY_TABLE, index=False, encoding="utf-8-sig")

    report = {
        "status": "pass",
        "input": str(DATASET.relative_to(ROOT)),
        "primary_metrics": primary_metrics,
        "outputs": [
            str(STABILITY_TABLE.relative_to(ROOT)),
            str(PRIMARY_SUMMARY_TABLE.relative_to(ROOT)),
        ],
    }
    REPORT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
